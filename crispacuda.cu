#include <fstream>
#include <getopt.h>
#include <inttypes.h>
#include <iostream>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include "crispacuda.h"

char* bits_to_string(uint64_t text, uint64_t length) {
    char *s = (char*)malloc(length + 1);
    memset(s, 0, length + 1);
    uint64_t shift = 2 * ( length - 1 ); //there are twice as many bits as there are characters

    //fill with N if its an error string (all bits set to 1)
    if ( text == ERROR_STR ) {
        memset(s, 'N', length);
    }

    //extract each character from the text
    for ( int i = 0; i < length; i++, shift -= 2 ) {
        //put the character we're interested in at the very end
        //of the integer, and switch all remaining bits to 0 with & 0x3
        uint8_t character = (text >> shift) & 0x3;
        switch ( character ) {
            case 0: s[i] = 'A'; break;
            case 1: s[i] = 'C'; break;
            case 2: s[i] = 'G'; break;
            case 3: s[i] = 'T'; break;
            default: break;
        }
    }

    return s;
}

void print_seq(uint64_t text, uint64_t length) {
    char *seq = bits_to_string(text, length);
    printf("%s\n", seq);
    free(seq);
}

uint64_t revcom(uint64_t text, int length) {
    unsigned int num_bits = sizeof(text) * CHAR_BIT;

    uint64_t mask = 0xFFFFFFFFFFFFFFFFull >> ( (num_bits - (length * 2)) - 1 );
    text = ~text & mask;
    uint64_t reversed = text >> (length * 2);
    int shift = 0;

    for ( int i = 0; i < length; i++, shift += 2 ) {
        reversed <<= 2;
        reversed |= ( text >> shift ) & 0x3;
    }
    return reversed;
}

__device__ targets_t targets;

__device__ void push_back(uint64_t id, int mm) {
    int insert_pt = atomicAdd(&targets.offc, 1);
    if ( insert_pt < max_off_list ) {
        targets.off[insert_pt] = id;
    }
    if ( mm == 0 ) {
        insert_pt = atomicAdd(&targets.onc, 1);
        if ( insert_pt < max_on_list ) {
            targets.on[insert_pt] = id;
        }
    }
}

__global__
void find_off_targets(uint64_t *crisprs, crispr_t query, int *summary, metadata_t metadata) {
    const uint64_t pam_on = 0x1ull << metadata.seq_length *2;
    const uint64_t pam_off = ~pam_on;

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    if(blockIdx.x == 0 && threadIdx.x == 0) {
        targets.offc = 0;
        targets.onc  = 0;
    }
    __syncthreads();
    for ( uint64_t j = index; j < metadata.num_seqs; j+= stride ) {
        uint64_t test_crispr = crisprs[j];
        if ( test_crispr == ERROR_STR ) continue;

        uint64_t match = query.seq ^ test_crispr;
        if ( match & pam_on ) {
            match = query.rev_seq ^ test_crispr;
        }

        match = match & pam_off;
        match = (match | (match >> 1)) & (0x5555555555555555ull);
        int mm = __popcll(match);

        if ( mm <= max_mismatches ) {
            atomicAdd(&summary[mm], 1);
            push_back(j + 1 + metadata.offset, mm);
        }
    }
}

void write_output(FILE *fp, crispr_t query, targets_t targets, int *summary,
        int species_id, bool store_offs) {
    int onc = std::min(targets.onc, max_on_list);
    int offc = std::min(targets.offc, max_off_list);
    if ( fp != NULL ) {
        fwrite(&query.id, sizeof(uint64_t), 1, fp); 
        fwrite(summary, sizeof(int), max_mismatches+1, fp);
        fwrite(&onc, sizeof(int), 1, fp);
        fwrite(&offc, sizeof(int), 1, fp);
        fwrite(targets.on, sizeof(uint64_t), onc, fp);
        fwrite(targets.off, sizeof(uint64_t), offc, fp);
    } else {
        thrust::sort(thrust::host, targets.off, targets.off + offc, thrust::less<uint64_t>());
        std::cout << query.id << "\t" << int(species_id);
        const char *sep = seps[0];
        if(!store_offs || targets.offc > max_off_list) {
            std::cout << "\t\\N\t{";
        } else {
            std::cout << "\t{";
            for( int j = 0; j < offc; j++ ) {
                std::cout << sep << targets.off[j];
                sep = seps[1];
            }
            std::cout << "}\t{";
        }
        sep = seps[0];
        for( int j = 0; j <= max_mismatches; j++ ) {
            std::cout << sep << j << ": " << summary[j];
            sep = seps[2];
        }
        std::cout << "}" << std::endl;
    }
}

void calc_off_targets(FILE *fp, uint64_t *crisprs, crispr_t query, metadata_t metadata,
        options_t options) {
    int summary_size = (max_mismatches+1)*sizeof(int);
    int *summary;
    cudaMalloc((void**)&summary, summary_size);
    cudaMemset(summary, 0, summary_size);

    const int blockSize = 128;
    const int numBlocks = (metadata.num_seqs + blockSize - 1) / blockSize;
    find_off_targets<<<numBlocks, blockSize>>>(crisprs, query, summary, metadata); 
    cudaDeviceSynchronize();

    targets_t h_targets;
    int *h_summary = (int*)malloc(summary_size);
    cudaMemcpyFromSymbol(&h_targets, targets, sizeof(targets_t));
    cudaMemcpy(h_summary, summary, summary_size, cudaMemcpyDeviceToHost);
    
    write_output(fp, query, h_targets, h_summary, metadata.species_id, options.store_offs);

    cudaFree(summary);
    free(h_summary);
}

metadata_t load_metadata(FILE *fp) {
    uint8_t endian_test;
    fread(&endian_test, sizeof(uint8_t), 1, fp);
    if ( endian_test != 1 ) {
        fprintf(stderr, "Endianess of the file does not match your hardware\n");
        exit(1);
    }

    uint32_t file_version;
    fread(&file_version, sizeof(uint32_t), 1, fp);
    if ( file_version != VERSION ) {
        fprintf(stderr, "Index file is the wrong version! Please regenerate!\n");
        exit(1);
    }
    fprintf(stderr, "Version is %d\n", file_version);

    metadata_t metadata;
    fread(&metadata, sizeof(metadata_t), 1, fp);
    fprintf(stderr, "Assembly is %s (%s)\n", metadata.assembly, metadata.species);
    fprintf(stderr, "File has %" PRIu64 " sequences\n", metadata.num_seqs);
    fprintf(stderr, "Sequence length is %" PRIu64 "\n", metadata.seq_length);
    fprintf(stderr, "Offset is %" PRIu64 "\n", metadata.offset);
    fprintf(stderr, "Species id is %d\n", metadata.species_id);
    return metadata;
}

uint64_t read_options(int argc, char *argv[], search_t *search, options_t *options) {
    int c = -1;
    uint64_t start = 0, num = 0;
    while ( ( c = getopt(argc, argv, "s:n:i:o:qm:") ) != -1 )  {
        switch(c) {
            case 's': start = atol(optarg); break;
            case 'n': num   = atol(optarg); break;
            case 'i': (*search).index_file = optarg; break;
            case 'o': (*search).output_file = optarg; break;
            case 'q': (*options).store_offs = false; break;
            case 'm': (*options).max_mismatches = atoi(optarg); break;

        }
    }
    if ( (*search).index_file == NULL ) {
        fprintf(stderr, "An index file must be specified with the -i option\n");
        return 0;
    }
    if ( start != 0 && num == 0 ) {
        fprintf(stderr, "If -s is specified, -n must be also");
        return 0;
    }
    if ( start == 0 && num != 0 ) {
        fprintf(stderr, "If -n is specified, -s must be also");
        return 0;
    }
    uint64_t num_queries = argc - optind + num;
    (*search).queries = (crispr_t*)malloc(num_queries * sizeof(crispr_t));
    for ( int i = 0; i < num; i++ ) {
        (*search).queries[i].id = start + i;
    }
    for ( int i = optind; i < argc; i++ ) {
        uint64_t id = atol(argv[i]);
        if ( id == 0 ) {
            fprintf(stderr, "Could not parse '%s' an ID\n", argv[i]);
            return false;
        }
        (*search).queries[i - optind + num].id = id;
    }
    return num_queries;
}

int main(int argc, char *argv[]) {
    search_t search = default_search;
    options_t options = default_options;
    uint64_t num_queries = read_options(argc, argv, &search, &options);
    if ( num_queries == 0 ) {
        return 1;
    }

    FILE *fp = fopen(search.index_file, "r");
    if ( fp == NULL ) {
        fprintf(stderr, "Could not open index\n");
        exit(1);
    }
    metadata_t metadata = load_metadata(fp);
    clock_t t = clock();
    const uint64_t data_size = metadata.num_seqs * sizeof(uint64_t);
    uint64_t *h_crisprs = (uint64_t*)malloc(data_size);
    fread(h_crisprs, sizeof(uint64_t), metadata.num_seqs, fp);
    t = clock() - t;
    fclose(fp);
    fprintf(stderr, "Loading took %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

    for ( int i = 0; i < num_queries; i++ ) {
        if ( search.queries[i].id < metadata.offset + 1
                || search.queries[i].id > metadata.offset + metadata.num_seqs ) {
            fprintf(stderr, "%" PRIu64 " is not a valid ID in this index\n",
                    search.queries[i].id);
            return 2;
        }
        search.queries[i].seq = h_crisprs[search.queries[i].id - metadata.offset - 1];
        search.queries[i].rev_seq = revcom(search.queries[i].seq, metadata.seq_length);
    }

    uint64_t *crisprs;
    cudaMalloc((void**)&crisprs, data_size);
    cudaMemcpy(crisprs, h_crisprs, data_size, cudaMemcpyHostToDevice);
    free(h_crisprs);

    FILE *fw = NULL;
    if ( search.output_file != NULL ) {
        fw = fopen(search.output_file, "w");
        if ( fw == NULL ) {
            fprintf(stderr, "Could not open output file\n");
            exit(1);
        }
    }
    for ( int i = 0; i < num_queries; i++ ) {
        calc_off_targets(fw, crisprs, search.queries[i], metadata, options);
    }
    if ( fw != NULL ) {
        fclose(fw);
    }

    cudaFree(crisprs);
    free(search.queries);
}
