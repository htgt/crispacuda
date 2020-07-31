#include <fstream>
#include <getopt.h>
#include <inttypes.h>
#include <iostream>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include <vector>
#include "crispacuda.h"
#include "devices.h"
#include "mongoose.h"
#include "seq.h"

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

uint64_t *crisprs;
uint64_t *h_crisprs;
metadata_t metadata;
options_t options;
__device__ targets_t targets;
__constant__ __device__ uint64_t pam_on;
__constant__ __device__ uint64_t pam_off;
__constant__ __device__ metadata_t d_metadata;
__constant__ __device__ options_t d_options;

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
void find_off_targets(uint64_t *crisprs, crispr_t query, int *summary) {
    int tid = threadIdx.x;
    int index = blockIdx.x * blockDim.x + tid;
    int stride = blockDim.x * gridDim.x;
    if(blockIdx.x == 0 && tid == 0 ) {
        targets.offc = 0;
        targets.onc  = 0;
    }
    __syncthreads();
    for ( uint64_t j = index; j < d_metadata.num_seqs; j+= stride ) {
        uint64_t test_crispr = crisprs[j];
        if ( test_crispr == ERROR_STR ) continue;

        uint64_t match = query.seq ^ test_crispr;
        if ( match & pam_on ) {
            match = query.rev_seq ^ test_crispr;
        }

        match = match & pam_off;
        match = (match | (match >> 1)) & (0x5555555555555555ull);
        int mm = __popcll(match);

        if ( mm <= d_options.max_mismatches ) {
            atomicAdd(&summary[mm], 1);
            push_back(j + 1 + d_metadata.offset, mm);
        }
    }
}

void write_output(std::ostream &stream, crispr_t query, int *summary, targets_t targets) {
    int onc = std::min(targets.onc, max_on_list);
    int offc = std::min(targets.offc, max_off_list);
    thrust::sort(thrust::host, targets.off, targets.off + offc, thrust::less<uint64_t>());
    stream << query.id << "\t" << int(metadata.species_id);
    const char *sep = seps[0];
    if(!options.store_offs || targets.offc > max_off_list) {
        stream << "\t\\N\t{";
    } else {
        stream << "\t{";
        for( int j = 0; j < offc; j++ ) {
            stream << sep << targets.off[j];
            sep = seps[1];
        }
        stream << "}\t{";
    }
    sep = seps[0];
    for( int j = 0; j <= options.max_mismatches; j++ ) {
        stream << sep << j << ": " << summary[j];
        sep = seps[2];
    }
    stream << "}" << std::endl;
}

void calc_off_targets(std::ostream &stream, crispr_t query) {
    int summary_size = (options.max_mismatches + 1) * sizeof(int);
    int *summary;
    cudaMalloc((void**)&summary, summary_size);
    cudaMemset(summary, 0, options.max_mismatches);
    const int blockSize = 128;
    const int numBlocks = (metadata.num_seqs + blockSize - 1) / blockSize;
    find_off_targets<<<numBlocks, blockSize>>>(crisprs, query, summary); 
    cudaDeviceSynchronize();

    targets_t h_targets;
    int *h_summary = (int*)malloc(summary_size);
    cudaMemcpyFromSymbol(&h_targets, targets, sizeof(targets_t), 0, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_summary, summary, summary_size, cudaMemcpyDeviceToHost);
    
    write_output(stream, query, h_summary, h_targets);

    free(h_summary);
}

metadata_t load_metadata(FILE *fp) {
    uint8_t endian_test;
    CHECK_FREAD(&endian_test, sizeof(uint8_t), 1, fp);
    if ( endian_test != 1 ) {
        fprintf(stderr, "Endianess of the file does not match your hardware\n");
        exit(1);
    }

    uint32_t file_version;
    CHECK_FREAD(&file_version, sizeof(uint32_t), 1, fp);
    if ( file_version != VERSION ) {
        fprintf(stderr, "Index file is the wrong version! Please regenerate!\n");
        exit(1);
    }
    fprintf(stderr, "Version is %d\n", file_version);

    metadata_t metadata;
    CHECK_FREAD(&metadata, sizeof(metadata_t), 1, fp);
    fprintf(stderr, "Assembly is %s (%s)\n", metadata.assembly, metadata.species);
    fprintf(stderr, "File has %" PRIu64 " sequences\n", metadata.num_seqs);
    fprintf(stderr, "Sequence length is %" PRIu64 "\n", metadata.seq_length);
    fprintf(stderr, "Offset is %" PRIu64 "\n", metadata.offset);
    fprintf(stderr, "Species id is %d\n", metadata.species_id);
    return metadata;
}

int64_t read_options(int argc, char *argv[], search_t *search, options_t *options) {
    int c = -1, device = 0;
    uint64_t start = 0, num = 0;
    bool show_help = false;
    while ( ( c = getopt(argc, argv, "s:n:i:m:d:p:hqz") ) != -1 )  {
        switch(c) {
            case 's': start = atol(optarg); break;
            case 'n': num   = atol(optarg); break;
            case 'i': (*search).index_file = optarg; break;
            case 'q': (*options).store_offs = false; break;
            case 'm': (*options).max_mismatches = atoi(optarg); break;
            case 'h': show_help = true; break;
            case 'd': device = atoi(optarg); break;
            case 'z': (*search).search_by_seq = true; break;
            case 'p': (*search).port = optarg; break;
        }
    }
    if ( show_help ) {
        printf("CRISPACUDA\n");
        printf("Searches for CRISPR off-targets on a GPU\n");
        printf("Usage: crispacuda [options] <ids...>\n");
        printf("Contact: Joel Rein joel.rein@sanger.ac.uk\n\n");
        printf("OPTIONS:\n");
        printf("-i <FILE>\n");
        printf("\tSpecify the index file to search. REQUIRED.\n");
        printf("-s <INT>\n");
        printf("\tThe index to start searches from\n");
        printf("-n <INT>\n");
        printf("\tHow many off-targets to calculate from start.\n");
        printf("\tRequired if, and only if, -s is specified.\n");
        printf("-m <INT>\n");
        printf("\tThe maximum number of mismatches to record results for. Defaults to 4.\n");
        printf("-q\n\tDo not report a list of off-targets.\n");
        printf("\t Has no effect if outputting to a binary file.\n");
        printf("-h\n\tPrint this help message and GPU information\n");
        printf("-d <INT>\n\tThe GPU device to use.\n");
        printf("-z\n\tSpecify CRISPRs as strings rather than IDs\n");
        printf("-p\n\tStart a webserver on port <number>\n");
        printf("Following these arguments you may specify individual CRISPRs to search.\n\n");
        printf("EXIT CODES\n");
        printf("\tNegative exit codes indicate errors in command line options.\n");
        printf("\tPositive exit codes indicate errors running the search.\n");
        printf("\tReturns 0 on success.\n\n");
        show_devices();
        return 0;
    }
    if ( device >= 1 && cudaSetDevice(device) != cudaSuccess ) {
        fprintf(stderr, "Could not use device %d\n", device);
        return -1;
    }
    if ( (*search).index_file == NULL ) {
        fprintf(stderr, "An index file must be specified with the -i option\n");
        return -2;
    }
    if ( start != 0 && num == 0 ) {
        fprintf(stderr, "If -s is specified, -n must be also");
        return -3;
    }
    if ( start == 0 && num != 0 ) {
        fprintf(stderr, "If -n is specified, -s must be also");
        return -4;
    }
    if ( start != 0 && (*search).search_by_seq ) {
        fprintf(stderr, "-s and -n are not compatible with search by sequence\n");
        return -5;
    }
    int64_t num_queries = argc - optind + num;
    (*search).queries = (crispr_t*)malloc(num_queries * sizeof(crispr_t));
    memset((*search).queries, 0, num_queries * sizeof(crispr_t));
    for ( int i = 0; i < num; i++ ) {
        (*search).queries[i].id = start + i;
    }
    for ( int i = optind; i < argc; i++ ) {
        uint64_t id = 0;
        if ((*search).search_by_seq) {
            id = i;
        } else {
            id = atol(argv[i]);
            if ( id == 0 ) {
                fprintf(stderr, "Could not parse '%s' an ID\n", argv[i]);
                return -6;
            }
        }
        (*search).queries[i - optind + num].id = id;
    }
    return num_queries;
}

void parse_request(struct mg_str req, std::vector<uint64_t> &ids) {
    char *query, *id;
    query = (char*)malloc(req.len + 1);
    memset(query, 0, req.len + 1);
    strncpy(query, req.p, req.len);
    while ( (id = strsep(&query, "\n")) != NULL ) {
        ids.push_back(atol(id));
    }
}

static void handle_request(struct mg_connection *c, int ev, void *p) {
    if ( ev == MG_EV_HTTP_REQUEST ) {
        struct http_message *hm = (struct http_message *)p;
        
        if( mg_vcmp(&hm->uri, "/search") == 0 ) {
            std::vector<uint64_t> ids;
            parse_request(hm->body, ids);
            std::stringstream results;
            for(uint64_t id : ids) {
                crispr_t crispr;
                if ( id <= metadata.offset || id > metadata.offset + metadata.num_seqs ) {
                    continue;
                }
                crispr.id = id;
                crispr.seq = h_crisprs[id - metadata.offset - 1];
                crispr.rev_seq = revcom(crispr.seq, metadata.seq_length);
                calc_off_targets(results, crispr);
            }
            const std::string tmp = results.str();
            struct mg_str response = mg_mk_str(tmp.c_str());
            mg_send_head(c, 200, response.len, "Content-Type: text/plain");
            mg_printf(c, "%.*s", (int)response.len, response.p);
        } else if( mg_vcmp(&hm->uri, "/favicon.ico") == 0 ) {
            mg_http_serve_file(c, hm, "favicon.ico", mg_mk_str("image/ico"), mg_mk_str(""));
        } else {
            mg_http_serve_file(c, hm, "index.htm", mg_mk_str("text/html"), mg_mk_str(""));
        }
    }
}

void run_server(char *port) {
    struct mg_mgr mgr;
    struct mg_connection *c;
    mg_mgr_init(&mgr, NULL);
    c = mg_bind(&mgr, port, handle_request);
    mg_set_protocol_http_websocket(c);
    for(;;) {
        mg_mgr_poll(&mgr, 1000000);
    }
    mg_mgr_free(&mgr);
}

int main(int argc, char *argv[]) {
    populate_cmap();
    search_t search = default_search;
    options = default_options;
    int64_t num_queries = read_options(argc, argv, &search, &options);
    if ( num_queries <= 0  && search.port == NULL) {
        return num_queries;
    }

    FILE *fp = fopen(search.index_file, "r");
    if ( fp == NULL ) {
        fprintf(stderr, "Could not open index\n");
        exit(1);
    }
    metadata = load_metadata(fp);
    clock_t t = clock();
    const uint64_t data_size = metadata.num_seqs * sizeof(uint64_t);
    h_crisprs = (uint64_t*)malloc(data_size);
    CHECK_FREAD(h_crisprs, sizeof(uint64_t), metadata.num_seqs, fp);
    t = clock() - t;
    fclose(fp);
    fprintf(stderr, "Loading took %f seconds\n", ((float)t)/CLOCKS_PER_SEC);
    
    uint64_t h_pam_on = 0x1ull << metadata.seq_length *2;
    uint64_t h_pam_off = ~h_pam_on;
    CHECK_CUDA(cudaMemcpyToSymbol(pam_on, &h_pam_on, sizeof(uint64_t)));
    CHECK_CUDA(cudaMemcpyToSymbol(pam_off, &h_pam_off, sizeof(uint64_t)));
    CHECK_CUDA(cudaMemcpyToSymbol(d_metadata, &metadata, sizeof(metadata_t)));
    CHECK_CUDA(cudaMemcpyToSymbol(d_options, &options, sizeof(options_t)));

    for ( int i = 0; i < num_queries; i++ ) {
        if ( search.search_by_seq ) {
            search.queries[i].seq = string_to_bits(argv[search.queries[i].id], metadata.seq_length, 1);
            if ( search.queries[i].seq == ERROR_STR ) {
                fprintf(stderr, "%s is not a valid sequence\n", argv[search.queries[i].id]);
                return 2;
            }
            search.queries[i].id = 0;
        } else {
            if ( search.queries[i].id < metadata.offset + 1
                    || search.queries[i].id > metadata.offset + metadata.num_seqs ) {
                fprintf(stderr, "%" PRIu64 " is not a valid ID in this index\n", search.queries[i].id);
                return 2;
            }
            search.queries[i].seq = h_crisprs[search.queries[i].id - metadata.offset - 1];
        }
        search.queries[i].rev_seq = revcom(search.queries[i].seq, metadata.seq_length);
    }

    size_t free_memory, total_memory;
    cudaMemGetInfo(&free_memory, &total_memory);
    fprintf(stderr, "Requires %" PRIu64 "mb of GPU memory, %" PRIu64 "mb is available\n",
            data_size >> 20, free_memory >> 20);
    if ( data_size > free_memory ) {
        fprintf(stderr, "Insufficient GPU memory, exiting.\n");
        return 3;
    }
    CHECK_CUDA(cudaMalloc((void**)&crisprs, data_size));
    CHECK_CUDA(cudaMemcpy(crisprs, h_crisprs, data_size, cudaMemcpyHostToDevice));

    for ( int i = 0; i < num_queries; i++ ) {
        calc_off_targets(std::cout, search.queries[i]);
    }
    if ( search.port != NULL ) {
        printf("Starting server on port %s...\n", search.port);
        run_server(search.port);
    }

    free(h_crisprs);
    cudaFree(crisprs);
    free(search.queries);
    return 0;
}
