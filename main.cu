#include <fstream>
#include <inttypes.h>

#define ERROR_STR 0xFFFFFFFFFFFFFFFFull
#define MAX_CHAR_SIZE 30
#define VERSION 3

struct metadata_t {
    uint64_t num_seqs;
    uint64_t seq_length;
    uint64_t offset; //so we can give real mouse db ids
    uint8_t species_id;
    char species[MAX_CHAR_SIZE]; //use fixed char arrays so we don't have to store size
    char assembly[MAX_CHAR_SIZE];
};

struct crispr_t {
    uint64_t id;
    uint64_t seq;
    uint64_t rev_seq;
};

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

__global__
void find_off_targets(uint64_t *crisprs, uint64_t fwd, uint64_t rev, int *summary,
        metadata_t metadata, int max_mismatches) {
    const uint64_t pam_on = 0x1ull << metadata.seq_length *2;
    const uint64_t pam_off = ~pam_on;

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for ( uint64_t j = index; j < metadata.num_seqs; j+= stride ) {
        uint64_t test_crispr = crisprs[j];
        if ( test_crispr == ERROR_STR ) continue;

        uint64_t match = fwd ^ test_crispr;
        if ( match & pam_on ) {
            match = rev ^ test_crispr;
        }

        match = match & pam_off;
        match = (match | (match >> 1)) & (0x5555555555555555ull);
        int mm = __popcll(match);

        if ( mm <= max_mismatches ) {
            atomicAdd(&summary[mm], 1);
        }
    }
}

void calcOffSeqs(uint64_t *crisprs, crispr_t query, metadata_t metadata, int max_mismatches) {
    int *summary;
    cudaMallocManaged(&summary, (max_mismatches+1)*sizeof(int));
    cudaMemset(summary, 0, sizeof(summary));
    const int blockSize = 128;
    const int numBlocks = (metadata.num_seqs + blockSize - 1) / blockSize;
    find_off_targets<<<numBlocks, blockSize>>>(crisprs, query.seq, query.rev_seq, summary,
            metadata, max_mismatches); 
    cudaDeviceSynchronize();

    printf("%10" PRIu64 " {0: %d; 1: %d; 2: %d; 3: %d 4: %d}\n",
            query.id, summary[0], summary[1], summary[2], summary[3], summary[4]);
    cudaFree(summary);
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
    printf("Version is %d\n", file_version);

    metadata_t metadata;
    fread(&metadata, sizeof(metadata_t), 1, fp);
    printf("Metadata size: %" PRIu64 "\n", sizeof(metadata_t));
    printf("Assembly is %s (%s)\n", metadata.assembly, metadata.species);
    printf("File has %" PRIu64 " sequences\n", metadata.num_seqs);
    printf("Sequence length is %" PRIu64 "\n", metadata.seq_length);
    printf("Offset is %" PRIu64 "\n", metadata.offset);
    printf("Species id is %d\n", metadata.species_id);
    return metadata;
}

int main(void) {
    FILE *fp = fopen("/nfs/team87/crispr_indexes/GRCh38_index.bin", "r");
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

    int num_queries = 1000;
    crispr_t *queries = (crispr_t*)malloc(num_queries * sizeof(crispr_t));
    for ( int i = 0; i < num_queries; i++ ) {
        queries[i].id = 1106710986 + i;
        queries[i].seq = h_crisprs[queries[i].id - metadata.offset - 1];
        queries[i].rev_seq = revcom(queries[i].seq, metadata.seq_length);
    }

    uint64_t *crisprs;
    cudaMallocManaged(&crisprs, data_size);
    cudaMemcpy(crisprs, h_crisprs, data_size, cudaMemcpyHostToDevice);
    free(h_crisprs);

    for ( int i = 0; i < num_queries; i++ ) {
        calcOffSeqs(crisprs, queries[i], metadata, 4);
    }

    free(queries);
    cudaFree(crisprs);
}
