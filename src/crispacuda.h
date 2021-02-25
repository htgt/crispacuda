#pragma once
#define MAX_CHAR_SIZE 30
#define VERSION 3

const int max_on_list = 2000;
const int max_off_list = 2000;
const char *seps[3] = { "", ",", ", " };

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

struct targets_t {
    uint64_t off[max_off_list];
    uint64_t on[max_on_list];
    int offc;
    int onc;
};

struct search_t {
    crispr_t *queries;
    char *index_file;
    char *output_file;
    bool search_by_seq;
} default_search = {NULL, NULL, NULL, false};

struct options_t {
    bool store_offs;
    int max_mismatches;
} default_options = {true, 4};

#ifndef CHECK_CUDA
static void checked_cuda(cudaError_t err, const char *file, int line) {
    if (err != cudaSuccess) {
        fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err), file, line );
        exit( EXIT_FAILURE );
    }
}
#define CHECK_CUDA( err ) (checked_cuda( err, __FILE__, __LINE__ ))
#endif

#ifndef CHECK_FREAD
static void checked_fread(void *ptr, size_t size, size_t count, FILE *fp,
        const char *file, int line) {
    size_t actual = fread(ptr, size, count, fp);
    if ( actual != count ) {
        fprintf(stderr, "Read %" PRIu64 " elements but expected %" PRIu64 " (%zu * %zu) in %s at line %d\n",
                actual, count, size, count, file, line);
        exit( EXIT_FAILURE );
    }
}
#define CHECK_FREAD(ptr,size,count,fp) (checked_fread(ptr, size, count, fp, __FILE__, __LINE__)) 
#endif
