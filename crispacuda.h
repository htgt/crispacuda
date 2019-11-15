#pragma once
#define ERROR_STR 0xFFFFFFFFFFFFFFFFull
#define MAX_CHAR_SIZE 30
#define VERSION 3

const int max_mismatches = 4;
const int max_on_list = 2000;
const int max_off_list = 2000;
const char *seps[2] = { "", ", " };

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

