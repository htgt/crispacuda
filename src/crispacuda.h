#pragma once
#define MAX_CHAR_SIZE 30
#define VERSION 3
#include "options.h"
#include <ostream>

const int max_on_list = 2000;
const int max_off_list = 2000;

struct metadata_t {
    uint64_t num_seqs;
    uint64_t seq_length;
    uint64_t offset; //so we can give real mouse db ids
    uint8_t species_id;
    char species[MAX_CHAR_SIZE]; //use fixed char arrays so we don't have to store size
    char assembly[MAX_CHAR_SIZE];
};

struct targets_t {
    uint64_t off[max_off_list];
    uint64_t on[max_on_list];
    int offc;
    int onc;
};

struct userdata_t {
    metadata_t metadata;
    options_t options;
    uint64_t *h_crisprs;
    uint64_t *d_crisprs;
};

void do_find_off_targets(std::ostream &stream, userdata_t *userdata, crispr_t query); 
void do_search_by_seq(std::ostream &stream, userdata_t *userdata, crispr_t query, short pam_right); 
