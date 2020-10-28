#pragma once
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>

struct crispr_t {
	uint64_t id;
	uint64_t seq;
	uint64_t rev_seq;
};

struct search_t {
	crispr_t *queries;
	char *index_file;
	bool search_by_seq;
	char *port;
};
extern const struct search_t default_search;

struct options_t {
	bool store_offs;
	int max_mismatches;
};
extern const struct options_t default_options;

int64_t read_options(int argc, char *argv[], search_t *search, options_t *options);
