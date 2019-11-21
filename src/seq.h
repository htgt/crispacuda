#pragma once
#define ERROR_STR 0xFFFFFFFFFFFFFFFFull
void print_seq(uint64_t text, uint64_t length);
char *bits_to_string(uint64_t text, uint64_t length);
void populate_cmap();
uint64_t string_to_bits(const char *seq, uint64_t seq_length, uint64_t bits = 0);
