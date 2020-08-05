#include <climits>
#include <cstdlib>
#include <cstring>
#include <inttypes.h>
#include <stdio.h>
#include "seq.h"
#define ERROR_STR 0xFFFFFFFFFFFFFFFFull

uint8_t cmap[256];

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

char* bits_to_string(uint64_t text, uint64_t length) {
    char *s = (char*)malloc(length + 1);
    memset(s, 0, length + 1);
    uint64_t shift = 2 * ( length - 1 ); //there are twice as many bits as there are characters

    //fill with N if its an error string (all bits set to 1)
    if ( text == ERROR_STR ) {
        memset(s, 'N', length);
    }

    //extract each character from the text
    for ( uint64_t i = 0; i < length; i++, shift -= 2 ) {
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

void populate_cmap() {
    memset(cmap, 4, 256);
    cmap['a'] = cmap['A'] = 0;
    cmap['c'] = cmap['C'] = 1;
    cmap['g'] = cmap['G'] = 2;
    cmap['t'] = cmap['T'] = 3;
}

uint64_t string_to_bits(const char *seq, uint64_t seq_length, uint64_t bits) {
	for ( uint64_t j = 0; j < seq_length; j++ ) {
		uint8_t const c = seq[j];

		if ( cmap[c] == 4 ) {
			bits = ERROR_STR;
			break;
		}
		else {
			bits <<= 2;
			bits |= cmap[c];
		}
	}

	return bits;
}
