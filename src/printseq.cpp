#include <cstdlib>
#include <cstring>
#include <inttypes.h>
#include <stdio.h>
#include "printseq.h"
#define ERROR_STR 0xFFFFFFFFFFFFFFFFull

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
