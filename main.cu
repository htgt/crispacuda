#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdint>
#include <inttypes.h>

#define OFFSET    900000001
#define NUM_SEQS  303795341
#define MAX_MISMATCHES 4
#define ERROR_STR 0xFFFFFFFFFFFFFFFFull
#define SEQ_LEN 20

char* bits_to_string(uint64_t text) {
    char *s = (char*)malloc(SEQ_LEN + 1);
    memset(s, 0, SEQ_LEN + 1);
    uint64_t shift = 2 * ( SEQ_LEN - 1 ); //there are twice as many bits as there are characters

	//fill with N if its an error string (all bits set to 1)
	if ( text == ERROR_STR ) {
        memset(s, 'N', SEQ_LEN);
	}

	//extract each character from the text
	for ( int i = 0; i < SEQ_LEN; i++, shift -= 2 ) {
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

void print_seq(uint64_t text) {
    char *seq = bits_to_string(text);
    printf("%s\n", seq);
    free(seq);
}

uint64_t revcom(uint64_t text) {
	unsigned int num_bits = sizeof(text) * CHAR_BIT;
	
    uint64_t mask = 0xFFFFFFFFFFFFFFFFull >> ( (num_bits - (SEQ_LEN * 2)) - 1 );
	text = ~text & mask;
	uint64_t reversed = text >> (SEQ_LEN * 2);
	int shift = 0;

	for ( int i = 0; i < SEQ_LEN; i++, shift += 2 ) {
		reversed <<= 2;
		reversed |= ( text >> shift ) & 0x3;
	}
	return reversed;
}

__global__
void find_off_targets(uint64_t *crisprs, uint64_t fwd, uint64_t rev, int *summary) {
    const uint64_t pam_on = 0x1ull << SEQ_LEN*2;
    const uint64_t pam_off = ~pam_on;

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for ( uint64_t j = index; j < NUM_SEQS; j+= stride ) {
        uint64_t test_crispr = crisprs[j];
        if ( test_crispr == ERROR_STR ) continue;
        
        uint64_t match = fwd ^ test_crispr;
        if ( match & pam_on ) {
            match = rev ^ test_crispr;
        }

        match = match & pam_off;
        match = (match | (match >> 1)) & (0x5555555555555555ull);
        int mm = __popcll(match);

        if ( mm <= MAX_MISMATCHES ) {
            atomicAdd(&summary[mm], 1);
        }
    }
}

void calcOffSeqs(uint64_t *crisprs, uint64_t *h_crisprs, uint64_t id) {
    int *summary;
    cudaMallocManaged(&summary, (MAX_MISMATCHES+1)*sizeof(int));
    cudaMemset(summary, 0, sizeof(summary));
    uint64_t fwd = h_crisprs[id - OFFSET];
    uint64_t rev = revcom(fwd);
    const int blockSize = 256;
    const int numBlocks = (NUM_SEQS + blockSize - 1) /blockSize;
    find_off_targets<<<numBlocks, blockSize>>>(crisprs, fwd, rev, summary); 
    cudaDeviceSynchronize();

    printf("%10" PRIu64 " {0: %d; 1: %d; 2: %d; 3: %d 4: %d}\n",
            id, summary[0], summary[1], summary[2], summary[3], summary[4]);

    cudaFree(summary);
}

int main(void) {
    const uint64_t data_size = NUM_SEQS * sizeof(uint64_t);
    FILE *fp = fopen("/nfs/team87/crispr_indexes/GRCh38_index.bin", "r");
    if ( fp == NULL ) {
        fprintf(stderr, "Could not open index\n");
        exit(1);
    }
    fseek(fp, 93, SEEK_SET);
    clock_t t = clock();
    uint64_t *seq_data = (uint64_t*)malloc(data_size);
    fread(seq_data, sizeof(uint64_t), NUM_SEQS, fp);
    t = clock() - t;
    fclose(fp);
    fprintf(stderr, "Loading took %f seconds\n", ((float)t)/CLOCKS_PER_SEC);
    
    uint64_t *crisprs;
    cudaMallocManaged(&crisprs, data_size);
    cudaMemcpy(crisprs, seq_data, data_size, cudaMemcpyHostToDevice);

    for ( int i = 0; i < 1000; i++ ) {
        calcOffSeqs(crisprs, seq_data, 1106710986 + i);
    }

    cudaFree(crisprs);
    free(seq_data);
}
