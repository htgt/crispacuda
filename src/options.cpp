#include "options.h"
#include "devices.h"
#include <cstring>

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
    if ( device >= 1 && set_device(device) != 0 ) {
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
    for ( uint64_t i = 0; i < num; i++ ) {
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

const struct search_t default_search = {NULL, NULL, false, NULL};
const struct options_t default_options = {true, 4};
