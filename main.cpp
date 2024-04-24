#include <iostream>
#include <zran.hpp>
#include "kseq++/seqio.hpp"
using namespace std;
using namespace klibpp;


int main(int argc, char **argv) {

    if (argc <= 1) {
        fprintf(stderr, "Command line arguments not provided\n");
    }
    if (strcmp(argv[1], "build") == 0) {
        // build mode
        build_index(argv[2], SPAN);
    } else {
        // use mode
        off_t offset = -1;
        off_t read_len = -1;
        if (argc > 2) {
            char *end;
            offset = strtoll(argv[4], &end, 0);
            if (*end) {
                fprintf(stderr, "zran: invalid offset\n");
                return 1;
            }

            char *end1;
            read_len = strtoll(argv[5], &end1, 0);
            if (*end1) {
                fprintf(stderr, "zran: invalid read_len\n");
                return 1;
            }
        }
        fprintf(stderr, "zran: extracting %d bytes at offset %ld\n", read_len, offset);
        read_index(argv[2], argv[3], offset, read_len);
    }
    return 0;
}