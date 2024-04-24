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
        off_t record_idx = -1;
        if (argc > 2) {
            char *end;
            record_idx = strtoll(argv[4], &end, 0);
            if (*end) {
                fprintf(stderr, "zran: invalid record_idx\n");
                return 1;
            }
        }
        read_index(argv[2], argv[3], record_idx);
    }
    return 0;
}