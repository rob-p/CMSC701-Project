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
        off_t span = 1048576L;
        char *end2;
        span = strtoll(argv[3], &end2, 0);
        if (*end2) {
            fprintf(stderr, "zran: invalid record_idx\n");
            return 1;
        }
        build_index(argv[2], span);

    } else {
        // use mode
        off_t record_idx = -1;
        off_t num_records = -1;
        if (argc > 2) {
            char *end;
            record_idx = strtoll(argv[4], &end, 0);
            if (*end) {
                fprintf(stderr, "zran: invalid record_idx\n");
                return 1;
            }

            char *end1;
            num_records = strtoll(argv[5], &end1, 0);
            if (*end1) {
                fprintf(stderr, "zran: invalid num_records\n");
                return 1;
            }
        }
        read_index(argv[2], argv[3], record_idx, num_records);
    }
    return 0;
}