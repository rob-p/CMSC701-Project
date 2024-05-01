#include <iostream>
#include <zran.hpp>
#include "kseq++/seqio.hpp"
#include "kseqcharstream.hpp"
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
        unsigned char* buf;
        int got;
        std::tie(buf, got) = read_index(argv[2], argv[3], record_idx, num_records);
        if (got < 0)
            fprintf(stderr, "zran: extraction failed: %s error\n",
                    got == Z_MEM_ERROR ? "out of memory" : "input corrupted");
        else {
            fwrite(buf, 1, got, stdout);
        }
        KseqCharStreamIn in(reinterpret_cast<const char*>(buf), got);

        klibpp::KSeq rec;
        while (in >> rec) {
          cout << rec.seq << endl;
        }
    }
    return 0;
}