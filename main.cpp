#include <iostream>
#include <zran.hpp>
using namespace std;


// Demonstrate the use of deflate_index_build() and deflate_index_extract() by
// processing the file provided on the command line, and extracting LEN bytes
// from 2/3rds of the way through the uncompressed output, writing that to
// stdout. An offset can be provided as the second argument, in which case the
// data is extracted from there instead.
int main(int argc, char **argv) {
    // Open the input file.
    if (argc < 2 || argc > 4) {
        fprintf(stderr, "usage: zran file.raw [offset]\n");
        return 1;
    }
    FILE *in = fopen(argv[1], "rb");
    if (in == NULL) {
        fprintf(stderr, "zran: could not open %s for reading\n", argv[1]);
        return 1;
    }

    // Get optional offset.
    off_t offset = -1;
    if (argc > 2) {
        char *end;
        offset = strtoll(argv[2], &end, 0);
        if (*end) {
            fprintf(stderr, "zran: invalid offset %s\n", argv[2]);
            fclose(in);
            return 1;
        }
    }
    fprintf(stderr, "zran: extracting %d bytes at offset %ld\n", LEN, offset);

    // Either build an index or read one from a file
    struct deflate_index *index = NULL;
    int len;

    // Read index from file if provided
    if (argc == 4) {
        fprintf(stderr, "zran: attempting to read index from %s\n", argv[3]);

        FILE *index_file = fopen(argv[3], "rb");
        if (index_file == NULL) {
            fprintf(stderr, "zran: could not open %s for reading\n", argv[3]);
            fclose(in);
            return 1;
        }

        len = deflate_index_load(index_file, &index);
        fclose(index_file);
    } else {
        // Build index.
        len = deflate_index_build(in, SPAN, &index);
    }

    if (len < 0) {
        fclose(in);
        switch (len) {
            case Z_MEM_ERROR:
                fprintf(stderr, "zran: out of memory\n");
                break;
            case Z_BUF_ERROR:
                fprintf(stderr, "zran: %s ended prematurely\n", argv[1]);
                break;
            case Z_DATA_ERROR:
                fprintf(stderr, "zran: compressed data error in %s\n", argv[1]);
                break;
            case Z_ERRNO:
                fprintf(stderr, "zran: read error on %s\n", argv[1]);
                break;
            default:
                fprintf(stderr, "zran: error %d while building index\n", len);
        }
        return 1;
    }

    if (argc == 4) {
        fprintf(stderr, "zran: read index with %d access points!\n", len);
//        print_index(index);
    } else {
        fprintf(stderr, "zran: built index with %d access points!\n", len);
        print_index(index);

        // Save index to file
        char *filename = (char *) malloc(strlen(argv[1]) + 6);
        if (filename == NULL) {
            fprintf(stderr, "zran: out of memory\n");
            deflate_index_free(index);
            fclose(in);
            return 1;
        }
        strcpy(filename, argv[1]);
        strcat(filename, ".index");

        fprintf(stderr, "zran: attempting to write index to %s\n", filename);

        // Open the index file for writing.
        FILE *idx = fopen(filename, "wb");
        if (idx == NULL) {
            fprintf(stderr, "zran: could not open %s for writing\n", filename);
            deflate_index_free(index);
            fclose(in);
            return 1;
        }

        // Write the index to the file.
        len = deflate_index_save(idx, index);
        if (len != 0) {
            fclose(idx);
            fprintf(stderr, "zran: write error on %s\n", filename);
            deflate_index_free(index);
            fclose(in);
            return 1;
        }
        fprintf(stderr, "zran: wrote index with %d access points to %s\n", index->have, filename);

        // Clean up and exit
        fclose(idx);
        free(filename);
        deflate_index_free(index);
        fclose(in);
        return 0;
    }

    // Use index by reading some bytes from an arbitrary offset.
    unsigned char buf[LEN];
    if (offset == -1)
        offset = ((index->length + 1) << 1) / 3;
    ptrdiff_t got = deflate_index_extract(in, index, offset, buf, LEN);
    if (got < 0)
        fprintf(stderr, "zran: extraction failed: %s error\n",
                got == Z_MEM_ERROR ? "out of memory" : "input corrupted");
    else {
        fwrite(buf, 1, got, stdout);
//        fprintf(stderr, "\nzran: extracted %ld bytes at %ld using index read from disk\n", got, offset);
    }

    // Clean up and exit.
    deflate_index_free(index);
    fclose(in);
    return 0;
}