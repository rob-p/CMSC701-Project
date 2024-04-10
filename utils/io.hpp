#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

// My own implementation of FASTA reader to sanity check the reads
void ReadFastaFile(string file_path, vector<string> &ids, vector<string> &genomes) {
    // Read each line of the FASTA file in the variable input_line
    // identify id and genome and save in corresponding vectors.
    // This is a generic implementation to support both FASTA file with reference genome and query FAST file
    ifstream input_stream(file_path);
    string id = "", genome = "", input_line;

    while(getline(input_stream, input_line)) {
        // to take care of end of file blank line
        if (int(input_line.size()) == 0) {
            continue;
        }
        // id always starts with `>`
        if (input_line[0] == '>') {
            // got a single-line id. save it and the genome.
            // first saved genome will be empty. will have to delete it later.
            id = input_line.substr(1);
            ids.push_back(id);
            genomes.push_back(genome);
            genome = "";
        } else {
            genome += input_line;
        }
    }
    if (!genome.empty()) {
        // save the last genome
        genomes.push_back(genome);
    }
    // first entry will be empty because of the way I have written while loop above
    if (genomes[0] == "") {
        genomes.erase(genomes.begin());
    }
    input_stream.close();
}