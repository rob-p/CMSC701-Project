#include "utils/io.hpp"
using namespace std;

int main(int argc, char* argv[]) {
    vector<string> ids, genomes;  //helpers to read the fasta file
    ReadFastaFile(argv[1], ids, genomes);
    int countA = 0;
    int countC = 0;
    int countG = 0;
    int countT = 0;

    for(auto &s: genomes) {
        for(auto &ch: s) {
            if (ch == 'A') {
                countA++;
            } else if (ch == 'C') {
                countC++;
            } else if (ch == 'G') {
                countG++;
            } else if (ch == 'T') {
                countT++;
            }
        }
    }
    cout <<"#A: " << countA << endl;
    cout <<"#C: " << countC << endl;
    cout <<"#G: " << countG << endl;
    cout <<"#T: " << countT << endl;

    return 0;
}