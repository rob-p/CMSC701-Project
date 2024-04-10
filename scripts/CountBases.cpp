#include "utils/io.hpp"
#include <chrono>
#include <kseq++/seqio.hpp>
using namespace std;
using namespace klibpp;

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();
    vector<string> ids, genomes;  //helpers to read the fasta file
    KSeq record;
    SeqStreamIn iss(argv[1]);
    while (iss >> record) {
        ids.push_back(record.name);
        genomes.push_back(record.seq);
    }
    auto end_read = std::chrono::high_resolution_clock::now();
    cout << ids.size() << endl;
    cout << genomes.size() << endl;
    long long int countA = 0;
    long long int countC = 0;
    long long int countG = 0;
    long long int countT = 0;

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
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration in milliseconds
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_read - start);
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Output the duration
    std::cout << "Time taken (to read): " << duration1.count() << " milliseconds" << std::endl;
    std::cout << "Time taken (totak): " << duration2.count() << " milliseconds" << std::endl;

    return 0;
}