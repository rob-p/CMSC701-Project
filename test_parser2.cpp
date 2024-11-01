#include "parser.hpp"
#include <iostream>
#include <thread>
#include <vector>
#include <chrono>
using namespace std;

struct Bases {
  uint64_t A, C, G, T;
};

int main(int argc, char* argv[]) {
  if (argc < 5) {
    std::cerr << "Command line arguments not provided\n";
    std::cerr << "Usage ./test_parser <fastq_file> <index_file> <num_threads>\n";
  }
  std::string fastqFile = argv[1];
  std::string indexFile = argv[2];
  size_t nt = stoi(argv[3]);  // number of consumer threads
  //size_t np = stoi(argv[4]);  // number of producer threads

  ParrFQParser parser;
  parser.init2(fastqFile, indexFile, nt);

  auto start = std::chrono::high_resolution_clock::now();
  cout << "Starting parsing" << endl;
  parser.start2();

  cout << "Parsers Started" << endl;

  std::vector<std::thread> readers;
  std::vector<Bases> counters(nt, {0, 0, 0, 0});
  std::atomic<size_t> ctr{0};
  for (size_t i = 0; i < nt; ++i) {
    readers.emplace_back([&, i]() {
      auto rg = parser.get_read_chunk();
      size_t lctr{0};
      size_t pctr{0};
      klibpp::KSeq seq;
      while (parser.refill(rg)) {
        auto& seq_stream = rg.get_seq_stream();
        while (seq_stream >> seq) {
            ++ctr; 
            for (size_t j = 0; j < seq.seq.length(); ++j) {
              char c = seq.seq[j];
              switch (c) {
              case 'A':
                counters[i].A++;
                break;
              case 'C':
                counters[i].C++;
                break;
              case 'G':
                counters[i].G++;
                break;
              case 'T':
                counters[i].T++;
                break;
              default:
                break;
              }
            }
          /*
          pctr = lctr;
          if (lctr > 1000000) {
              lctr = 0;
              pctr = 0;
              std::cout << "parsed " << ctr << " reads.\n";
          }
          */
        }
      }
    });
  }

  for (auto& t : readers) {
    t.join();
  }

  parser.stop2();

  Bases b = {0, 0, 0, 0};
  for (size_t i = 0; i < nt; ++i) {
    b.A += counters[i].A;
    b.C += counters[i].C;
    b.G += counters[i].G;
    b.T += counters[i].T;
  }
  std::cerr << "\n";
  std::cerr << "Parsed " << ctr << " total read pairs.\n";
  std::cerr << "\n#A = " << b.A << '\n';
  std::cerr << "#C = " << b.C << '\n';
  std::cerr << "#G = " << b.G << '\n';
  std::cerr << "#T = " << b.T << '\n';
  auto end = std::chrono::high_resolution_clock::now();

  // Calculate the duration in milliseconds
  auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  // Output the duration
  std::cout << "Time taken (total): " << duration2.count() << " milliseconds" << std::endl;
  return 0;
}
