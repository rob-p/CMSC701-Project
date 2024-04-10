#include "FastxParser.hpp"
#include <iostream>
#include <thread>
#include <vector>
#include <chrono>
using namespace std;

struct Bases {
  uint32_t A, C, G, T;
};

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

   std::vector<std::string> files;
   files.push_back(argv[1]);

   cout << files[0] << endl;

  size_t nt = stoi(argv[2]);  // number of consumer threads
  size_t np = stoi(argv[3]);  // number of producer threads
  fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(files, nt, np);
  cout << "Starting parsing" << endl;
  parser.start();
  cout << "Ending parsing" << endl;
  auto end_read = std::chrono::high_resolution_clock::now();

  std::vector<std::thread> readers;
  std::vector<Bases> counters(nt, {0, 0, 0, 0});
  std::atomic<size_t> ctr{0};
  for (size_t i = 0; i < nt; ++i) {
    readers.emplace_back([&, i]() {
      auto rg = parser.getReadGroup();
      size_t lctr{0};
      size_t pctr{0};
      while (true) {
        if (parser.refill(rg)) {
          auto chunk_frag_offset = rg.chunk_frag_offset();
          std::cerr << "chunk_offset_info: [file_idx: " << chunk_frag_offset.file_idx 
                    << ", frag_idx:" << chunk_frag_offset.frag_idx << ", chunk_size: " << rg.size() <<"]\n";
          for (auto& seq : rg) {
            ++lctr;

//            auto& seq = seqPair.first;
//            auto& seq2 = seqPair.second;

            size_t j = 0;
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
            //for (size_t j = 0; j < seq2.seq.length(); ++j) {
//              c = seq2.seq[j];
//              switch (c) {
//              case 'A':
//                counters[i].A++;
//                break;
//              case 'C':
//                counters[i].C++;
//                break;
//              case 'G':
//                counters[i].G++;
//                break;
//              case 'T':
//                counters[i].T++;
//                break;
//              default:
//                break;
//              }
           // }
          }
          ctr += (lctr - pctr);
          pctr = lctr;
          if (lctr > 1000000) {
              lctr = 0;
              pctr = 0;
              //std::cout << "parsed " << ctr << " read pairs.\n";
          }

        } else {
          break;
        }
      }
    });
  }

  for (auto& t : readers) {
    t.join();
  }

  parser.stop();

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
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_read - start);
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Output the duration
    std::cout << "Time taken (to read): " << duration1.count() << " milliseconds" << std::endl;
    std::cout << "Time taken (totak): " << duration2.count() << " milliseconds" << std::endl;
  return 0;
}
