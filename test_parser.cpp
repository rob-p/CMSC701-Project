#include "parser.hpp"
#include "indicators.hpp"
#include "CLI11.hpp"
#include <iostream>
#include <thread>
#include <vector>
#include <chrono>
using namespace std;

struct Bases {
  uint64_t A, C, G, T;
};

int main(int argc, char* argv[]) {
  CLI::App app{"test program for ffparser"};
  argv = app.ensure_utf8(argv);
 
  std::string fastqFile;
  std::string indexFile;
  size_t nt{4};
  app.add_option<std::string>("fastq-path", fastqFile, "path to input fastq file.")->required();
  app.add_option<std::string>("fastq-index-path", indexFile, "path to input fastq file index.")->required();
  app.add_option<size_t>("num-threads", nt, "number of parsing threads to use.")->required();
  CLI11_PARSE(app, argc, argv);

  ParrFQParser parser;
  parser.init(fastqFile, indexFile, nt);

  auto start = std::chrono::high_resolution_clock::now();
  cout << "Starting parsing" << endl;
  parser.start();

  cout << "Parsers Started" << endl;

  using namespace indicators;

  // Hide cursor
  show_console_cursor(false);

  indicators::ProgressBar bar{
    option::BarWidth{50},
    option::Start{" ["},
    option::Fill{"█"},
    option::Lead{"█"},
    option::Remainder{"-"},
    option::End{"]"},
    option::MaxProgress{parser.get_num_chunks()},
    option::ForegroundColor{Color::yellow},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };


  std::vector<std::thread> readers;
  std::vector<Bases> counters(nt, {0, 0, 0, 0});
  std::atomic<size_t> ctr{0};
  for (size_t i = 0; i < nt; ++i) {
    readers.emplace_back([&, i]() {
      auto rg = parser.get_read_chunk();
      size_t lctr{0};
      size_t pctr{0};
      klibpp::KSeq seq;
      uint64_t cur_rec{0};
      while (parser.refill(rg)) {
        bar.tick();
        auto& seq_stream = rg.get_seq_stream();
        while (seq_stream >> seq) { 
            //std::cerr << "rec : " << j << " / " << expected_rec << "\n";
            ++cur_rec;
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
        }
        if (cur_rec != rg.expected_rec) {
          std::cerr << "CHUNK NUM : " << rg.chunk_num << ", observed rec : " << cur_rec << ", expected rec : " << rg.expected_rec << "\n";
          std::cerr << "last_name: " << seq.name << ", seq : " << seq.seq << "\n";
        }

        ctr += cur_rec; 
        cur_rec = 0;

      }
    });
  }

  for (auto& t : readers) {
    t.join();
  }
  bar.mark_as_completed();
  parser.stop();
  // Show cursor
  indicators::show_console_cursor(true);

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
