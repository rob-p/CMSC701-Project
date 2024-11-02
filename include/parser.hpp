#pragma once
#include "concurrentqueue/concurrentqueue.h"
#include "kseq++/seqio.hpp"
#include "kseqcharstream.hpp"
#include "zran.hpp"
#include <atomic>
#include <functional>
#include <memory>
#include <stdio.h>
#include <thread>
#include <vector>

class ReadChunk {
public:
  ReadChunk(std::string fastq_file_name, struct deflate_index* idx) {
    FILE* fp = fopen(fastq_file_name.c_str(), "rb");
    if (fp != nullptr) {
      file_ptr_ = std::unique_ptr<FILE, std::function<void(FILE*)>>(fp, [](FILE* f) {fclose(f);});
    }
    idx_ptr_ = new struct deflate_index(*idx);
  }

  KseqCharStreamIn& get_seq_stream() { return *in_stream_; }

  std::unique_ptr<FILE, std::function<void(FILE*)>> file_ptr_{nullptr}; // this reader's local file ptr
  struct deflate_index* idx_ptr_{nullptr}; // this reader's local file ptr
  std::unique_ptr<KseqCharStreamIn> in_stream_{nullptr};
  std::vector<unsigned char> data_;
  // the number of records we expect to parse from
  // this chunk.
  uint64_t expected_rec;
  // the underlying assigned chunk num
  uint64_t chunk_num;
};

class ParrFQParser {
public:
  ParrFQParser()
      : m_index(nullptr,
                [](struct deflate_index *p) { deflate_index_free(p); }) {}

  ~ParrFQParser();

  int init(const std::string &fastqFilename, const std::string &indexFileName,
           uint64_t num_consumers);


  // Start and stop the parser
  int start();
  int stop();
  uint64_t get_num_chunks();

  // Consumer functions
  ReadChunk get_read_chunk();
  bool refill(ReadChunk& tlc);

private:
  std::unique_ptr<struct deflate_index,
                  std::function<void(struct deflate_index *)>>
      m_index;

  uint64_t m_num_consumers; // number of threads that will consume reads
  uint64_t m_perThreadReads;
  uint64_t m_numThreads;
  std::string m_fastqFilename;
  std::string m_indexFileName;
  bool m_isRunning = false;
  std::atomic_uint64_t chunk_counter_{0};

  // Helper functions
  int loadIndex(const std::string &indexFileName);
  uint64_t getMaxBufLen();
};

#include "parser.inl"
