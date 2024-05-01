#pragma once
#include "zran.hpp"
#include "kseq++/seqio.hpp"
#include "kseqcharstream.hpp"
#include "concurrentqueue/concurrentqueue.h"
#include <stdio.h>
#include <atomic>
#include <thread>
#include <vector>
#include <memory>
#include <functional>

class ParrFQParser {
 public:
  ParrFQParser() : m_index(nullptr, [](struct deflate_index* p) { deflate_index_free(p); }) {}

  ~ParrFQParser();

  int init(const std::string& fastqFilename, const std::string& indexFileName, uint64_t perThreadReads, uint64_t numThreads);

  // Main function that will be called by each thread to parse the reads
  int parse_reads(uint64_t threadId, moodycamel::ProducerToken* token, uint64_t maxBufLen);

  // Start and stop the parser
  int start();
  int stop();

  // Consumer functions
  moodycamel::ConsumerToken getConsumerToken();
  bool getRead(moodycamel::ConsumerToken& token, klibpp::KSeq& rec);
  bool checkFinished();

 private:
  // Each parser thread will check and update the current offset to claim a chunk of reads
  std::atomic<uint64_t> m_currMaxOffset;
  std::vector<std::unique_ptr<std::thread>> m_workers;
  std::unique_ptr<moodycamel::ConcurrentQueue<klibpp::KSeq>> m_readQueue;
  std::vector<std::unique_ptr<moodycamel::ProducerToken>> m_producerTokens;
  std::unique_ptr<struct deflate_index, std::function<void(struct deflate_index*)>> m_index;

  uint64_t m_perThreadReads;
  uint64_t m_numThreads;
  std::string m_fastqFilename;
  std::string m_indexFileName;
  bool m_isRunning = false;
  std::atomic<uint32_t> m_numActiveThreads = 0;

  // Helper functions
  int loadIndex(const std::string& indexFileName);
  uint64_t getMaxBufLen();
};

#include "parser.inl"
