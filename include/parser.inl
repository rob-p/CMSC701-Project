ParrFQParser::~ParrFQParser() {
  if (m_isRunning) {
    stop();
  }
}

int ParrFQParser::init (const std::string& fastqFilename, const std::string& indexFileName, uint64_t perThreadReads, uint64_t numThreads) {
  m_fastqFilename = fastqFilename;
  m_indexFileName = indexFileName;
  m_perThreadReads = perThreadReads;
  m_numThreads = numThreads;
  m_currMaxOffset = 0;

  m_readQueue = std::make_unique<moodycamel::ConcurrentQueue<klibpp::KSeq>>();

  for (uint64_t i = 0; i < m_numThreads; ++i) {
    m_producerTokens.emplace_back(std::make_unique<moodycamel::ProducerToken>(*m_readQueue));
  }

  return 0;
}

int ParrFQParser::parse_reads(uint64_t threadId, moodycamel::ProducerToken* token, uint64_t maxBufLen) {
  // Load the index file
  // Create a shallow copy of the index object for each thread
  // Do not call deflate_index_free on this object - TODO: This might be a memory leak
  struct deflate_index* indexPerThread = new struct deflate_index(*m_index);

  unsigned char* buf = new unsigned char[maxBufLen];
  int got;

  while (true) {
    uint64_t startRecordIdx = m_currMaxOffset.fetch_add(this->m_perThreadReads);
    if (startRecordIdx >= indexPerThread->num_records) {
      // All records in the file have been or are being processed
      // This thread has nothing more to do
      break;
    }

    off_t numRecords = this->m_perThreadReads;
    std::tie(buf, got) = read_index(m_fastqFilename.c_str(), indexPerThread, startRecordIdx, numRecords, buf);
    if (got < 0) {
      fprintf(stderr, "[%llu] Parsing failed failed: %s error\n", threadId,
              got == Z_MEM_ERROR ? "out of memory" : "input corrupted");
      --m_numActiveThreads;
      return -1;
    }
    KseqCharStreamIn in(reinterpret_cast<const char*>(buf), got);

    klibpp::KSeq rec;
    while (in >> rec) {
      m_readQueue->enqueue(*token, rec);
    }
  }

  delete[] buf;
  --m_numActiveThreads;
  return 0;
}

int ParrFQParser::start() {
  if (m_isRunning == true) {
    std::cout << "ParrFQParser is already running" << std::endl;
    return -1;
  }

  // Load the index
  int ret = loadIndex(m_indexFileName);
  if (ret != 0) return ret;


  // Get the maximum buffer length required to store the reads based on the index
  uint64_t maxBufLen = getMaxBufLen();
  if (maxBufLen == 0) {
    std::cout << "Error: Could not get the maximum buffer length" << std::endl;
    return -1;
  }

  // TODO: Save the result of each thread in a vector and return it
  for (uint64_t i = 0; i < m_numThreads; ++i) {
    ++m_numActiveThreads;
    m_workers.emplace_back(new std::thread([this, i, maxBufLen]() {
      this->parse_reads(i, m_producerTokens[i].get(), maxBufLen);
    }));
  }
  m_isRunning = true;
  return 0;
}

int ParrFQParser::stop() {
  if (m_isRunning == false) {
    std::cout << "ParrFQParser is not running" << std::endl;
    return -1;
  }
  for (uint64_t i = 0; i < m_numThreads; ++i) {
    m_workers[i]->join();
  }
  m_isRunning = false;
  return 0;
}

moodycamel::ConsumerToken ParrFQParser::getConsumerToken() {
  return moodycamel::ConsumerToken(*m_readQueue);
}

bool ParrFQParser::getRead(moodycamel::ConsumerToken& token, klibpp::KSeq& rec) {
  // Returns true if a read was found, false otherwise. false indicates that the queue is empty
  bool found = m_readQueue->try_dequeue(token, rec);
  return found;
}

bool ParrFQParser::checkFinished() {
  // Checks if the parser has finished parsing all the reads and all threads have enqued all the reads
  return m_isRunning == true && m_numActiveThreads == 0;
}

int ParrFQParser::loadIndex(const std::string& indexFileName) {
  struct deflate_index* index = NULL;
  FILE *indexFile = fopen(indexFileName.c_str() , "rb");
  if (indexFile == NULL) {
    fprintf(stderr, "Could not open index for reading\n");
    return -1;
  }
  int len = deflate_index_load(indexFile, &index);
  fclose(indexFile);
  if (len < 0) {
    fprintf(stderr, "Could not load index %d\n", len);
    return -1;
  }
  m_index.reset(index);
  return 0;
}

uint64_t ParrFQParser::getMaxBufLen() {
  if (m_index == nullptr) {
    return 0;
  }
  uint64_t maxBufLen = 0;
  for (uint64_t i = 0; i < m_index->num_records; i+=m_perThreadReads) {
    // Find all possible lengths of the buffer given the m_perThreadReads and the index to pre-allocate the buffer per thread
    uint64_t bufLen = get_read_len(m_index.get(), i, m_perThreadReads);
    if (bufLen > maxBufLen) {
      maxBufLen = bufLen;
    }
  }
  return maxBufLen;
}