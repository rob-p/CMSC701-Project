ParrFQParser::~ParrFQParser() {
  if (m_isRunning) {
    stop();
  }
}

int ParrFQParser::init(const std::string& fastqFilename, const std::string& indexFileName, uint64_t num_consumers) {
  m_fastqFilename = fastqFilename;
  m_indexFileName = indexFileName;
  m_perThreadReads = 0;
  m_numThreads = 0;
  return 0;
}

uint64_t ParrFQParser::get_num_chunks() {
  return m_index->have;
}

ReadChunk ParrFQParser::get_read_chunk() {
  struct deflate_index* idx = m_index.get();
  return ReadChunk(m_fastqFilename, idx);
}

// takes a read chunk, given to us by an underlying 
// consumer thread, and fills it.
bool ParrFQParser::refill(ReadChunk& tlc) {
  uint64_t next_chunk = chunk_counter_++;
  tlc.chunk_num = next_chunk;
  if (next_chunk >= m_index->have) {
    // no more chunks!
    return false;
  }
  
  bool got_chunk = get_uncompressed_chunk(tlc.file_ptr_.get(), tlc.idx_ptr_, next_chunk, tlc.data_, tlc.expected_rec);
  if (got_chunk) {
    tlc.in_stream_.reset(new KseqCharStreamIn(reinterpret_cast<const char*>(tlc.data_.data()), tlc.data_.size()));
  } else {
    std::cerr << "failed to get chunk";
  }
  return got_chunk;
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
  chunk_counter_ = 0;
  m_isRunning = true;
  return 0;
}

int ParrFQParser::stop() {
  if (m_isRunning == false) {
    std::cout << "ParrFQParser is not running" << std::endl;
    return -1;
  }
  m_isRunning = false;
  return 0;
}

int ParrFQParser::loadIndex(const std::string& indexFileName) {
  struct deflate_index* index = NULL;
  FILE *indexFile = fopen(indexFileName.c_str() , "rb");
  if (indexFile == NULL) {
    fprintf(stderr, "Could not open index for reading\n");
    return -1;
  }
  //int len = deflate_index_load(indexFile, &index);
  gzFile index_file_gzip = gzopen(indexFileName.c_str(), "rb");
  int len = deflate_index_load_gzip(index_file_gzip, &index);
  fclose(indexFile);
  if (len < 0) {
    fprintf(stderr, "Could not load index %d\n", len);
    return -1;
  }
  m_index.reset(index);
  return 0;
}
