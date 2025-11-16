
ParrFQParser::~ParrFQParser() {
  if (m_isRunning) {
    stop();
  }
}

int ParrFQParser::init(const std::string& fastqFilename, const std::string& indexFileName, uint64_t num_consumers) {
  m_fastqFilename = fastqFilename;
  m_indexFileName = indexFileName;
  m_perThreadReads = 0;
  m_numThreads = num_consumers;
  return 0;
}

uint64_t ParrFQParser::get_num_chunks() {
  return m_index->have;
}

uint64_t ParrFQParser::get_num_reads() {
  return m_index->total_record_count;
}

std::unique_ptr<ReadChunk> ParrFQParser::get_read_chunk() {
  uint64_t curr_token = token_counter_.fetch_add(1);
  struct deflate_index* idx = m_index.get();
  return (curr_token >= m_numThreads) ? nullptr : std::make_unique<ReadChunk>(m_fastqFilename, idx, curr_token, chunk_ranges_[curr_token]);
}

// takes a read chunk, given to us by an underlying 
// consumer thread, and fills it.
bool ParrFQParser::refill(ReadChunk& tlc) {
  if (tlc.token_ >= m_numThreads) {
    std::cerr << "invalid worker token (" << tlc.token_ << "), but only " << m_numThreads << " workers were registered.\n";
    return false;
  }
  return tlc.refill();
}

// Distributes N items among M threads in contiguous chunks
// Minimizes the maximum chunk size (load balancing)
std::vector<Chunk> distribute_chunks(size_t N, size_t M) {
    assert(M > 0 && "Number of chunks must be positive");
    assert(N > 0 && "Number of items must be positive");
    
    std::vector<Chunk> chunks;
    chunks.reserve(M);
    
    // If we have fewer items than chunks, some threads get nothing
    if (N < M) {
        for (size_t i = 0; i < N; ++i) {
            chunks.push_back({i, i + 1});
        }
        for (size_t i = N; i < M; ++i) {
            chunks.push_back({N, N});  // Empty chunk
        }
        return chunks;
    }
    
    // Base size for each chunk and number of chunks that get +1
    size_t base_size = N / M;
    size_t remainder = N % M;
    
    // First 'remainder' chunks get (base_size + 1) items
    // Remaining chunks get base_size items
    size_t current_start = 0;
    
    for (size_t i = 0; i < M; ++i) {
        size_t chunk_size = base_size + (i < remainder ? 1 : 0);
        chunks.push_back({current_start, current_start + chunk_size});
        current_start += chunk_size;
    }
    
    return chunks;
}

int ParrFQParser::start() {
  if (m_isRunning == true) {
    std::cout << "ParrFQParser is already running" << std::endl;
    return -1;
  }

  // Load the index
  int ret = loadIndex(m_indexFileName);
  if (ret != 0) return ret;

  // determine a work plan for the different threads 
  // we create a vector of cursors, each keyed on the 
  // consumer token.
  const size_t total_blocks = m_index->list->size() + 1;
  chunk_ranges_ = distribute_chunks(total_blocks, m_numThreads);
  std::cerr << "chunk ranges: {";
  for (size_t i = 0; i < chunk_ranges_.size(); ++i) {
    std::cerr << "[" << chunk_ranges_[i].start << ", " << chunk_ranges_[i].end << ")";
    if (i != chunk_ranges_.size() - 1) { std::cerr << ", "; }
  }
  std::cerr << "}\n";

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
