#pragma once
#include "concurrentqueue/concurrentqueue.h"
#include "kseq++/seqio.hpp"
#include "zran.hpp"
#include "kseqcharstream.hpp"
#include <atomic>
#include <functional>
#include <memory>
#include <stdio.h>
#include <thread>
#include <vector>
#include <optional>


// Represents a contiguous chunk [start, end)
struct Chunk {
    size_t start;
    size_t end;
    
    size_t size() const { return end - start; }
    bool empty() const { return size() == 0; }
};


class ReadChunk {
public:
  ReadChunk(std::string fastq_file_name, struct deflate_index* idx, uint64_t token, Chunk chunk) : data_(std::make_unique<std::vector<char>>(1, 0)), token_(token), chunk_range_(chunk) {
    /*
    FILE* fp = fopen(fastq_file_name.c_str(), "rb");
    if (fp != nullptr) {
      file_ptr_ = std::unique_ptr<FILE, std::function<void(FILE*)>>(fp, [](FILE* f) {fclose(f);});
    }
    */
    idx_ptr_ = new struct deflate_index(*idx);
    size_t start_chunk = chunk_range_.start;
    size_t end_chunk = chunk_range_.end;
    reader_ = open_gzip_at_checkpoint(fastq_file_name.c_str(), idx_ptr_, start_chunk);

    // get the byte offset of the first read in this chunk
    auto& rec_boundaries = *idx_ptr_->record_boundaries;
    uint64_t first_read_byte = rec_boundaries.at(start_chunk).byte_offset;
    uint64_t first_record_id = rec_boundaries.at(start_chunk).first_record_in_chunk;

    // we want to discard bytes from the current uncompressed_offset until the first_read_byte
    off_t amount_to_discard = first_read_byte - reader_->uncompressed_offset;
    // if we don't have to discard anything, avoid all of this work
    if (amount_to_discard != 0) { 
      std::vector<char> discard_buffer(amount_to_discard); 
      ptrdiff_t bytes_read = gzip_read(reader_.get(), const_cast<char*>(discard_buffer.data()), amount_to_discard);
      if (bytes_read < 0) {
        std::cerr << "WORKER " << token_ << ", ATTEMPTING TO DISCARD INITIAL " << amount_to_discard 
                  << " BYTES; ERROR OPENING AND READING FROM CHECKPOINT " << start_chunk << "\n";
        return;
      }
    } 
    current_rec_ = first_record_id;

    uint64_t last_record_id = (end_chunk+1 < rec_boundaries.size()) ? rec_boundaries.at(end_chunk+1).first_record_in_chunk : idx_ptr_->total_record_count;

    last_rec_ = last_record_id;
    //in_stream_.reset(new KseqIndexedGzipStreamIn(reader_.get()));
    std::cerr << "SUCCESSFULLY CREATED READ CHUNK " << token_ << ", starting and checkpoint " << start_chunk << "\n";
  }
 
  ReadChunk(ReadChunk&& other) = delete;
  ReadChunk(const ReadChunk& other) = delete;
  ReadChunk& operator=(const ReadChunk& other) = delete;
  ReadChunk& operator=(ReadChunk&& other) = delete;
  /*
  ~ReadChunk() {
    std::cerr << "running ReadChunk destructor\n";
    reader_.reset(nullptr);
    idx_ptr_ = nullptr;
    in_stream_.reset(nullptr);
    data_.reset(nullptr);
    chunk_range_.start = chunk_range_.end = 0;
  }

  ReadChunk(ReadChunk&& other)  {
    std::cerr << "using move!\n";
    idx_ptr_ = other.idx_ptr_; other.idx_ptr_=nullptr;
    std::swap(in_stream_, other.in_stream_);
    std::swap(data_, other.data_);
    std::swap(reader_, other.reader_);
    other.in_stream_.reset(nullptr);
    other.reader_.reset(nullptr);
    other.data_.reset(nullptr);
    chunk_range_ = other.chunk_range_;
    expected_rec = other.expected_rec;
    chunk_num = other.chunk_num;
    token_ = other.token_;
  }
  */
  /*
  ReadChunk(ReadChunk& other) = delete;
  ReadChunk() = delete;
  ReadChunk& operator=(ReadChunk& other) = delete;

  ReadChunk& operator=(ReadChunk&& other)  {
    idx_ptr_ = other.idx_ptr_; other.idx_ptr_=nullptr;
    in_stream_ = std::move(other.in_stream_); other.in_stream_.reset();
    std::swap(data_, other.data_);
    std::swap(reader_, other.reader_);
    expected_rec = other.expected_rec;
    chunk_num = other.chunk_num;
    token_ = other.token_;
    return *this;
  }
  */
 
  /*
  inline KseqIndexedGzipStreamIn& operator>>( klibpp::KSeq& rec )  {
    if (current_rec_ < last_rec_) {
      (*in_stream_) >> rec;
      current_rec_++;
    } else {
      in_stream_->set_eof();
    }
    return *in_stream_;
  }
  */
  inline KseqCharStreamIn& operator>>( klibpp::KSeq& rec )  {
      (*in_stream_) >> rec;
      return *in_stream_;
  }
  //KseqCharStreamIn& kget_seq_stream() { return *in_stream_; }
  
  inline bool refill() {
    // first check that there are chunks remaining 
    if (chunk_range_.start < chunk_range_.end) {
      // advance the chunk
      chunk_range_.start++;
      // the current buffer should span from the current uncompressed_offset until the 
      // first read of the next chunk.
      // get the next record boundary byte offset
      std::vector<record_checkpoint>* rec_boundaries = idx_ptr_->record_boundaries;
      if (rec_boundaries == nullptr) { std::cerr << "unexpected!\n"; return false; }
      uint64_t next_read_byte = (chunk_range_.start < rec_boundaries->size()) ? rec_boundaries->at(chunk_range_.start).byte_offset : idx_ptr_->length;
      size_t bytes_to_read = next_read_byte - reader_->uncompressed_offset;
      //std::cerr << "ATTEMPTING TO REFILL CHUNK; reading bytes " << reader_->uncompressed_offset << " to " << next_read_byte << "\n";
      data_->resize(bytes_to_read);
      ptrdiff_t bytes_read = gzip_read(reader_.get(), const_cast<char*>(data_->data()), bytes_to_read);
      if (bytes_read == bytes_to_read) {
        in_stream_.reset(new KseqCharStreamIn(reinterpret_cast<const char*>(data_->data()), data_->size()));
        return true;
      } else {
        std::cerr << "worker " << token_ << " expected to read " << bytes_to_read << " bytes, but read " << bytes_read  << "\n";
      }
    }
    return false;
  }
  
 
  struct deflate_index* idx_ptr_{nullptr}; // this reader's local file ptr
  //std::unique_ptr<KseqIndexedGzipStreamIn> in_stream_{nullptr};
  std::unique_ptr<KseqCharStreamIn> in_stream_{nullptr};
  std::unique_ptr<std::vector<char>> data_{nullptr};
  // the number of records we expect to parse from
  // this chunk.
  uint64_t current_rec_;
  uint64_t last_rec_;
  // the underlying assigned chunk num
  uint64_t chunk_num;
  // the unique token identifying the worker
  // to which this ReadChunk belongs
  uint64_t token_;
  // this thread's local reader
  std::unique_ptr<GzipStreamReader> reader_{nullptr};
  // 
  Chunk chunk_range_;
};

/*
class ParrFQPairParser {
public:
  ParrFQPairParser()
      : m_index(nullptr,
                [](struct deflate_index *p) { deflate_index_free(p); }) {}

  ~ParrFQPairParser();

  int init(const std::string &fastqFilename, const std::string &indexFileName,
           uint64_t num_consumers);


  // Start and stop the parser
  int start();
  int stop();
  uint64_t get_num_chunks();

  // Consumer functions
  ReadPairChunk get_read_chunk();
  bool refill(ReadPairChunk& tlc);

private:
  std::unique_ptr<struct deflate_index,
                  std::function<void(struct deflate_index *)>>
      m_index;
  std::unique_ptr<struct deflate_index,
                  std::function<void(struct deflate_index *)>>
      m_index2;


  uint64_t m_num_consumers; // number of threads that will consume reads
  uint64_t m_perThreadReads;
  uint64_t m_numThreads;
  std::string m_fastqFilename;
  std::string m_indexFileName;
  std::string m_fastqFilename2;
  std::string m_indexFileName3;
  bool m_isRunning = false;
  std::atomic_uint64_t chunk_counter_{0};

  // Helper functions
  int loadIndex(const std::string &indexFileName);
  uint64_t getMaxBufLen();
};
*/

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
  uint64_t get_num_reads();

  // Consumer functions
  std::unique_ptr<ReadChunk> get_read_chunk();
  bool refill(ReadChunk& tlc);
  // called ONCE by each thread to get it's token
  std::optional<uint64_t> get_token();
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
  std::atomic_uint64_t token_counter_{0};
  std::vector<Chunk> chunk_ranges_;
  // Helper functions
  int loadIndex(const std::string &indexFileName);
  uint64_t getMaxBufLen();
};

#include "parser.inl"
