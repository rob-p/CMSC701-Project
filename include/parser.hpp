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
  ReadChunk(std::string fastq_file_name, struct deflate_index* idx, uint64_t token, Chunk chunk) : token_(token), chunk_range_(chunk) {
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
    uint64_t last_record_id = (end_chunk+1 < rec_boundaries.size()) ? rec_boundaries.at(end_chunk).first_record_in_chunk : idx_ptr_->total_record_count;
    last_rec_ = last_record_id;
    in_stream_.reset(new KseqIndexedGzipStreamIn(reader_.get()));
  }
 
  ReadChunk(ReadChunk&& other) = delete;
  ReadChunk(const ReadChunk& other) = delete;
  ReadChunk& operator=(const ReadChunk& other) = delete;
  ReadChunk& operator=(ReadChunk&& other) = delete;
  
  inline KseqIndexedGzipStreamIn& operator>>( klibpp::KSeq& rec )  {
    if (current_rec_ < last_rec_) {
      (*in_stream_) >> rec;
      current_rec_++;
    } else {
      in_stream_->set_eof();
    }
    return *in_stream_;
  }
 
  struct deflate_index* idx_ptr_{nullptr}; // this reader's local file ptr
  std::unique_ptr<KseqIndexedGzipStreamIn> in_stream_{nullptr};
  // the number of records we expect to parse from
  // this chunk.
  uint64_t current_rec_;
  uint64_t last_rec_;
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
  std::optional<ReadChunk> get_read_chunk();
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
