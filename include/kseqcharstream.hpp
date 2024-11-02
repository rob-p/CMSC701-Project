#pragma once
#include <zlib.h>
#include <sstream>

#include "kseq++/kseq++.hpp"

struct CharBuffer {
  const char* buffer;
  size_t size;
  int index;

  CharBuffer() = delete;
  CharBuffer(const CharBuffer& other) = default;
  CharBuffer& operator=(const CharBuffer& other) = default;
  CharBuffer& operator=(CharBuffer&& other) = default;

  CharBuffer(CharBuffer&& other) {
    buffer = other.buffer;
    other.buffer = nullptr;
    size = other.size;
    index = other.index;
  }
  CharBuffer( const char* buffer, size_t size ) : buffer(buffer), size(size), index(0) {}

  int copy(void* data, unsigned int bytesToCopy) {
    if (size - index < bytesToCopy) {
      bytesToCopy = size - index;
    }
    std::memcpy(data, buffer + index, bytesToCopy);
    index += bytesToCopy;
    return bytesToCopy;
  }
};

// Custom stream reader for reading from a character buffer instead of a file
class KseqCharStreamIn : public klibpp::KStreamIn<CharBuffer, int(*)(CharBuffer&, void*, unsigned int)> {
 public:
  using Base = klibpp::KStreamIn<CharBuffer, int(*)(CharBuffer&, void*, unsigned int)>;

  KseqCharStreamIn(const char* buffer, size_t size) : Base(CharBuffer(buffer, size), readChar, close) {}

  static int readChar(CharBuffer& buf, void* data, unsigned int size) {
    return buf.copy(data, size);
  }

  static int close(CharBuffer buf) {
    return 0;
  }
};
