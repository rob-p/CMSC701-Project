#include <iostream>
#include <zlib.h>
#include "kseq++/kseq++.hpp"

using namespace klibpp;

int main(int argc, char* argv[])
{
  KSeq record;
  gzFile fp = gzopen(argv[1], "r");
  auto ks = make_kstream(fp, gzread, mode::in);
  //auto ks = KStream(fp, gzread, mode::in);  // C++17
  // auto ks = KStreamIn(fp, gzread);  // C++17
  int nr = 0;
  while ((ks >> record) && (nr < 5)) {
    std::cout << record.name << std::endl;
    if (!record.comment.empty()) std::cout << record.comment << std::endl;
    std::cout << record.seq << std::endl;
    if (!record.qual.empty()) std::cout << record.qual << std::endl;
    ++nr;
  }
  gzclose(fp);
}
