#ifndef MINDEX_COMMON_HPP
#define MINDEX_COMMON_HPP

#include <cassert>
#include <algorithm>
#include <stdint.h>


#define MINDEX_VERSION "0.3"

static const char alpha[4] = {'A','C','G','T'};

inline size_t rndup(size_t v) {

  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  v++;

  return v;
}

inline uint32_t rndup(uint32_t v) {

  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;

  return v;
}



#endif // MINDEX_COMMON_HPP
