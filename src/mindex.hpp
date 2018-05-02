#ifndef MINDEX_MINDEX_H
#define MINDEX_MINDEX_H

#include <iostream>
#include <thread>
#include <atomic>
#include <zlib.h>


#include "KmerIterator.hpp"
#include "KmerHashTable.hpp"
#include "minHashIterator.hpp"
#include "RepHash.hpp"
#include "Kmer.hpp"
#include "TinyVector.hpp"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
#endif


using namespace std;

struct Mindex_opt {

  bool verbose;
  double ratio;

  size_t k;
  size_t w;
  size_t threads;

  vector<string> files;

  Mindex_opt() : threads(1), k(31), w(49), ratio(0.01) {}
};

class Mindex {

public:

  Mindex(const size_t k_, const size_t w_);

  bool build(const Mindex_opt& opt);

  inline uint64_t scramble(uint64_t x) {
    return 11400714785074694791ULL * x + 14029467366897019727ULL;
  }
private:

  int k;
  int w;

  bool invalid;

  static const int tiny_vector_sz = 1;

  // <Number of ids, list of ids> -> Number of id is always right but the list of ids might have been cut off
  typedef pair<size_t, tiny_vector<size_t, tiny_vector_sz>> p_id_t;
  MinimizerHashTable<vector<uint32_t>> min_table;
  vector<vector<Minimizer>> minimizers;

};

#endif
