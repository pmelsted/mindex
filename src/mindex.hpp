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

  size_t K;
  size_t maxdeg;
  size_t k;
  size_t w;
  size_t threads;

  vector<string> files;
  vector<string> input;
  string index;

  Mindex_opt() : threads(1), k(31), K(100), maxdeg(0), w(49), ratio(0.01) {}
};

class Mindex {

public:

  Mindex(const size_t k_, const size_t w_);
  Mindex();

  bool build(const Mindex_opt& opt);
  bool countData(const Mindex_opt& opt);
  bool findExistent(const Mindex_opt& opt);
  bool probData(vector<double> &prob, const Mindex_opt& opt);
  bool writeToFile(string fn, const Mindex_opt& opt) const;
  bool loadFromFile(string fn, Mindex_opt& opt);


  inline uint64_t scramble(uint64_t x) {
    return 11400714785074694791ULL * x + 14029467366897019727ULL;
  }


  int k;
  int w;
  uint64_t limit;
  bool invalid;

  static const int tiny_vector_sz = 1;

  // maps minimizer -> list of ids that contain k-mer
  KmerHashTable<vector<uint32_t>> min_table;
  // maps minimizer -> number of times seen
  KmerHashTable<uint32_t> counts;
  // list (ids) of list of minimizers for each id
  vector<vector<Kmer>> minimizers;

};

#endif
