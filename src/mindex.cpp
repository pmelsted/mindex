#include "mindex.hpp"
#include <algorithm>

Mindex::Mindex(const size_t k_, const size_t w_): k(k_), w(w_), invalid(false) {}

bool Mindex::build(const Mindex_opt& opt) {
  size_t file_id = 0;
  size_t num_min = 0;
  string s,name;
  int nseq = 0, sseq = 0;

  uint64_t limit = (uint64_t)( ((double)numeric_limits<uint64_t>::max()) * opt.ratio);
  minimizers.clear();
  minimizers.resize(opt.files.size());  

  // Main worker thread
  auto minz_worker_function = [&](const vector<string>& vseq, vector<Minimizer> &out) {

    for (const auto &seq : vseq) {
      const char* str = seq.c_str();
      const int len = seq.size();
      
      minHashIterator<RepHash> it_min(str, len, w, k, RepHash(), false), it_min_end;

      // for every minimizer in the sequence

      for (; it_min != it_min_end; ++it_min) {
        // downsample by a factor of r      
        if (scramble(it_min.getHash()) <= limit) {             
          out.push_back(Minimizer(seq.c_str() + it_min.getKmerPosition()).rep());
        }
      }
    }
  };

  // goes through every file in batched mode
  size_t nfiles = opt.files.size();
  for (size_t i = 0; i < nfiles; ) {
    size_t batch = std::min(nfiles-i, opt.threads);
    vector<thread> workers;

    for (size_t b = 0; b < batch; b++,i++) {
      // for each file, read all the sequences
      vector<string> readv;
      gzFile fp = gzopen(opt.files[i].c_str(), "r");
      kseq_t* kseq = kseq_init(fp);

      int r;
      while ((r = kseq_read(kseq)) >= 0) {
        readv.push_back(kseq->seq.s);
        nseq++;
        sseq += kseq->seq.l;
      }
      kseq_destroy(kseq);
      gzclose(fp);
      kseq = nullptr;
      // process in thread
      workers.push_back(thread(minz_worker_function, readv, std::ref(minimizers[i])));
    }
    
    for (auto &t : workers) {
      t.join();
    }
  }

  // figure out what to do with output
  for (auto &x : minimizers) {
    num_min += x.size();
  }
  cerr << "Number of files: " << nfiles <<", number of sequences: " << nseq << ", total bp: " << sseq  << std::endl;
  cerr << "Number of minimizers inserted: " << num_min << endl;

  // count the number of occurrances for each minimizer
  MinimizerHashTable<uint32_t> occ_table;
  occ_table.reserve(num_min);
  for (const auto &minv : minimizers) {
    vector<Minimizer> cp(minv);
    std::sort(cp.begin(), cp.end());
    auto last = unique(cp.begin(), cp.end());
    cp.erase(last, cp.end());
    for (const auto x : cp) {
      auto it = occ_table.find(x);
      if (it != occ_table.end()) {
        (*it)++;
      } else {
        occ_table.insert(x,1);
      }
    }    
  }

  uint32_t maxocc = 0;
  for (const auto &x : occ_table) {
    if (maxocc < x) {
      maxocc = x;
    }
  }

  /*
  vector<int> histo(maxocc+1);
  for (const auto &x : occ_table) {
    ++histo[x];
  }  

  cerr << "Histogram:" << std::endl;
  for (uint32_t i =0; i <= maxocc; i++) {
    cerr << i << "\t" << histo[i] << "\n";
  }
  cerr << std::endl;
  */

  size_t K = 1000;
  
  


  return true;
}
