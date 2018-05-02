#include "mindex.hpp"

Mindex::Mindex(const size_t k_, const size_t w_): k(k_), w(w_), invalid(false) {}

bool Mindex::build(const Mindex_opt& opt) {
  size_t file_id = 0;
  size_t num_min = 0;
  string s,name;
  int nseq = 0;

  uint64_t limit = (uint64_t)( ((double)numeric_limits<uint64_t>::max()) * opt.ratio);
  vector<vector<uint64_t>> minimizers(opt.files.size());

  // Main worker thread
  auto minz_worker_function = [&](const vector<string>& vseq, vector<uint64_t> &out) {

    for (const auto &seq : vseq) {
      const char* str = seq.c_str();
      const int len = seq.size();
      
      minHashIterator<RepHash> it_min(str, len, w, k, RepHash(), false), it_min_end;

      // for every minimizer in the sequence

      for (; it_min != it_min_end; ++it_min) {
        // downsample by a factor of r      
        if (scramble(it_min.getHash()) <= limit) {   
          out.push_back(it_min.getHash());
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
  cerr << "Number of files: " << nfiles <<", number of sequences: " << nseq << std::endl;
  cerr << "Number of minimizers inserted: " << num_min << endl;

  return true;
}
