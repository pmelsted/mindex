#include "mindex.hpp"
#include <fstream>
#include <algorithm>
#include <sstream>

Mindex::Mindex(const size_t k_, const size_t w_): k(k_), w(w_), invalid(false) {}
Mindex::Mindex(): invalid(true) {}

bool Mindex::build(const Mindex_opt& opt) {
  size_t file_id = 0;
  size_t num_min = 0;
  string s,name;
  int nseq = 0, sseq = 0;

  limit = (uint64_t)( ((double)numeric_limits<uint64_t>::max()) * opt.ratio);
  minimizers.clear();
  minimizers.resize(opt.files.size());  

  // Main worker thread
  auto minz_worker_function = [&](const vector<string>& vseq, vector<Minimizer> &out) {

    for (const auto &seq : vseq) {
      const char* str = seq.c_str();
      const int len = seq.size();
      
      minHashIterator<RepHash> it_min(str, len, w, k, RepHash(), false), it_min_end;
      uint64_t lasthash = it_min.getHash()+1;
      // for every minimizer in the sequence

      for (; it_min != it_min_end; ++it_min) {
        // downsample by a factor of r      
        uint64_t h = it_min.getHash();
        if (lasthash != h && scramble(h) <= limit) {             
          out.push_back(Minimizer(seq.c_str() + it_min.getPosition()).rep());
          lasthash = h;
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

  // TODO: make this an option and figure out better defaults
  size_t K = 100;
  size_t max_deg = maxocc;

  vector<size_t> degree(nfiles);
  vector<vector<pair<size_t,Minimizer>>> full_minimizers(nfiles);
  auto deg_sorter = [&](const pair<size_t,Minimizer>& a, const pair<size_t,Minimizer>& b) -> bool {
      return a.first < b.first; // sort by degree
    };

  for (size_t i = 0; i < nfiles; i++) {
    for (const auto &x : minimizers[i]) {
      const auto px = occ_table.find(x);
      if (px != occ_table.end()) {
        if ((*px) <= max_deg) {
          full_minimizers[i].push_back({*px,x});
        }
      }
    }
    std::sort(full_minimizers[i].begin(), full_minimizers[i].end(), deg_sorter);
    minimizers[i].clear();
  }

  // minimizers is empty, full_minimizers contains all info needed
  
  // basic strategy, start with unique
  int num_ins = 0;
  for (size_t i = 0; i < nfiles; i++) {
    for (auto &x : full_minimizers[i]) {
      if (x.first == 1) {
        minimizers[i].push_back(x.second);
        degree[i]++;
        num_ins++;
        min_table.insert(x.second, {(uint32_t)i});
      }
    }
  }
  cerr << "degree " << 1 << ", inserted " << num_ins << " minimizers" << endl;


  // greedy strategy, prioritize by degree, pick all that are needed
  bool done = false;
  size_t deg = 2; // start at degree 2
  while (!done && deg <= max_deg) {
    num_ins = 0;
    done = true;
    for (size_t i = 0; i < nfiles; i++) {
      if (degree[i] < K) {
        done = false;
        const auto &v = full_minimizers[i];
        // for each minimizer of degree deg
        auto start = lower_bound(v.begin(), v.end(), make_pair(deg,Minimizer()), deg_sorter);
        for (; start != v.end(); ++start) {
          if (start->first > deg) {
            break;
          }
          auto x = start->second;
          degree[i]++;
          num_ins++;
          minimizers[i].push_back(x);
          auto px = min_table.find(x);
          if (px == min_table.end()) {
            // newly inserted
            min_table.insert(x, {(uint32_t)i});
          } 
        }
      }
    }

    cerr << "degree " << deg << ", inserted " << num_ins << " minimizers" << endl;
    deg++; // increase the degree
  }
  if (!done && num_ins == 0) {
    cerr << "Impossible to fill table, increase max_deg and retry" << endl;
  }
  
  size_t sdeg = 0;
  for (size_t i = 0; i < nfiles; i++) {
    //cerr << i << "\t" << degree[i] << "\n";
    sdeg += degree[i];
  }
  cerr << endl;
  cerr << "Inserted " << min_table.size() <<   " minimizers into the graph " << endl;


  for (auto &x : min_table) {
    x.clear();
  }

  for (size_t i = 0; i < nfiles; i++) {
    for (const auto &x : minimizers[i]) {
      auto px = min_table.find(x);
      if (px == min_table.end()) {
        cerr << "Error" << endl;
        return false;
      }
      px->push_back((uint32_t)i);
    }
  }

  size_t sdeg2 = 0;
  for (auto x : min_table) {
    sdeg2 += x.size();
  }
  cerr << "Left degree = " << sdeg << ", right degree " << sdeg2 << endl;
  
  // TODO write to index and take output as option

  

  return writeToFile("index.txt", opt);;
}

bool Mindex::writeToFile(string fn, const Mindex_opt& opt) const{
  size_t nfiles = opt.files.size();
  ofstream of(fn);
  of << "#k,w,limit,nfiles\n";
  of << k << "\n" << w << "\n" << limit << "\n" << nfiles << "\n";
  of << "#files\n";
  for (size_t i = 0; i < nfiles; i++) {
    of << i << "\t" << opt.files[i] << "\n";
  }
  of << "#minimizers\n";
  for (size_t i = 0; i < nfiles; i++) {
    const auto &v = minimizers[i];
    of << i << "\t" << v.size() << "\n";
    for (const auto &x : v) {
      of << x.toString() << "\n";
    }
  }
  of << "#min_table\n";
  of << min_table.size() << "\n";

  const auto end = min_table.end();
  for (auto x = min_table.begin(); x != end; ++x) {
    of << x.getKey().toString() << "\t" << x->size() << "\n";
    const auto v = *x;
    for (const auto t : v) {
      of << t << "\n";
    }
  }
  of.close();
  return true;
}

bool Mindex::loadFromFile(string fn, Mindex_opt& opt) {
  ifstream in(fn);
  string line,tl; // temporary
  int stage = 0;
  size_t nfiles;
  while (in) {
    if (in.peek() == '#') {
      getline(in,line); // skip
    } else {
      if (stage == 0) {
        // read basic params
        in >> opt.k >> opt.w >> limit >> nfiles;
        getline(in,line);
        Kmer::set_k(opt.k);
        Minimizer::set_g(opt.k);
        stage = 1;        
      } else if (stage == 1) {
        // read files names
        opt.files.clear();
        
        minimizers.resize(nfiles);
        for (size_t i = 0; i < nfiles; i++) {
          size_t j;
          getline(in,line);
          stringstream ss(line);
          ss >> j >> tl;
          opt.files.push_back(tl);          
        }
        stage = 2;        
      } else if (stage == 2) {
        // read minimizers
        for (size_t i = 0; i < nfiles; i++) {
          size_t tmp, num;
          getline(in,line);
          stringstream ss(line);
          ss >> tmp >> num;
          auto &v = minimizers[i];
          for (size_t j = 0; j < num; j++) {
            getline(in,line);
            v.push_back(Minimizer(line.c_str()));
          }
        }
        stage = 3;
      } else if (stage == 3) {
        size_t num_min,sz;
        getline(in,line);
        stringstream ss(line);
        ss >> num_min;
        min_table.reserve(num_min);

        for (size_t i = 0; i < num_min; i++) {
          getline(in,line);
          stringstream ss(line);
          ss >> tl >> sz;
          Minimizer x(tl.c_str());
          auto p = min_table.insert(x, {});
          auto &v = *(p.first);
          v.reserve(sz);
          uint32_t t;
          for (size_t j = 0; j < sz; j++) {
            in >> t;
            v.push_back(t);
          }          
        }
        stage = 4;
        break;
        // break here
      }
    }  
  }
  if (stage == 4) {
    invalid = false;
  }
  return !invalid;
}
