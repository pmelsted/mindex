#include "mindex.hpp"
#include <fstream>
#include <algorithm>
#include <sstream>
#include <unordered_set>

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
  auto minz_worker_function = [&](const vector<string>& vseq, vector<Kmer> &out) {

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
          auto x = Kmer(seq.c_str() + it_min.getPosition()).rep();
          out.push_back(x);
          lasthash = h;
        }
      }      
    }

    unordered_set<Kmer, KmerHash> tmp(out.begin(), out.end());
    out.clear();
    out.assign(tmp.begin(), tmp.end());
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

  for (auto &x : minimizers) {
    num_min += x.size();
  }
  cerr << "Number of files: " << nfiles <<", number of sequences: " << nseq << ", total bp: " << sseq  << std::endl;
  cerr << "Number of minimizers considered: " << num_min << endl;

  // count the number of occurrances for each minimizer
  KmerHashTable<uint32_t> occ_table;
  occ_table.reserve(num_min);
  for (const auto &minv : minimizers) {
    for (const auto x : minv) {      
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

  size_t K = opt.K;
  size_t max_deg =  (opt.maxdeg <= 0) ? maxocc : opt.maxdeg;

  if (max_deg < maxocc) {
    // maybe: clean up unneccessary minimizers, maybe remove options
  }

  vector<size_t> degree(nfiles);
  vector<vector<pair<size_t,Kmer>>> full_minimizers(nfiles);
  auto deg_sorter = [&](const pair<size_t,Kmer>& a, const pair<size_t,Kmer>& b) -> bool {
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

  occ_table.clear(); // not needed again

  // minimizers is empty, full_minimizers contains all info needed
  // basic strategy, start with unique
  int num_ins = 0;
  for (size_t i = 0; i < nfiles; i++) {
    for (auto &x : full_minimizers[i]) {
      if (x.first == 1) {
        degree[i]++;
        num_ins++;
        min_table.insert(x.second, {(uint32_t)i});
      } else {
        break; // exit early
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
    // go throught every file and see if we need to insert
    for (size_t i = 0; i < nfiles; i++) {
      if (degree[i] < K) {
        done = false;
        const auto &v = full_minimizers[i];
        // for each minimizer of degree deg
        int goal = (int) (K - degree[i]);
        auto start = lower_bound(v.begin(), v.end(), make_pair(deg,Kmer()), deg_sorter);
        for (; start != v.end() && goal > 0 ; ++start) {
          if (start->first > deg) {
            break;
          }
          auto x = start->second;
          goal--;
          auto px = min_table.find(x);
          if (px == min_table.end()) {
            // newly inserted
            min_table.insert(x, {});
          }
        }
      }
    }

    for (size_t i = 0; i < nfiles; i++) {
      const auto &v = full_minimizers[i];
      // for each minimizer of degree deg    
      auto start = lower_bound(v.begin(), v.end(), make_pair(deg,Kmer()), deg_sorter);
      for (; start != v.end(); ++start) {
        if (start->first > deg) {
          break;
        }
        auto x = start->second;
        auto px = min_table.find(x);
        if (px != min_table.end()) {
          minimizers[i].push_back(x);
          px->push_back(i);
          degree[i]++;
          num_ins++;
        }
      }
    }

    cerr << "degree " << deg << ", inserted " << num_ins << " minimizers" << endl;
    deg++; // increase the degree
  }
  if (!done && num_ins == 0) {
    cerr << "Impossible to fill table, increase max_deg and retry" << endl;
  }
  
  // check statistics
  size_t sdeg = 0;
  for (size_t i = 0; i < nfiles; i++) {
    sdeg += degree[i];
  }
  size_t sdeg2 = 0;
  for (auto x : min_table) {
    sdeg2 += x.size();
  }

  // fix the minimizers
  full_minimizers.clear();

  cerr << endl;
  cerr << "Inserted " << min_table.size() <<   " minimizers into the graph " << endl;

  
  cerr << "Left degree = " << sdeg << ", right degree " << sdeg2 << endl;
  
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
        k = opt.k;
        w = opt.w;
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
            v.push_back(Kmer(line.c_str()));
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
          Kmer x(tl.c_str());
          auto p = min_table.insert(x, {});
          auto &v = *(p.first);
          v.reserve(sz);
          uint32_t t;
          for (size_t j = 0; j < sz; j++) {
            in >> t;
            getline(in,line);
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

bool Mindex::countData(const Mindex_opt& opt) {
  // each file is independent for now

  size_t found = 0;
  // clear the counts
  counts.clear();
  counts.reserve(min_table.size());
  auto it_end = min_table.end();
  for (auto it = min_table.begin(); it != it_end; ++it) {
    counts.insert(it.getKey(), 0);
  }

  auto count_end = counts.end();
  minHashIterator<RepHash> it_min(w,k,RepHash()), it_min_end;
  for (auto fn : opt.input) {
    gzFile fp = gzopen(fn.c_str(), "r");
    kseq_t* kseq = kseq_init(fp);

    int r;
    while ((r = kseq_read(kseq)) >= 0) {
      const char *seq = kseq->seq.s;
      it_min.initString(kseq->seq.s, kseq->seq.l);
      uint64_t lasthash = it_min.getHash()+1;
      // for every minimizer in the sequence
      for (; it_min != it_min_end; ++it_min) {
        // downsample by a factor of r      
        uint64_t h = it_min.getHash();
        if (lasthash != h && scramble(h) <= limit) {             
          // check if we find this
          lasthash = h;
          auto it = counts.find(Kmer(seq + it_min.getPosition()).rep());
          if (it != count_end) {
            ++(*it);   
            ++found;       
          }          
        }
      }
    }
    kseq_destroy(kseq);
    gzclose(fp);
    kseq = nullptr;
  }

  cerr << "Found " << found << " hits for " << counts.size() << " minimizers" << endl;

  //print out assignment table
  for (auto it = counts.begin(); it != count_end; ++it) {
    if ((*it) > 1) {
      auto x = it.getKey();

      cerr << x.toString() << "\t";
      auto it2 = min_table.find(x);
      for (auto t : *(it2)) {
        cerr << t << ",";
      }
      cerr << "\t" << *it << "\n";
    }
  }

  return true;
}

bool Mindex::findExistent(const Mindex_opt& opt) {

  {
    vector<Kmer> rem;
    for (auto it = counts.begin(); it != counts.end(); ++it) {
      if (*it <= 1) {
        rem.push_back(it.getKey());
      }
    }
    for (auto &x : rem) {
      counts.erase(x);
    }
    rem.clear();
  }

  if (counts.empty()) {
    return false;
  }
  size_t nfiles = minimizers.size();
  vector<int> support(nfiles, 0);
  vector<int> num_uniq(nfiles,0);
  vector<int> count_uniq(nfiles,0);

  auto cnt_end = counts.end();
  auto min_end = min_table.end();

  for (auto it = counts.begin(); it != cnt_end; ++it) {
    auto x = it.getKey();
    auto px = min_table.find(x);
    if (px != min_end) {
      for (auto i : *px) {
        support[i] = 1;
      }
    }  
  }

  int round = 1;
  bool done = false;
  const int min_num_uniq = 2;
  const double min_ratio = 0.5;
  while (!done) {
    done = true;
    vector<size_t> rem;
    // count the number of unique hits
    for (auto it = min_table.begin(); it != min_end; ++it) {
      auto x = it.getKey();
      if (it->size() == 1) {
        size_t i = (*it)[0];
        ++num_uniq[i];
        auto px = counts.find(x);
        if (px != counts.end()) {
          if (*px > 1) {
            ++count_uniq[i];
          }
        }
      }
    }
    
    for (size_t i = 0; i < nfiles; i++) {
      if (support[i] > 0) {
        size_t c = count_uniq[i];
        size_t nm = num_uniq[i];
        cerr << i << "\t" << c << "\t" << nm << "\t" << c / double(nm) << endl;
        if (nm >= 4) {
          double ru = c / double(nm);
          if (ru < min_ratio) {
            rem.push_back(i);
          }
        }
      }
    }

    for (auto i : rem) {
      support[i] = 0;
      auto v = minimizers[i]; // explicit copy
      minimizers[i].clear();
      for (auto x : v) {        
        auto px = min_table.find(x);
        if (px != min_table.end()) {
          if (px->size() == 1) {
            min_table.erase(px);
            auto cx = counts.find(x);
            if (cx != counts.end()) {
              counts.erase(cx);
            }
          } else {
            auto &t = *px;
            // remove from min_table;
            auto ipos = find(t.begin(), t.end(), (uint32_t) i);
            if (ipos != t.end()) {
              t.erase(ipos); // linear, but who cares
            }
          }
        }
      }
    }

    done = rem.empty();
    ++round;
    if (!done) {
      cerr << "Round " << round << ", removed " << rem.size() << " entries, count size = " << min_table.size() << endl;
    }
  }
  cerr << endl;


  cerr << "Entries with support" << endl;
  int supp= 0;
  for (size_t i = 0; i < nfiles; i++) {
    if (support[i] > 0) {
      supp++;
      cerr << i << "\n";
    }
  }
  cerr << endl;
  cerr << "Total of " << supp << " entries" << endl;

  //print out assignment table
  cnt_end = counts.end();
  for (auto it = counts.begin(); it != cnt_end; ++it) {
    if ((*it) > 1) {
      auto x = it.getKey();

      cerr << x.toString() << "\t";
      auto it2 = min_table.find(x);
      for (auto t : *(it2)) {
        cerr << t << ",";
      }
      cerr << "\t" << *it << "\n";
    }
  }

  return true;
}

bool Mindex::probData(vector<double> &prob, const Mindex_opt& opt) {
  // returns the probability of presence for each id
  prob.clear();
  prob.resize(opt.files.size());
  if (counts.empty()) {
    return false;
  }

  // figure out how to estimate the probability

  return true;
}