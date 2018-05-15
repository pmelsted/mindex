#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>

#include <iostream>
#include <thread>
#include <atomic>

#include "mindex.hpp"

void parse_ProgramOptions_build(int argc, char **argv, Mindex_opt& opt) {

  const char* opt_string = "t:k:w:r:";

  static struct option long_options[] = {

    {"threads",         required_argument,  0, 't'},
    {"kmer-length",     required_argument,  0, 'k'},
    {"window",          required_argument,  0, 'w'}, 
    {"max-deg",         required_argument,  0, 'm'}, 
    {"left-deg",        required_argument,  0, 'l'}, 
    {"ratio-gmers",     required_argument,  0, 'r'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {

    case 't':
      opt.threads = atoi(optarg);
      break;
    case 'k':
      opt.k = atoi(optarg);
      break;
    case 'w':
      opt.w = atoi(optarg);
      break;
    case 'r':
      opt.ratio = atof(optarg);
      break;
    case 'l':
      opt.K = atoi(optarg);
      break;
    case 'm':
      opt.maxdeg = atoi(optarg);
      break;
    default:
      break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) opt.files.push_back(argv[optind++]);
}

bool check_ProgramOptions_build(Mindex_opt& opt) {

  bool ret = true;

  size_t max_threads = std::thread::hardware_concurrency();

  if (opt.threads <= 0) {

    cerr << "Error: Number of threads cannot be less than or equal to 0" << endl;
    ret = false;
  }
  else if (opt.threads > max_threads) {

    cerr << "Warning: Number of threads cannot be greater than or equal to " << max_threads 
    << ". Setting number of threads to " << max_threads << endl;
    opt.threads = max_threads;
  }

  if (opt.k <= 0) {

    cerr << "Error: Length k of k-mers cannot be less than or equal to 0" << endl;
    ret = false;
  }

  if (opt.k >= MAX_KMER_SIZE) {

    cerr << "Error: Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
    ret = false;
  }

  if (opt.w < opt.k) {

    cerr << "Error: window, " << opt.w << ", cannot be less than the k-mer size, " << opt.k << endl;
    ret = false;
  }


  if (opt.ratio <= 0) {

    cerr << "Error: Ratio of g-mers to conserve cannot be less than or equal to 0" << endl;
    ret = false;
  }

  if (opt.ratio > 1) {

    cerr << "Error: Ratio of g-mers to conserve cannot be more than 1.0" << endl;
    ret = false;
  }

  if (opt.files.size() == 0) {

    cerr << "Error: Missing FASTA/FASTQ input files" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    vector<string>::const_iterator it;
    int intStat;

    for (const auto& it : opt.files) {  
      intStat = stat(it.c_str(), &stFileInfo);

      if (intStat != 0) {
        cerr << "Error: File not found, " << it << endl;
        ret = false;
      }
    }
  }

  return ret;
}

void parse_ProgramOptions_detect(int argc, char **argv, Mindex_opt& opt) {

  const char* opt_string = "t:i:";

  static struct option long_options[] = {

    {"threads",         required_argument,  0, 't'},
    {"index",           required_argument,  0, 'i'},
    {0,                 0,                  0,  0 }
  };

  int option_index = 0, c;

  while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

    switch (c) {

    case 't':
      opt.threads = atoi(optarg);
      break;
    case 'i':
      opt.index = optarg;
      break;
    default:
      break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  while (optind < argc) opt.input.push_back(argv[optind++]);
}

bool check_ProgramOptions_detect(Mindex_opt& opt) {

  bool ret = true;

  size_t max_threads = std::thread::hardware_concurrency();

  if (opt.threads <= 0) {

    cerr << "Error: Number of threads cannot be less than or equal to 0" << endl;
    ret = false;
  }
  else if (opt.threads > max_threads) {

    cerr << "Warning: Number of threads cannot be greater than or equal to " << max_threads 
    << ". Setting number of threads to " << max_threads << endl;
    opt.threads = max_threads;
  }

  if (opt.index.empty()) {
    cerr << "Error: Missing index" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    int intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << "Error: Index is not found " << opt.index << endl;
      ret = false;
    }
  }

  if (opt.input.size() == 0) {

    cerr << "Error: Missing FASTA/FASTQ input files" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    vector<string>::const_iterator it;
    int intStat;

    for (const auto& it : opt.input) {  
      intStat = stat(it.c_str(), &stFileInfo);

      if (intStat != 0) {
        cerr << "Error: File not found, " << it << endl;
        ret = false;
      }
    }
  }

  return ret;
}

void Mindex_Usage() {
  std::cout << "Usage: mindex build [options] fasta-files" << std::endl << std::endl
  << "Options: " << std::endl
  << "-t, --threads         Number of threads to use" << std::endl
  << "-k, --kmer-length     Size of k-mers" << std::endl
  << "-w, --window          Window for minimizers" << std::endl
  << "-r, --ratio           Sampling ratio used" << std::endl << std::endl;
}

int main(int argc, char **argv) {

  if (argc < 2) {
    // Print error message, function?
    Mindex_Usage();
    exit(1);
  } else {

    std::string cmd(argv[1]);
    Mindex_opt opt;


    if (cmd == "build") {
      parse_ProgramOptions_build(argc-1, argv+1, opt);
      if (check_ProgramOptions_build(opt)) { //Program options are valid
        Kmer::set_k(opt.k);
        //Minimizer::set_g(opt.k);

        Mindex mindex(opt.k, opt.w);

        mindex.build(opt);
      }
    } else if (cmd == "detect") {
      parse_ProgramOptions_detect(argc-1,argv+1,opt);
      if (check_ProgramOptions_detect(opt)) {
        Mindex mindex;
        if (mindex.loadFromFile(opt.index, opt)) {
          cerr << "index loaded" << endl;
          // possibly quantify all the data
          cerr << "Minimizers " << mindex.minimizers.size() << endl; 
          mindex.countData(opt);
        }
      }
    }     
  }
}
