#include <vector>

#include <mlalign.h>
#include <celero/Celero.h>
#include <fstream>
#include <assert.h>
#include <iostream>
#include <chrono>

CELERO_MAIN

struct Sequence {
  std::string description;
  std::vector<mlalign_Nucleotide> bases;
};

std::vector<Sequence> parse_FASTA(const std::string& filename);

inline void benchmark(const std::vector<Sequence>& sequences);


#define BENCHMARK_PROTEIN(protein, baseline) \
  BASELINE(protein, Baseline, 1, 1) { std::this_thread::sleep_for(std::chrono::microseconds(baseline)); } \
  const std::vector<Sequence> protein = parse_FASTA("sequences/" #protein "_unaligned.fas"); \
  BENCHMARK(protein, Impl, 5, 10) { benchmark(protein); }

BENCHMARK_PROTEIN(BDNF, 126494);
BENCHMARK_PROTEIN(cytb, 419193);
BENCHMARK_PROTEIN(RAG1, 355899);
BENCHMARK_PROTEIN(RAG2, 666819);
BENCHMARK_PROTEIN(RBP3, 720429);
BENCHMARK_PROTEIN(vWF, 953576);

#undef BENCHMARK_PROTEIN


std::vector<Sequence> parse_FASTA(const std::string& filename) {
  std::vector<Sequence> sequences;
  std::ifstream f(filename);

  if (!f.is_open()) {
    std::cerr << "ERROR: Could not open file " << filename << std::endl;
    assert(0);
  }

  // Would love to use std::regex, but GCC.

  Sequence cur = Sequence();
  std::string line;
  bool first = true;

  while (std::getline(f, line)) {

    if (line[0] == '>') {

      if (!first) {
        sequences.push_back(cur);
      }

      first = false;
      cur.description = line.substr(1, line.length() - 1);
      cur.bases.clear();

    } else {

      for (const char c : line) {
        // essentially what we did in main.cpp
        switch (c) {
          case 'A':
          case 'N': // Actually this is a wildcard, but w/e
          case 'R': // A or G
          case 'M': // A or C
            cur.bases.push_back(mlalign_A);
            break;
          case 'C':
          case 'Y': // C or T
          case 'S': // C or G
          case 'B': // not A
            cur.bases.push_back(mlalign_C);
            break;
          case 'G':
          case 'K': // G or T
          case 'D': // not C
          case 'V': // neither T nor U
            cur.bases.push_back(mlalign_G);
            break;
          case 'T':
          case 'U': // Uracil
          case 'W': // A or T
          case 'H': // not G
            cur.bases.push_back(mlalign_T);
            break;
          default:
            std::cerr << "Invalid character '" << c << "'" << std::endl;
            assert(0);
            break;
        }

      }

    }
  }

  if (!first) {
    sequences.push_back(cur);
  }

  return sequences;
}

inline void align_non_uniform(const std::vector<mlalign_Nucleotide>& as,
                              const std::vector<mlalign_Nucleotide>& bs,
                              std::vector<mlalign_Site>& result) {
  double pi[4] = {0.27, 0.24, 0.26, 0.23 };
  result.clear();
  result.resize(as.size() + bs.size());
  size_t size = result.size();
  double score = 0.0;
  mlalign_TKF91(const_cast<mlalign_Nucleotide*>(as.data()), as.size(),
                const_cast<mlalign_Nucleotide*>(bs.data()), bs.size(),
                1.0, 2.0, 0.1, pi, result.data(), &size, &score);
}

inline void benchmark(const std::vector<Sequence>& sequences) {
    const size_t N = 10;
    // Pick out N sequences and align them pairwise (= (N*(N-1))/2 Alignments)
    size_t stride = sequences.size() / N;
    std::vector<mlalign_Site> alignment; 
 
    for (size_t i = 0; i < N; ++i) { 
      for (size_t j = i + 1; j < N; ++j) { 
        auto& a = sequences[i*stride].bases;
        auto& b = sequences[j*stride].bases;
 
        align_non_uniform(a, b, alignment); 
      } 
    } 
}