#include <cassert>
#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
#include <chrono>
#include "optionparser.h"

#include "mlalign.h"

static std::vector<mlalign_Nucleotide> parseNucleotideString(const std::string &ns) {
  std::vector<mlalign_Nucleotide> ret;
  ret.reserve(ns.size());
  std::transform(ns.begin(), ns.end(), std::back_inserter(ret), [](char c) {
    switch (tolower(c)) {
      case 'a':
        return mlalign_A;
      case 'c':
        return mlalign_C;
      case 'g':
        return mlalign_G;
      case 't':
        return mlalign_T;
      default:
        // This doesn't happen because we already check the input when parsing the args.
        throw std::domain_error("The nucleotide sequence was in a bad format.");
    }
  });

  return ret;
}

static char toChar(mlalign_Nucleotide n) {
  switch (n) {
    case mlalign_A:
      return 'A';
    case mlalign_C:
      return 'C';
    case mlalign_G:
      return 'G';
    case mlalign_T:
      return 'T';
    case mlalign_Gap:
      return '-';
    default:
      // This would mean that the output of the algorithm was corrupt.
      throw std::domain_error("Invalid mlalign_Nucleotide value.");
  }
}

enum OptionIndex {
  UNKNOWN, HELP, SEQUENCE1, SEQUENCE2, LAMBDA, MU, TAU, PI_ADENINE, PI_CYTOSINE, PI_GUANINE, PI_THYMINE, RUNS
};

const option::Descriptor usage[] = {
        {UNKNOWN,     0, "",  "",           option::Arg::None,             "USAGE: mlalign [options]\n\n"
                                                                                   "mlalign - Compute a maximum likelihood alignment between two sequences.\n\n"
                                                                                   "Options:"},
        {HELP,        0, "h", "help",       option::Arg::None,             "  --help, -h  \tPrint usage and exit."},
        {SEQUENCE1,   0, "a", "sequence-1", option::Arg::NucleotideString, "  --sequence-1, -a  \tFirst nucleotide sequence, string encoded."
                                                                                   "                   Format: (A|C|G|T)*"},
        {SEQUENCE2,   0, "b", "sequence-2", option::Arg::NucleotideString, "  --sequence-2, -b  \tSecond nucleotide sequence, string encoded."
                                                                                   "                  Format: (A|C|G|T)*"},
        {LAMBDA,      0, "l", "lambda",     option::Arg::Float,            "  --lambda, -l  \tBirth rate parameter"},
        {MU,          0, "m", "mu",         option::Arg::Float,            "  --mu, -mu  \tDeath rate parameter"},
        {TAU,         0, "t", "tau",        option::Arg::Float,            "  --tau, -t  \tTime delta"},
        {PI_ADENINE,  0, "A", "pa",         option::Arg::Float,            "  --pa, -A  \tEquilibrium frequency for adenine"},
        {PI_CYTOSINE, 0, "C", "pc",         option::Arg::Float,            "  --pc, -C  \tEquilibrium frequency for cytosine"},
        {PI_GUANINE,  0, "G", "pg",         option::Arg::Float,            "  --pg, -G  \tEquilibrium frequency for guanine"},
        {PI_THYMINE,  0, "T", "pt",         option::Arg::Float,            "  --pt, -T  \tEquilibrium frequency for thymine"},
        {RUNS,        0, "r", "runs",       option::Arg::Numeric,          "  --runs, -r  \tNumber of alignment runs"},
        {0,           0, "",  "",           0,                             "\nExample:\n"
                                                                                   "  mlalign --sequence-1 ACGATA --sequence-2 AGTGGTA --lambda 1 "
                                                                                   "--mu 2 --tau 0.1 --pa 0.25 --pc 0.25 --pg 0.25 --pt 0.25\n"},
        {0,           0, 0,   0,            0,                             0}
};

#define REQUIRED_OPTION(idx) \
   if (!options[idx] || options[idx].arg == nullptr) { \
     std::cerr << "Option '--" << usage[idx].longopt << "' is required." << std::endl; \
     return 1; \
   }

int main(int argc, char **argv) {
  argc -= (argc > 0);
  argv += (argc > 0); // skip program name argv[0] if present
  option::Stats stats(usage, argc, argv);
  std::vector<option::Option> options(stats.options_max), buffer(stats.buffer_max);
  option::Parser parse(usage, argc, argv, options.data(), buffer.data());

  if (parse.error()) {
    return 1;
  }

  if (argc == 0 || options[HELP]) {
    option::printUsage(std::cout, usage);
    return 0;
  }

  REQUIRED_OPTION(SEQUENCE1);
  REQUIRED_OPTION(SEQUENCE2);
  REQUIRED_OPTION(LAMBDA);
  REQUIRED_OPTION(MU);
  REQUIRED_OPTION(TAU);
  REQUIRED_OPTION(PI_ADENINE);
  REQUIRED_OPTION(PI_CYTOSINE);
  REQUIRED_OPTION(PI_GUANINE);
  REQUIRED_OPTION(PI_THYMINE);

  int runs = (!options[RUNS] || options[RUNS].arg == nullptr) ? 1 : atoi(options[RUNS].arg);

  std::vector<mlalign_Nucleotide> A, B;
  std::string sequence1(options[SEQUENCE1].arg);
  std::string sequence2(options[SEQUENCE2].arg);

  A = parseNucleotideString(sequence1);
  B = parseNucleotideString(sequence2);

  double lambda = atof(options[LAMBDA].arg);
  double mu = atof(options[MU].arg);
  double tau = atof(options[TAU].arg);
  std::array<double, 4> pi;
  pi[mlalign_A] = atof(options[PI_ADENINE].arg);
  pi[mlalign_C] = atof(options[PI_CYTOSINE].arg);
  pi[mlalign_G] = atof(options[PI_GUANINE].arg);
  pi[mlalign_T] = atof(options[PI_THYMINE].arg);

  std::vector<mlalign_Site> alignment(A.size() + B.size());
  size_t n;
  double score;

  for (int i = 0; i < runs; i++) {
    auto startTime = std::chrono::system_clock::now();
    mlalign_TKF91(A.data(), A.size(), B.data(), B.size(), lambda, mu, tau, pi.data(), alignment.data(), &n, &score);
    alignment.resize(n);
    auto endTime = std::chrono::system_clock::now();

    std::cout << "Run time: " << std::chrono::duration_cast<std::chrono::microseconds>(endTime-startTime).count() << " µs" << std::endl;

    std::cout << "Found alignment of length " << n << " with score " << score << ":" << std::endl;
    for (auto site : alignment) {
      std::cout << toChar(site.a);
    }
    std::cout << std::endl;
    for (auto site : alignment) {
      std::cout << toChar(site.b);
    }
    std::cout << std::endl;
  }

  return 0;
}

