//
// Created by Bernardo Clavijo (TGAC) on 27/06/2017.
//

#ifndef SKIPMERSPECTRA_H
#define SKIPMERSPECTRA_H

#include <vector>
#include <string>
#include "deps/kseqcpp/kseq.hpp"

#include "SkipMer.h"

class SkipMerEntry{
public:
    SkipMer skipmer;
    uint16_t count;
    bool operator<(const SkipMerEntry & rhs) const { return skipmer<rhs.skipmer;};
};

class SkipMerSpectrum {
public:
    SkipMerSpectrum(uint8_t _m, uint8_t _n, uint8_t _k, uint64_t _alloc_block=10000000) : m(_m), n(_n), k(_k), alloc_block(_alloc_block) {
        S=_k;
        S=S+((S-1)/(int)_m)*((int)_n-(int)_m);
    };
    void count_from_file(std::string filename, bool fasta_input);
    void add_from_string(const std::string & s);
    void sort_and_collapse();
    void sort_and_collapse(uint16_t min_freq, uint16_t max_freq);
    void dump(std::string filename); //dumps m,n,k, count and then just the binary skipmer
    void generate_stats(std::string filename); //writes both a .txt stats file and a histogram
    uint8_t m,n,k;
    int S;
    std::vector<SkipMerEntry> skipmers;
    uint64_t count=0;
    uint64_t alloc_block;
};


#endif //SKIPMERSPECTRA_H
