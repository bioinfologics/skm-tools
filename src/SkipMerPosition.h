//
// Created by Bernardo Clavijo (TGAC) on 27/06/2017.
//

#ifndef SKIPMERPOSITION_H
#define SKIPMERPOSITION_H

#include <vector>
#include <string>
#include "deps/kseqcpp/kseq.hpp"

#include "SkipMer.h"

class SkipMerPositionEntry{
public:
    SkipMer skipmer;
    uint64_t position;
    bool rev;
    const bool operator<(const SkipMerPositionEntry & rhs) const { return skipmer<rhs.skipmer;};
};

class ReferenceSeq{
public:
    std::string name;
    uint64_t offset;
    uint64_t size;
};

class ReferenceFeature{
public:
    ReferenceFeature(std::string & n,uint64_t s,uint64_t e ) : name(n),start(s),end(e){};
    std::string name;
    uint64_t start;
    uint64_t end;
    std::vector<uint32_t> votes_from_refs;
};

class SkipMerPosition {
public:
    uint8_t m,n,k;
    int S;
    SkipMerPosition(uint8_t _m, uint8_t _n, uint8_t _k, uint64_t _alloc_block=10000000) : m(_m), n(_n), k(_k), alloc_block(_alloc_block) {
        S=_k;
        S=S+((S-1)/(int)_m)*((int)_n-(int)_m);
    };
    void create_from_reference_file(std::string filename, bool fasta_input);
    void mark_with_gff3_feature(std::string filename, std::string feature, uint64_t min_feature_size=1000);
    void add_from_string(const std::string & s,uint64_t offset);
    void sort();
    void filter_repeated();
    void dump(std::string filename); //dumps m,n,k, count and then just the binary skipmer
    void load(std::string filename);
    void generate_stats(std::string filename); //writes both a .txt stats file and a histogram
    void merge();
    void intersect();//kat comp's equivalent
    std::vector<SkipMerPositionEntry> skipmer_positions;
    std::vector<ReferenceSeq> refseqs;
    std::vector<ReferenceFeature> reffeatures;
    std::vector<bool> in_feature;
    std::vector<uint8_t> shared_count;

private:

    uint64_t alloc_block;
};


#endif //SKIPMERPOSITION_H
