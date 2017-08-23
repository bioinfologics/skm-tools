//
// Created by Bernardo Clavijo (EI) on 24/07/2017.
//

#ifndef SKIPMER_SKIPMERMULTIWAYCOVERAGEANALISER_H
#define SKIPMER_SKIPMERMULTIWAYCOVERAGEANALISER_H

#include "SkipMerPosition.h"
#include "SkipMerSpectrum.h"



class SkipMerMultiWayCoverageAnalyser {
public:
    explicit SkipMerMultiWayCoverageAnalyser(SkipMerPosition & _skmp, uint8_t maxcoverage) : ref(_skmp) {
        count_reference_coverage(maxcoverage);
    };
    void count_reference_coverage(uint8_t maxcoverage);
    void add_coverage_track(SkipMerSpectrum &skms, uint8_t maxcoverage);
    void dump_features_conservation_scores(std::string & filename);
    void dump_all_tracks(std::string & filename);
    void dump_features_tracks_stats(std::string & filename);
    void project_max_coverage_to_bases();
    std::vector<uint64_t> anchor_stats();
    std::vector<uint64_t> covered_by_all_stats();
    std::vector<uint64_t> covered_by_all_projected_stats();
    SkipMerPosition ref;
    std::vector<std::vector<uint8_t>> coverage_tracks;
};


#endif //SKIPMER_SKIPMERMULTIWAYCOVERAGEANALISER_H
