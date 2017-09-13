//
// Created by Bernardo Clavijo (EI) on 24/07/2017.
//

#include "SkipMerMultiWayCoverageAnalyser.h"
#include <iostream>

#if defined(__clang__)
#include <algorithm>
#elif defined(__GNUC__) || !defined(__GNUG__)
#include <parallel/algorithm>
#endif

#include <fstream>
#include "deps/zstr/src/zstr.hpp"


void SkipMerMultiWayCoverageAnalyser::add_coverage_track(SkipMerSpectrum &skms,uint8_t maxcoverage) {
    coverage_tracks.emplace_back(ref.refseqs.back().offset+ref.refseqs.back().size,0);
    auto current_pos=ref.skipmer_positions.begin();
    for (auto &c:skms.skipmers){
        while (current_pos<ref.skipmer_positions.end() and current_pos->skipmer < c.skipmer) ++current_pos;
        for (;current_pos<ref.skipmer_positions.end() and current_pos->skipmer == c.skipmer;++current_pos)
            if (coverage_tracks[0][current_pos->position]>0) coverage_tracks.back()[current_pos->position]=(c.count<255 ? c.count:255);
    }
    for (auto &c:coverage_tracks.back()){
        if (maxcoverage != 0 && c>maxcoverage) c=0;
    }
}


void SkipMerMultiWayCoverageAnalyser::count_reference_coverage(uint8_t maxcoverage) {
    coverage_tracks.emplace_back(std::vector<uint8_t>(ref.refseqs.back().offset+ref.refseqs.back().size,0));
    for (auto fis=ref.skipmer_positions.begin();fis<ref.skipmer_positions.end();){
        auto lis=fis;
        uint64_t cv=0;
        while (lis->skipmer==fis->skipmer) {
            ++lis;
            ++cv;
        }
        lis=fis;
        while (lis->skipmer==fis->skipmer) {
            coverage_tracks.back()[lis->position]=(cv<255? cv:255);
            ++lis;
        }
        fis=lis;
    }
    for (auto &c:coverage_tracks.back()){
        if (c>maxcoverage) c=0;
    }
}

void SkipMerMultiWayCoverageAnalyser::dump_features_conservation_scores(std::string &filename) {
    auto fn=filename+"_detailed_coverage.csv";
    //dump_all_tracks(fn);
    fn=filename+"_detailed_feat.csv";
    dump_features_tracks_stats(fn);
    project_max_coverage_to_bases();

    fn=filename+"_projected_coverage.csv";
    //dump_all_tracks(fn);
    fn=filename+"_projected_feat.csv";
    dump_features_tracks_stats(fn);
}

void SkipMerMultiWayCoverageAnalyser::dump_all_tracks(std::string &filename){
    std::ofstream detcov(filename);

    detcov<<"Ref";
    for (auto i=1;i<coverage_tracks.size();++i) detcov<<",Set"<<i;
    detcov<<std::endl;
    for (uint64_t i=0;i<coverage_tracks[0].size();++i){
        detcov<<(int)coverage_tracks[0][i];
        for (auto j=1;j<coverage_tracks.size();++j) detcov<<","<<(int)coverage_tracks[j][i];
        detcov<<std::endl;
    }
}

void SkipMerMultiWayCoverageAnalyser::dump_features_tracks_stats(std::string &filename){
    std::ofstream detcov(filename);
    size_t ct=coverage_tracks.size()-1;
    detcov<<"Feature,Ref_0,Ref_1,Ref_2,Ref_M";
    for (auto i=1;i<coverage_tracks.size();++i) detcov<<",Set"<<i<<"_0,Set"<<i<<"_1,Set"<<i<<"_2,Set"<<i<<"_M,Set"<<i<<"P";
    detcov<<",Any,PAny,All,PAll"<<std::endl;
    for (auto &f:ref.reffeatures){
        detcov<<f.name;
        bool refcounts=true;
        uint64_t denom=0;
        for (auto &c:coverage_tracks) {
            uint64_t f0 = 0, f1 = 0, f2 = 0, fm = 0;
            for (auto p = f.start; p < f.end; ++p) {
                if (0==c[p]) ++f0;
                else if (1==c[p]) ++f1;
                else if (2==c[p]) ++f2;
                else ++fm;
            }
            detcov<<","<<f0<<","<<f1<<","<<f2<<","<<fm;
            if (refcounts){
                denom=f1+f2+fm;
                refcounts=false;
            } else detcov<<","<<(denom ? ((float) (f1+f2+fm)/ (float)denom) : 0);
        }
        //Stupid way, it goes again through all the sets, but I can't be bothered now.
        uint64_t any=0,all=0;
        for (auto p = f.start; p < f.end; ++p) {
            unsigned t=0;

            for (auto i=1;i<=ct;++i) if (coverage_tracks[i][p]) ++t;
            if (t) {
                ++any;
                if (t==ct) ++all;
            }
        }
        detcov<<","<<any<<","<<(denom ? ((float) (any)/ (float)denom) : 0)<<","<<all<<","<<(denom ? ((float) (all)/ (float)denom) : 0);
        detcov<<std::endl;
    }
}

void SkipMerMultiWayCoverageAnalyser::project_max_coverage_to_bases(){
    //in-place, goes from the back because skip-mers go fwd.
    uint64_t p;
    int cp,t;
    for (auto &c:coverage_tracks){
        for (int64_t i=c.size()-ref.S-1;i>=0;--i){
            if (c[i]){
                for (p=i, cp=0,t=0;t<ref.S;++p,cp=(cp == ref.n-1 ? 0 : cp+1),++t){
                    if (cp<ref.m and c[i]>c[p]) c[p]=c[i];
                }
            }
        }
    }
}

std::vector<uint64_t> SkipMerMultiWayCoverageAnalyser::covered_by_all_stats(){
    //Total Ref,Total Any,Total All,Feature Ref,Feature Any,Feature All
    std::vector<bool> in_feature_skm(coverage_tracks[0].size());
    for (auto &f:ref.reffeatures){
        for (uint64_t p=f.start;p<=f.end-ref.S;++p) in_feature_skm[p]=true;
    }

    uint64_t total_ref=0,total_any=0,total_all=0,feature_ref=0,feature_any=0,feature_all=0;

    for (uint64_t p=0;p<coverage_tracks[0].size();++p) {
        if (coverage_tracks[0][p]>0) {
            ++total_ref;
            if (in_feature_skm[p]) ++feature_ref;
            uint64_t sc=0;
            for (auto si=1;si<coverage_tracks.size();++si) if (coverage_tracks[si][p]>0) ++sc;
            if (sc>0){
                ++total_any;
                if (in_feature_skm[p]) ++feature_any;
                if (sc==coverage_tracks.size()-1) {
                    ++total_all;
                    if (in_feature_skm[p]) ++feature_all;
                }
            }
        }
    }
    return {total_ref,total_any,total_all,feature_ref,feature_any,feature_all};
}

std::vector<uint64_t> SkipMerMultiWayCoverageAnalyser::covered_by_all_projected_stats(){
    //Total Ref,Total Any,Total All,Feature Ref,Feature Any,Feature All
    std::vector<bool> projected_ref(coverage_tracks[0].size(),false);
    std::vector<bool> projected_all(coverage_tracks[0].size(),false);
    std::vector<bool> projected_any(coverage_tracks[0].size(),false);
    std::vector<bool> in_feature(coverage_tracks[0].size(),false);

    for (auto &f:ref.reffeatures){
        for (uint64_t p=f.start;p<=f.end;++p) in_feature[p]=true;
    }

    uint64_t total_ref=0,total_any=0,total_all=0,feature_ref=0,feature_any=0,feature_all=0;
    uint64_t xp,xcp,xt;

    for (uint64_t p=0;p<coverage_tracks[0].size()-ref.S;++p) {
        if (coverage_tracks[0][p]>0) {
            for (xp=p, xcp=0,xt=0;xt<ref.S;++xp,xcp=(xcp == ref.n-1) ? 0 : xcp+1,++xt)
                if (!projected_ref[xp] && xcp<ref.m){
                    ++total_ref;
                    if (in_feature[xp]) ++feature_ref;
                    projected_ref[xp]=true;
                }


            uint64_t sc=0;//count of datasets that HAVE the current skip-mer
            for (auto si=1;si<coverage_tracks.size();++si) if (coverage_tracks[si][p]>0) ++sc;
            if (sc>0){
                for (xp=p, xcp=0,xt=0;xt<ref.S;++xp,xcp=(xcp == ref.n-1) ? 0 : xcp+1,++xt)
                    if (!projected_any[xp] and xcp<ref.m){
                        ++total_any;
                        if (in_feature[p]) ++feature_any;
                        projected_any[xp]=true;
                    }

                if (sc==coverage_tracks.size()-1)
                    for (xp=p, xcp=0,xt=0;xt<ref.S;++xp,xcp=(xcp == ref.n-1) ? 0 : xcp+1,++xt)
                        if (!projected_all[xp] and xcp<ref.m){
                            ++total_all;
                            if (in_feature[p]) ++feature_all;
                            projected_all[xp]=true;
                        };

            }
        }
    }
    return {total_ref,total_any,total_all,feature_ref,feature_any,feature_all};
}