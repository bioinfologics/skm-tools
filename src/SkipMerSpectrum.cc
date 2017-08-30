//
// Created by Bernardo Clavijo (TGAC) on 27/06/2017.
//

#include <iostream>
#include <fstream>
#include <parallel/algorithm>
#include <fcntl.h>

#include "SkipMerSpectrum.h"
#include "deps/zstr/src/zstr.hpp"

void SkipMerSpectrum::count_from_file(std::string filename, bool fasta_input){
    gzFile fp;

    uint64_t curr_offset=0;

    int l = 0;

    kseq fasta_entry;
    if ( (filename.substr(filename.find_last_of('.') + 1) == "gz") or (filename.substr(filename.find_last_of('.') + 1) == "bz2")) {
#ifdef HAVE_ZLIB
        if (filename.substr(filename.find_last_of('.') + 1) == "gz"){
            gzFile fp = gzopen(filename.c_str(), "r");
            FunctorZlib gzr;
            kstream<gzFile, FunctorZlib> ks(fp, gzr);
            while((l = ks.read(fasta_entry)) >= 0) {
                if (fasta_entry.seq.size() < S) continue;
                add_from_string(fasta_entry.seq);
            }
            gzclose(fp);
        }
#endif
#ifdef HAVE_BZIP2
        if (filename.substr(filename.find_last_of('.') + 1) == "bz2") {
            BZFILE * fp = BZ2_bzopen(filename.c_str(), "r");
            FunctorBZlib2 bzr;
            kstream<BZFILE*, FunctorBZlib2>ks(fp, bzr);
            while((l = ks.read(fasta_entry)) >= 0) {
                if (fasta_entry.seq.size() < S) continue;
                add_from_string(fasta_entry.seq);
            }
            BZ2_bzclose(fp);
        }
#endif
    } else {
        int fp = open(filename.c_str(), O_RDONLY);
        FunctorRead r;
        kstream<int, FunctorRead> ks(fp, r);
        while((l = ks.read(fasta_entry)) >= 0) {
            if (fasta_entry.seq.size() < S) continue;
            add_from_string(fasta_entry.seq);
        }
        close(fp);
    }

    std::cout<<skipmers.size()<<" skipmers. ";
};

void SkipMerSpectrum::add_from_string(const std::string & seq){
    const SkipMer KMER_MASK=( ((SkipMer)1)<<(k*2) )-1;
    const SkipMer KMER_FIRSTOFFSET=(k-1)*2;
    const char * s=seq.c_str();
    int64_t last_unknown[n];
    SkipMer fkmer[n];
    SkipMer rkmer[n];
    uint8_t cycle_pos[n];
    uint8_t fi;

    for (auto i=0;i<n;++i){
        cycle_pos[i]=n-i-1;
        last_unknown[i]=-1;
        fkmer[i]=0;
        rkmer[i]=0;
    }
    if (skipmers.capacity()<skipmers.size()+seq.size()) skipmers.reserve(skipmers.size()+(seq.size()>alloc_block? seq.size():alloc_block));

    for (int64_t p=0;p<seq.size();p++) {
        //fkmer: grows from the right (LSB)
        //rkmer: grows from the left (MSB)
        for (auto ni=0;ni<n;++ni) {

            cycle_pos[ni]++;
            if (cycle_pos[ni]==n)cycle_pos[ni]=0;

            if (cycle_pos[ni]<m){
                switch (s[p]) {
                    case 'A':
                    case 'a':
                        fkmer[ni] = ((fkmer[ni] << 2) + 0) & KMER_MASK;
                        rkmer[ni] = (rkmer[ni] >> 2) + (((SkipMer) 3) << KMER_FIRSTOFFSET);
                        break;
                    case 'C':
                    case 'c':
                        fkmer[ni] = ((fkmer[ni] << 2) + 1) & KMER_MASK;
                        rkmer[ni] = (rkmer[ni] >> 2) + (((SkipMer) 2) << KMER_FIRSTOFFSET);
                        break;
                    case 'G':
                    case 'g':
                        fkmer[ni] = ((fkmer[ni] << 2) + 2) & KMER_MASK;
                        rkmer[ni] = (rkmer[ni] >> 2) + (((SkipMer) 1) << KMER_FIRSTOFFSET);
                        break;
                    case 'T':
                    case 't':
                        fkmer[ni] = ((fkmer[ni] << 2) + 3) & KMER_MASK;
                        rkmer[ni] = (rkmer[ni] >> 2) + (((SkipMer) 0) << KMER_FIRSTOFFSET);
                        break;
                    default:
                        fkmer[ni] = ((fkmer[ni] << 2) + 0) & KMER_MASK;
                        rkmer[ni] = (rkmer[ni] >> 2) + (((SkipMer) 3) << KMER_FIRSTOFFSET);
                        last_unknown[ni] = p;
                        break;
                }
            }
        }
        //if we are at p, the skip-mer started at p-S is now done
        if (p>=S-1){
            if (p==S-1) fi=0;
            else {
                ++fi;
                if (fi==n) fi=0;
            }
            if (last_unknown[fi] + S <= p ) {
                if (fkmer[fi] <= rkmer[fi]) {
                    skipmers.emplace_back(fkmer[fi],1);
                } else {
                    skipmers.emplace_back(rkmer[fi],1);
                }
            }
        }
    }
};

void SkipMerSpectrum::sort_and_collapse(){
    __gnu_parallel::sort(skipmers.begin(),skipmers.end());

    size_t wi=0;
    for (size_t ri=1; ri<skipmers.size();++ri){
        if (skipmers[ri].skipmer==skipmers[wi].skipmer) {
            uint64_t new_count=skipmers[wi].count+skipmers[ri].count;
            skipmers[wi].count=(uint16_t)(new_count < UINT16_MAX ? new_count : UINT16_MAX);
        }
        else skipmers[++wi]=skipmers[ri];
    }
    skipmers.resize(wi+1);
    std::cout << skipmers.size() << " unique. ";
}

void SkipMerSpectrum::sort_and_collapse(uint16_t min_freq,uint16_t max_freq){
    uint64_t pass=0,fail=0;
    __gnu_parallel::sort(skipmers.begin(),skipmers.end());
    size_t wi=0;
    for (size_t ri=1; ri<skipmers.size();++ri){
        if (skipmers[ri].skipmer==skipmers[wi].skipmer) {
            uint64_t new_count=skipmers[wi].count+skipmers[ri].count;
            skipmers[wi].count=(uint16_t)(new_count < UINT16_MAX ? new_count : UINT16_MAX);
        }
        else {
            if (skipmers[wi].count>=min_freq and skipmers[wi].count<=max_freq){
                ++wi;
                ++pass;
            } else {
                ++fail;
            }
            skipmers[wi]=skipmers[ri];
        }
    }
    if (skipmers[wi].count>=min_freq and skipmers[wi].count<=max_freq){
        ++wi;
        ++pass;
    } else {
        ++fail;
    }
    skipmers.resize(wi);
    std::cout<<"Filtering "<< pass+fail<<" skip-mers by frequency in ["<<min_freq<<"-"<<max_freq<<"], pass="
            <<pass<<", fail="<<fail<<std::endl;
}

void SkipMerSpectrum::dump(std::string filename){
    zstr::ofstream output_file(filename+".skm_spectrum.gz");
    output_file.write((char *)&m, sizeof(m));
    output_file.write((char *)&n, sizeof(n));
    output_file.write((char *)&k, sizeof(n));
    uint64_t c=skipmers.size();
    output_file.write((char *)&c, sizeof(c));
    output_file.write((char *)skipmers.data(), c*sizeof(SkipMerEntry));
}; //dumps m,n,k, count and then just the binary skipmer

void SkipMerSpectrum::generate_stats(std::string filename){
    // First part: create the histogram
    size_t u=0,m=0, tm=0;
    uint64_t hist[100001]={0};
    for (auto skm:skipmers) {
        ++hist[std::min((int)skm.count,(int)100000)];
        if (skm.count==1) ++u;
        else {
            ++m;
            tm+=skm.count;
        }
    }
    std::cout << "Total  skipmers:       " << u+tm << std::endl;
    std::cout << "Distinct skipmers:     " << u+m << std::endl;
    std::cout << "Single-copy skipmers:  " << u << std::endl;
    std::cout << "Multi-copy skipmers:   " << m << std::endl;
    std::ofstream hist_file(filename+"_coveragehist.csv");
    hist_file<<"Frequency,Count"<<std::endl;
    std::ofstream mdist_file(filename+"_maxdistinct.hist");
    mdist_file<<"Total,Distinct"<<std::endl<<"0,0"<<std::endl;
    uint64_t pct_count=(u+tm)/100;
    uint64_t next_pct_count=pct_count;
    uint64_t count=0,tcount=0;
    for (auto i=0;i<100001;++i){
        if (0!=hist[i]) hist_file<<i<<","<<hist[i]<<std::endl;
        while (tcount+i*hist[i]>=next_pct_count){
            mdist_file<<next_pct_count<<","<<count+(next_pct_count-tcount)/i<<std::endl;
            next_pct_count+=pct_count;
        }
        tcount+=i*hist[i];
        count+=hist[i];
    }
    if (next_pct_count==100*pct_count) mdist_file<<tcount<<","<<count<<std::endl;

    hist_file.close();
    if (0!=hist[100000]) std::cout<<"WARNING: skip-mers with coverage larger than 100,000 will be imprecisely reported"<<std::endl;

}; //writes both a .txt stats file and a histogram