//
// Created by Bernardo Clavijo (TGAC) on 27/06/2017.
//

#include <iostream>
#include <fstream>
#include <parallel/algorithm>
#include <map>
#include <fcntl.h>
#include "SkipMerPosition.h"
#include "deps/zstr/src/zstr.hpp"

void SkipMerPosition::create_from_reference_file(std::string filename, bool fasta_input) {
    std::ifstream seqfile(filename.c_str());

    uint64_t curr_offset=0;

    int l = 0;

    kseq fasta_entry;
    if ( (filename.substr(filename.find_last_of('.') + 1) == "gz") or (filename.substr(filename.find_last_of('.') + 1) == "bz2")) {
#ifdef HAVE_ZLIB
        if (filename.substr(filename.find_last_of('.') + 1) == "gz") {
            gzFile fp = gzopen(filename.c_str(), "r");
            FunctorZlib gzr;
            kstream<gzFile, FunctorZlib> ks(fp, gzr);
            while ((l = ks.read(fasta_entry)) >= 0) {

                if (fasta_entry.seq.size() < S) continue;
                ReferenceSeq rs;
                rs.name = fasta_entry.name;
                rs.offset = curr_offset;
                rs.size = fasta_entry.seq.size();
                add_from_string(fasta_entry.seq, curr_offset);
                refseqs.push_back(rs);
                curr_offset += fasta_entry.seq.size();
            }
            gzclose(fp);
        }
#endif
#ifdef HAVE_BZIP2
        if (filename.substr(filename.find_last_of('.') + 1) == "bz2") {
            BZFILE *fp = BZ2_bzopen(filename.c_str(), "r");
            FunctorBZlib2 bzr;
            kstream<BZFILE *, FunctorBZlib2> ks(fp, bzr);
            while ((l = ks.read(fasta_entry)) >= 0) {

                if (fasta_entry.seq.size() < S) continue;
                ReferenceSeq rs;
                rs.name = fasta_entry.name;
                rs.offset = curr_offset;
                rs.size = fasta_entry.seq.size();
                add_from_string(fasta_entry.seq, curr_offset);
                refseqs.push_back(rs);
                curr_offset += fasta_entry.seq.size();
            }
            BZ2_bzclose(fp);
        }
#endif
    } else {
        int fp = open(filename.c_str(), O_RDONLY);
        FunctorRead r;
        kstream<int, FunctorRead> ks(fp, r);
        while ((l = ks.read(fasta_entry)) >= 0) {
            if (fasta_entry.seq.size() < S) continue;
            ReferenceSeq rs;
            rs.name = fasta_entry.name;
            rs.offset = curr_offset;
            rs.size = fasta_entry.seq.size();
            add_from_string(fasta_entry.seq, curr_offset);
            refseqs.push_back(rs);
            curr_offset += fasta_entry.seq.size();
        }

        close(fp);
    }

    shared_count.resize(refseqs.back().offset+refseqs.back().size);
    //std::cout<<"Counting finished with "<<skipmers.size()<<" skipmers"<<std::endl;
};

void SkipMerPosition::add_from_string(const std::string & seq, uint64_t offset){
    const SkipMer KMER_MASK=( ((SkipMer)1)<<(k*2) )-1;
    const SkipMer KMER_FIRSTOFFSET=(k-1)*2;
    const char * s=seq.c_str();
    int64_t last_unknown[n];
    SkipMer fkmer[n];
    SkipMer rkmer[n];
    uint8_t cycle_pos[n];
    uint8_t fi(0);

    for (auto i=0;i<n;++i){
        cycle_pos[i]=n-i-1;
        last_unknown[i]=-1;
        fkmer[i]=0;
        rkmer[i]=0;
    }

    if (skipmer_positions.capacity()<skipmer_positions.size()+seq.size()) skipmer_positions.reserve(skipmer_positions.size()+(seq.size()>alloc_block? seq.size():alloc_block));


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
        if (p>=S-1) {
            if (p == S - 1) fi = 0;
            else {
                ++fi;
                if (fi == n) fi = 0;
            }
            if (last_unknown[fi] + S <= p) {
                if (fkmer[fi] <= rkmer[fi]) {
                    skipmer_positions.emplace_back(fkmer[fi],offset + p,false);
                } else {
                    skipmer_positions.emplace_back(rkmer[fi],offset + p,true);
                }
            }
        }
    }
};

void SkipMerPosition::sort(){
    __gnu_parallel::sort(skipmer_positions.begin(),skipmer_positions.end());
}

void SkipMerPosition::mark_with_gff3_feature(std::string filename, std::string featurename, uint64_t min_feature_size){
    std::cout<<"Parsing GFF3 on "<<filename<<" for features of type '"<<featurename<<"'"<<std::endl;
    std::ifstream gff3_file(filename.c_str());
    std::map<std::string,ReferenceSeq> seqsmap;
    for (auto &s:refseqs) seqsmap[s.name]=s;
    in_feature.resize(refseqs.back().offset+refseqs.back().size);
    //read a line, check it is the relevant feature
    std::string line,sname,dummy,fname;
    uint64_t start,end;
    uint64_t all_lines=0,noseq=0,toosmall=0,used=0,tooshort=0;
    while (!gff3_file.eof()) {
        std::getline(gff3_file, line);
        if (line.empty()) continue;
        ++all_lines;
        std::stringstream lss(line);
        lss>>sname;
        lss>>dummy;
        lss>>fname;
        if (fname!=featurename) continue;
        lss>>start;
        lss>>end;
        auto s=seqsmap.find(sname);
        if (s==seqsmap.end()) {
            ++noseq;
            continue;
        }
        if (end<start) {
            auto t=start;
            start=end;
            end=t;
        }
        if (end-start<min_feature_size) {
            ++tooshort;
            continue;
        }
        if (end>s->second.size) {
            ++toosmall;
            continue;
        }
        //Insert feature in references
        //std::cout<<"TODO: insert feature with name in references"<<std::endl;
        lss>>dummy;//score
        lss>>dummy;//strand
        lss>>dummy;//phase

        //parse the attributes

        std::string attr_id,attr_name;
        std::string attkv,attk,attv;
        while ( std::getline( lss, attkv, ';' ) ){
            std::stringstream iss(attkv);
            std::getline( iss, attk, '=' );
            if (attk == "ID") std::getline( iss, attr_id);
            if (attk == "Name") std::getline( iss, attr_name);
        }
        if (!attr_id.empty()) attr_name=attr_id+"|"+attr_name;
        if (attr_name.empty()) attr_name=featurename+":"+std::to_string(start)+"-"+std::to_string(end);

        reffeatures.emplace_back(ReferenceFeature(attr_name,s->second.offset+start+1,s->second.offset+end+1-S+1));


        //make the flag true from beginning to end of that feature using offset from the feature sequences
        for (uint64_t i=s->second.offset+start+1;i<s->second.offset+end+1-S+1;++i) in_feature[i]=true;
        ++used;
    }
    uint64_t fc=0;
    for (auto b:in_feature) if (b) ++fc;
    std::cout<< all_lines<<" lines, "<<noseq+toosmall+used<<" '"<<featurename
             <<"' features"<<std::endl;
    std::cout<<used<<" used, "<<tooshort<<" too short, "<<noseq<<" on nonexistent sequences, "<<toosmall<<" on invalid coordinates"<<std::endl;
    std::cout<<fc<<" / "<<in_feature.size()-(refseqs.size()*(S-1))<<" positions covered by feature"<<std::endl;
};

void SkipMerPosition::filter_repeated() {
    auto wi=skipmer_positions.begin();
    auto ri=skipmer_positions.begin();
    if ((ri+1)->skipmer!=ri->skipmer) ++wi;
    ++ri;
    for (;ri!=skipmer_positions.end()-1;++ri){
        if (ri->skipmer!=(ri-1)->skipmer and ri->skipmer!=(ri+1)->skipmer) *wi++=*ri;
    }
    if ((ri-1)->skipmer!=ri->skipmer) *wi++=*ri;
    skipmer_positions.resize(wi-skipmer_positions.begin());
    skipmer_positions.shrink_to_fit();
}

void SkipMerPosition::dump(std::string filename){
    zstr::ofstream output_file(filename+".skm_posindex.gz");
    output_file.write((char *)&m, sizeof(m));
    output_file.write((char *)&n, sizeof(n));
    output_file.write((char *)&k, sizeof(n));
    uint64_t c=refseqs.size();
    output_file.write((char *)&c, sizeof(c));
    uint64_t totalsize=0;
    for (auto &rs:refseqs) {
        c = rs.name.size();
        output_file.write((char *) &c, sizeof(c));
        output_file.write((char *) rs.name.data(), c);
        output_file.write((char *) &rs.offset, sizeof(rs.offset));
        output_file.write((char *) &rs.size, sizeof(rs.size));
        totalsize+=rs.size;
    }
    std::cout<<"Dumped header for "<<refseqs.size()<<" sequences totaling "<<totalsize<<"bp"<<std::endl;
    c=skipmer_positions.size();
    output_file.write((char *)&c, sizeof(c));
    output_file.write((char *)skipmer_positions.data(), c*sizeof(SkipMerPositionEntry));
    //output_file.close();
}; //dumps m,n,k, count and then just the binary skipmer

void SkipMerPosition::load(std::string filename){};

void SkipMerPosition::generate_stats(std::string filename){
    std::cout<<"Generating stats not yet FULLY supported"<<std::endl;
    /*size_t u=0,m=0, tm=0;
    uint64_t hist[10001]={0};
    for (auto skm:skipmers) {
        ++hist[std::min((int)skm.count,(int)10000)];
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
    std::ofstream output_file(filename+".hist");
    for (auto i=0;i<10001;++i){
        if (hist[i]) output_file<<i<<", "<<(int)hist[i]<<std::endl;
    }
    output_file.close();
    */
}; //writes both a .txt stats file and a histogram

void SkipMerPosition::merge(){};

void SkipMerPosition::intersect(){};//kat comp's equivalent