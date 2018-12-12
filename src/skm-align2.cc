#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <chrono>
#include <numeric>
#include <iomanip>
#include "deps/cxxopts/include/cxxopts.hpp"
#include "SkipMer.h"
#include "SkipMerPosition.h"
#include "SkipMerSpectrum.h"
#include "SkipMerMultiWayCoverageAnalyser.h"


struct Block {
    Block() {}
    //Block(sgNodeID_t nodeID, uint32_t readID, int32_t nstart = 0, int32_t nend = 0, int32_t qstart = 0, int32_t qend = 0, int32_t score = 0) :
    //        node(nodeID), read_id(readID), nStart(nstart), nEnd(nend), qStart(qstart), qEnd(qend), score(score) {}

    bool operator==(const Block &other) const {
        return std::tie(nStart,nEnd,qStart,qEnd)
               == std::tie(other.nStart,other.nEnd,other.qStart,other.qEnd);
    }


    bool operator<(const Block &other) const {
        return std::tie(qStart,qEnd,nStart,nEnd)
               < std::tie(other.qStart,other.qEnd,other.nStart,other.nEnd);
    }

    //sgNodeID_t node = 0;        /// Node ID, sign represents direction
    //uint32_t read_id = 0;       /// ID of the read from the Datastore   (this is never negative!)
    int64_t nStart = 0;         /// Position of the starting node kmer of this mapping
    int64_t nEnd = 0;           /// Position of the ending node kmer of this mapping
    int64_t qStart = 0;         /// Query start position
    int64_t qEnd = 0;           /// Query end position
    int32_t score = 0;          /// Alignment score
};

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

bool file_exists(const std::string &fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

void print_workdir(char *const *argv) {
    char wdir[2048];
    getcwd(wdir, sizeof(wdir));

    std::cout << "Working directory: " << wdir << std::endl;
    std::cout << "Executable: " << argv[0] << std::endl;
}


int main(int argc, char * argv[]) {

    print_workdir(argv);

    std::string target_filename,query_filename;
    std::string output_prefix;

    uint8_t m=1;
    uint8_t n=3;
    uint8_t k=21;

    int max_jump=2000;
    int max_delta_change=100;
    int min_chain=20;
    int min_size=2000;

    uint64_t alloc_block=10000000;
    int max_freq=1;

    try
    {
        cxxopts::Options options("skm-align2", "Even simpler aligner based on Skipmer direct matches");

        options.add_options()
                ("help", "Print help")
                ("t,target", "target sequence set", cxxopts::value<std::string>(target_filename))
                ("q,query", "query sequence set", cxxopts::value<std::string>(query_filename))
                ("f,max_freq", "max frequency of skip-mers to use as anchors", cxxopts::value<int>(max_freq))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Skip-mer shape (m every n, total k)")
                ("m,used_bases", "m", cxxopts::value<uint8_t>(m))
                ("n,skipped_bases", "n", cxxopts::value<uint8_t>(n))
                ("k,total_bases", "k", cxxopts::value<uint8_t>(k));

        options.add_options("Performance")
                ("alloc_block", "allocation block for the skip-mer sets", cxxopts::value<uint64_t>(alloc_block));

        options.parse(argc, argv);

        if (0 != options.count("help"))
        {
            std::cout << options.help({"","Skip-mer shape (m every n, total k)","Performance"}) << std::endl;
            exit(0);
        }

        if (options.count("o")!=1 or options.count("t")!=1 or options.count("q")!=1) {
            std::cout << "Error: please specify input files and output prefix"<<std::endl
                <<" Use option --help to check command line arguments." << std::endl;
            exit(1);
        }


    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }

#ifdef USE_LARGE_KMERS
    //option validation
    if (n<1 or n<m or k<m or k>63 or k%m!=0){
        std::cout << "Error: invalid skip-mer shape! m="<<m<<" n="<<n<<" k="<<k<<std::endl
                  << "Conditions: 0 < m <= n, k <= 63 , k must multiple of m." << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
#else
    //option validation
    if (n<1 or n<m or k<m or k>31 or k%m!=0 /*or k%2==0*/){
        std::cout << "Error: invalid skip-mer shape! m="<<m<<" n="<<n<<" k="<<k<<std::endl
                  << "Conditions: 0 < m <= n, k <= 31 , k must multiple of m." << std::endl
                  //<< "Condition for alignment: k must be odd." << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
#endif
    std::cout<< "Welcome to skm-align"<<std::endl<<std::endl;
    print_skipmer_shape(m,n,k);
    std::cout<<std::endl;


    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::time_point<std::chrono::system_clock> ini,fin;

    std::chrono::duration<double> elapsed_seconds;
    start = std::chrono::system_clock::now();

    std::cout<<"Target file: "<<target_filename<<std::endl;
    std::cout<<"Counting... "<<std::flush;
    SkipMerPosition target(m,n,k,alloc_block);
    ini = std::chrono::system_clock::now();
    target.create_from_reference_file(target_filename,true);
    uint64_t target_last=target.skipmer_positions.back().position;
    std::cout<< target.skipmer_positions.size() << " skip-mer positions, last position "<<target_last<<"."<<std::endl;
    std::cout<<"Sorting... "<<std::flush;
    target.sort();
    fin = std::chrono::system_clock::now();
    elapsed_seconds = fin-ini;
    std::cout << "Target counting time: " << elapsed_seconds.count() << "s\nDONE!\n\n";
    std::ofstream coords_file(output_prefix+".coords");
    std::cout << "Now processing queries one-by-one" <<std::endl;
    std::cout<<"Query file: "<<query_filename<<std::endl;
    std::ifstream queries_file(query_filename);
    while (!queries_file.eof()) {
        std::string sequence_string,sequence_name;
        uint64_t sequence_size;
        while(!queries_file.eof()){
            queries_file>>sequence_string;
            if (queries_file.eof()) break;
            if (sequence_string.size()>0) {
                if (sequence_string[0]=='>') {
                    sequence_name=sequence_string.substr(1);
                }
                else {
                    sequence_size=sequence_string.size();
                    break;
                }
            }
        }
        if (queries_file.eof()) break;
        std::cout << "Counting... " << std::flush;
        SkipMerPosition query(m, n, k, alloc_block);
        ini = std::chrono::system_clock::now();
        query.add_from_string(sequence_string,0);
        uint64_t query_last = query.skipmer_positions.back().position;
        std::cout << query.skipmer_positions.size() << " skip-mer positions, last position " << query_last << "."
                  << std::endl;


        std::cout << "Sorting... " << std::flush;
        query.sort();
        fin = std::chrono::system_clock::now();
        elapsed_seconds = fin - ini;
        std::cout << "Target counting time: " << elapsed_seconds.count() << "s\nDONE!\n\n";


        std::cout << "Looking for matches with <" << max_freq << " copies on both sequence sets... " << std::endl;

        std::vector<std::pair<uint64_t, uint64_t>> fwmatches;
        std::vector<std::pair<uint64_t, uint64_t>> revmatches;
        for (uint64_t ti = 0, qi = 0;
             ti < target.skipmer_positions.size() and qi < query.skipmer_positions.size(); ++ti) {
            //std::cout<<"ti: "<<ti<<"  qi: "<<qi<<std::endl;
            while (query.skipmer_positions[qi].skipmer < target.skipmer_positions[ti].skipmer and
                   qi < query.skipmer_positions.size())
                ++qi;
            if (qi == query.skipmer_positions.size()) break;
            if (query.skipmer_positions[qi].skipmer == target.skipmer_positions[ti].skipmer) {
                //TODO check f<max_freq on both sets
                //std::cout<<"hit"<<std::endl;
                uint64_t ft = 0;
                for (auto te = ti; te < target.skipmer_positions.size() and
                                   target.skipmer_positions[te].skipmer == target.skipmer_positions[ti].skipmer; ++te)
                    ++ft;
                //std::cout<<"ft: "<<ft<<std::endl;
                uint64_t fq = 0;
                for (auto qe = qi; qe < query.skipmer_positions.size() and
                                   query.skipmer_positions[qe].skipmer == query.skipmer_positions[qi].skipmer; ++qe)
                    ++fq;
                //std::cout<<"fq: "<<fq<<std::endl;
                if (ft <= max_freq and fq <= max_freq) {
                    for (auto tii = ti; tii < ti + ft; ++tii) {
                        for (auto qii = qi; qii < query.skipmer_positions.size() and
                                            query.skipmer_positions[qii].skipmer ==
                                            target.skipmer_positions[tii].skipmer; ++qii) {
                            if (query.skipmer_positions[qii].rev == target.skipmer_positions[tii].rev)
                                fwmatches.emplace_back(query.skipmer_positions[qii].position,
                                                       target.skipmer_positions[tii].position
                                );
                            else
                                revmatches.emplace_back(query_last - query.skipmer_positions[qii].position,
                                                        target.skipmer_positions[tii].position);
                        }
                    }
                }
                ti = ti + ft - 1;
                qi = qi + fq - 1;
            }
        }
        std::sort(fwmatches.begin(), fwmatches.end());
        std::sort(revmatches.begin(), revmatches.end());
        {
            std::cout << "Chaining " << fwmatches.size() << " FW matches... " << std::endl;
            int64_t start_p = 0;
            int64_t start_t = 0;
            int64_t last_delta = 0;
            int64_t last_t = -1000;
            int64_t last_p = -1000;
            int64_t chain = -1;
            auto &chits = fwmatches;
            std::vector<bool> used(chits.size());
            //std::sort(chits.begin(),chits.end());
            //Filter hits on overused positions
            for (auto starti = 0; starti < chits.size(); ++starti) {
                auto endi = starti;
                while (endi < chits.size() - 1 and chits[starti].first == chits[endi + 1].first) ++endi;
                bool discard = false;
                if (endi > starti + 4) discard = true; //more than 5 hits to candidate ->discard
                else {
                    for (auto x = starti; x < endi; ++x) {
                        if (chits[x + 1].second - chits[x].second < max_jump) discard = true;
                    }
                }
                if (discard) for (auto x = starti; x <= endi; ++x) used[x] = true;
            }

            for (auto starti = 0; starti < chits.size(); ++starti) {
                if (used[starti]) continue;
                chain = 1;
                start_p = chits[starti].first;
                start_t = chits[starti].second;
                last_p = chits[starti].first;
                last_t = chits[starti].second;
                last_delta = last_t - last_p;
                used[starti] = true;
                //std::cout<<"New chain started with "<<start_p<<"->"<<start_t<<std::endl;
                for (auto i = starti + 1; i < chits.size() and chits[i].first - last_p <= max_jump; ++i) {
                    auto &h = chits[i];
                    //if not in chain, continue;
                    if (used[i] or last_p > h.first or last_t > h.second or h.second - last_t > max_jump
                        or llabs(last_delta - (h.second - h.first)) > max_delta_change) {
                        //std::cout<<" skipping "<<h.first<<"->"<<h.second<<std::endl;
                        continue;
                    }
                    ++chain;
                    last_p = h.first;
                    last_t = h.second;
                    last_delta = last_t - last_p;
                    used[i] = true;
                    //std::cout<<" chain is now length "<<chain<<" after adding "<<h.first<<"->"<<h.second<<" (delta "<<h.first-h.second<<")"<<std::endl;
                }
                if (chain >= min_chain and last_p - start_p >= min_size) {
                    Block b;
                    b.qStart = start_p;
                    b.qEnd = last_p;
                    b.nStart = start_t;
                    b.nEnd = last_t;
                    b.score = chain;
                    int i=0;
                    for (; start_t > target.refseqs[i].offset+target.refseqs[i].size; ++i );
                    coords_file<<start_t-target.refseqs[i].offset<<"\t"<<last_t-target.refseqs[i].offset<<"\t"<<start_p<<"\t"<<last_p<<"\t"<<last_t-start_t<<"\t"<<last_p-start_p
                            <<"\t"<<std::fixed<<std::setprecision(2)<<100.0*chain/(last_t-start_t)
                            <<"\t"<< target.refseqs[i].size << "\t" << query_last << "\t" << target.refseqs[i].name << "\t" << sequence_name <<std::endl;
                }

            }
        }
        //std::cout << "There are now " << blocks.size() << " blocks" << std::endl;
        {
            std::cout << "Chaining " << fwmatches.size() << " REV matches... " << std::endl;
            int64_t start_p = 0;
            int64_t start_t = 0;
            int64_t last_delta = 0;
            int64_t last_t = -1000;
            int64_t last_p = -1000;
            int64_t chain = -1;
            auto &chits = revmatches;
            std::vector<bool> used(chits.size());
            //std::sort(chits.begin(),chits.end());
            //Filter hits on overused positions
            for (auto starti = 0; starti < chits.size(); ++starti) {
                auto endi = starti;
                while (endi < chits.size() - 1 and chits[starti].first == chits[endi + 1].first) ++endi;
                bool discard = false;
                if (endi > starti + 4) discard = true; //more than 5 hits to candidate ->discard
                else {
                    for (auto x = starti; x < endi; ++x) {
                        if (chits[x + 1].second - chits[x].second < max_jump) discard = true;
                    }
                }
                if (discard) for (auto x = starti; x <= endi; ++x) used[x] = true;
            }

            for (auto starti = 0; starti < chits.size(); ++starti) {
                if (used[starti]) continue;
                chain = 1;
                start_p = chits[starti].first;
                start_t = chits[starti].second;
                last_p = chits[starti].first;
                last_t = chits[starti].second;
                last_delta = last_t - last_p;
                used[starti] = true;
                //std::cout<<"New chain started with "<<start_p<<"->"<<start_t<<std::endl;
                for (auto i = starti + 1; i < chits.size() and chits[i].first - last_p <= max_jump; ++i) {
                    auto &h = chits[i];
                    //if not in chain, continue;
                    if (used[i] or last_p > h.first or last_t > h.second or h.second - last_t > max_jump
                        or llabs(last_delta - (h.second - h.first)) > max_delta_change) {
                        //std::cout<<" skipping "<<h.first<<"->"<<h.second<<std::endl;
                        continue;
                    }
                    ++chain;
                    last_p = h.first;
                    last_t = h.second;
                    last_delta = last_t - last_p;
                    used[i] = true;
                    //std::cout<<" chain is now length "<<chain<<" after adding "<<h.first<<"->"<<h.second<<" (delta "<<h.first-h.second<<")"<<std::endl;
                }
                if (chain >= min_chain and last_p - start_p >= min_size) {
                    Block b;
                    b.qStart = query_last - start_p;
                    b.qEnd = query_last - last_p;
                    b.nStart = start_t;
                    b.nEnd = last_t;
                    b.score = chain;
                    int i=0;
                    for (; start_t > target.refseqs[i].offset+target.refseqs[i].size; ++i );
                    coords_file<<start_t-target.refseqs[i].offset<<"\t"<<last_t-target.refseqs[i].offset<<"\t"<<query_last - start_p<<"\t"<<query_last - last_p<<"\t"<<last_t-start_t<<"\t"<<last_p-start_p
                             <<"\t"<<std::fixed<<std::setprecision(2)<<100.0*chain/(last_t-start_t)
                             <<"\t"<< target.refseqs[i].size << "\t" << query_last << "\t" << target.refseqs[i].name << "\t" << sequence_name <<std::endl;
                    //blocks.push_back(b);
                }

            }
        }

        //std::cout << "There are now " << blocks.size() << " blocks" << std::endl;
    }
    /*std::ofstream mf(output_prefix+".matches");
    std::ofstream fgpdf(output_prefix+"_f.gpd");
    std::ofstream rgpdf(output_prefix+"_r.gpd");
    for (auto b:blocks){
        mf<<b.nStart<<":"<<b.nEnd<<" -> "<<b.qStart<<":"<<b.qEnd<<" ( "<<b.score<<" hits, "<<(b.score*100/(b.nEnd-b.nStart+1))<<"% )"<<std::endl;
        if (b.qStart<b.qEnd) fgpdf<<b.nStart<<" "<<b.qStart<<std::endl<<b.nEnd<<" "<<b.qEnd<<std::endl<<std::endl;
        else rgpdf<<b.nStart<<" "<<b.qStart<<std::endl<<b.nEnd<<" "<<b.qEnd<<std::endl<<std::endl;
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << "Total time: " << elapsed_seconds.count() << "s\n";

    std::cout<<"DONE!"<<std::endl;*/
    return 0;
}
