#include <iostream>
#include <fstream>
#include "deps/cxxopts/include/cxxopts.hpp"
#include "SkipMer.h"
#include "SkipMerPosition.h"
#include "SkipMerSpectrum.h"

void print_workdir(char *const *argv) {
    char wdir[2048];
    getcwd(wdir, sizeof(wdir));

    std::cout << "Working directory: " << wdir << std::endl;
    std::cout << "Executable: " << argv[0] << std::endl;
}

void multiway_intersection(std::vector<SkipMerSpectrum> & _skmp, std::string out_prefix){
    //first point all skmp to firs element;

    //populate is_multi (i.e. for each element check if the next element is still the same one)
    unsigned total_sets=_skmp.size();
    unsigned active_sets=_skmp.size();
    bool is_active[_skmp.size()];
    for (auto &x:is_active) x=true;
    uint64_t unique_shares[_skmp.size()][_skmp.size()];
    for (auto &x: unique_shares) for (auto &y:x) y=0;
    uint64_t anchor_shares[_skmp.size()][_skmp.size()];
    for (auto &x: anchor_shares) for (auto &y:x) y=0;
    uint64_t anchor_feature_shares[_skmp.size()][_skmp.size()];
    for (auto &x: anchor_feature_shares) for (auto &y:x) y=0;
    uint64_t all_shares[_skmp.size()][_skmp.size()];
    for (auto &x: all_shares) for (auto &y:x) y=0;
    uint64_t idx[_skmp.size()]={0};
    for (auto &x:idx) x=0;
    SkipMer min;
    std::vector<uint8_t> smallers;
    uint64_t grandtotal=0;
    //uint64_t itr=0,next_point=1000000;
    //std::ofstream anchorsf(out_prefix+"_anchors.csv");
    //anchorsf<<"Pos1";
    //for (auto i=1;i<total_sets;++i) anchorsf<<",Pos"<<i+1;
    //anchorsf<<std::endl;
    while(active_sets){

        //find the smaller (s)
        smallers.clear();
        for (uint8_t i=0;i<total_sets;++i){
            if (is_active[i]){
                if (smallers.size()==0 or min>_skmp[i].skipmers[idx[i]].skipmer){
                    min=_skmp[i].skipmers[idx[i]].skipmer;
                    smallers.clear();
                    smallers.push_back(i);
                }else if (min==_skmp[i].skipmers[idx[i]].skipmer){
                    smallers.push_back(i);
                }
            }
        }

        bool real_anchor=true;
        for (auto si:smallers) {
            if (_skmp[si].skipmers[idx[si]].count>1){
                real_anchor=false;
                break;
            }

        }


        for (auto si:smallers) {
            //update count on shared skip-mers.
            all_shares[si][smallers.size()-1]+=_skmp[si].skipmers[idx[si]].count;

            if (_skmp[si].skipmers[idx[si]].count==1)
                ++unique_shares[si][smallers.size()-1];
            if (real_anchor) {
                ++anchor_shares[si][smallers.size()-1];

            }

            //update index
            ++idx[si];

            if (idx[si]==_skmp[si].skipmers.size()) {
                is_active[si]=false;
                --active_sets;
            }
        }
        grandtotal++;

    }
    std::ofstream statsf(out_prefix+"_mwstats.csv");
    statsf<<"Sample,Total,Unique,Non-unique";
    for (auto i=0;i<total_sets;++i) statsf<<",T"<<i+1;
    for (auto i=0;i<total_sets;++i) statsf<<",U"<<i+1;
    for (auto i=0;i<total_sets;++i) statsf<<",A"<<i+1;
    statsf<<std::endl;
    for (auto i=0;i<total_sets;++i){
        uint64_t st=0,su=0;
        for (auto j=0;j<total_sets;++j) {
            st+=all_shares[i][j];
            su+=unique_shares[i][j];
        }
        statsf<<i<<","<<st<<","<<su<<","<<st-su;
        for (auto j=0;j<total_sets;++j) statsf<<","<<all_shares[i][j];
        for (auto j=0;j<total_sets;++j) statsf<<","<<unique_shares[i][j];
        for (auto j=0;j<total_sets;++j) statsf<<","<<anchor_shares[i][j];
        statsf<<std::endl;
    }

    std::cout<<"Grand total of analysed skip-mers: "<<grandtotal<<std::endl;
};

int main(int argc, char * argv[]) {

    print_workdir(argv);

    std::vector<std::string> seq_filenames;

    std::string output_prefix;
    uint32_t minfreq=0;

    uint8_t m=1;
    uint8_t n=3;
    uint8_t k=21;

    uint64_t alloc_block=10000000,minfeature=1000;
    uint8_t maxcoverage=1;
    bool fasta_input=true,stats_only=false;

    try
    {
        cxxopts::Options options("skm-multiway-shared", "Multi-way shared skip-mer intersection");

        options.add_options()
                ("help", "Print help")
                ("s,sequence", "sequence set (use as many as needed)", cxxopts::value<std::vector<std::string>>(seq_filenames))
                //("max_coverage", "max coverage (both on ref and sets - default 1)", cxxopts::value<uint8_t >(maxcoverage))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Skip-mer shape (m every n, total k)")
                ("m", "m", cxxopts::value<uint8_t>(m))
                ("n", "n", cxxopts::value<uint8_t>(n))
                ("k", "k", cxxopts::value<uint8_t>(k));

        options.add_options("Performance")
                ("alloc_block", "allocation block for the skip-mer sets", cxxopts::value<uint64_t>(alloc_block));

        options.parse(argc, argv);

        if (0 != options.count("help"))
        {
            std::cout << options.help({"","Skip-mer shape (m every n, total k)","Performance"}) << std::endl;
            exit(0);
        }

        if (options.count("o")!=1 or options.count("s")<2) {
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
    if (n<1 or n<m or k<m or k>31 or k%m!=0){
        std::cout << "Error: invalid skip-mer shape! m="<<m<<" n="<<n<<" k="<<k<<std::endl
                  << "Conditions: 0 < m <= n, k <= 31 , k must multiple of m." << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
#endif
    std::cout<< "Welcome to skm-multiway-coverage"<<std::endl<<std::endl;
    print_skipmer_shape(m,n,k);
    std::cout<<std::endl;
    unsigned i=0;
    auto refcount=seq_filenames.size();
    std::vector<uint64_t[6]> final_counts(refcount);


    std::cout<<"Performing "<<refcount<<"-way coverage projection."<<std::endl;
    std::vector<SkipMerSpectrum> skm_counts;
    for (auto ref_id=0;ref_id<refcount;++ref_id) {
        std::cout<<"Processing sequence file #"<<ref_id+1<<": "<<seq_filenames[ref_id]<<std::endl;
        skm_counts.emplace_back(m,n,k,alloc_block);
        std::cout<<"Counting... "<<std::flush;
        skm_counts.back().count_from_file(seq_filenames[ref_id],true);
        std::cout<<"Sorting and collapsing... "<<std::flush;
        skm_counts.back().sort_and_collapse();
        std::cout<<"DONE!"<<std::endl;
    }
//    for (auto &fc:final_counts) fc[0]=fc[1]=fc[2]=fc[3]=fc[4]=fc[5]=0;
//
//    std::cout<<std::endl<<"Scanning through sorted sets..."<<std::flush;
//    std::ofstream of(output_prefix+"_shared.csv");
//    of<<"Set#,Total,TotalExclusive,Distinct,DistinctExclusive,Unique,UniqueExclusive"<<std::endl;
    multiway_intersection(skm_counts,output_prefix);
    return 0;

}
