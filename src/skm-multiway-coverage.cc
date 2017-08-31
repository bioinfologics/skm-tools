#include <iostream>
#include <fstream>
#include "deps/cxxopts/include/cxxopts.hpp"
#include "SkipMer.h"
#include "SkipMerPosition.h"
#include "SkipMerSpectrum.h"
#include "SkipMerMultiWayCoverageAnalyser.h"


void print_workdir(char *const *argv) {
    char wdir[2048];
    getcwd(wdir, sizeof(wdir));

    std::cout << "Working directory: " << wdir << std::endl;
    std::cout << "Executable: " << argv[0] << std::endl;
}


int main(int argc, char * argv[]) {

    print_workdir(argv);

    std::vector<std::string> seq_filenames;
    std::string ref_filename;
    std::string output_prefix;
    std::string gff3_filename, gff3_feature;
    uint32_t minfreq=0;

    uint8_t m=1;
    uint8_t n=3;
    uint8_t k=21;

    uint64_t alloc_block=10000000,minfeature=1000;
    uint8_t maxcoverage=1;
    bool fasta_input=true,stats_only=false;

    try
    {
        cxxopts::Options options("skm-multiway-coverage", "Multi-way coverage projection");

        options.add_options()
                ("help", "Print help")
                ("r,reference", "reference to map coverage to", cxxopts::value<std::string>(ref_filename))
                ("gff3_file", "gff3 file for features on the first genome", cxxopts::value<std::string>(gff3_filename))
                ("gff3_feature", "gff3 feature name to split classification", cxxopts::value<std::string>(gff3_feature))
                ("gff3_min_size", "gff3 feature minimum size (default 1000)", cxxopts::value<uint64_t >(minfeature))
                ("c,coverage_set", "coverage set (use as many as needed)", cxxopts::value<std::vector<std::string>>(seq_filenames))
                ("max_coverage", "max coverage (both on ref and sets - default 1)", cxxopts::value<uint8_t >(maxcoverage))
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

        if (options.count("o")!=1 /*or options.count("i")<2*/) {
            std::cout << "Error: please specify input files and output prefix"<<std::endl
                <<" Use option --help to check command line arguments." << std::endl;
            exit(1);
        }

        if (options.count("gff3_file")!=options.count("gff3_feature")) {
            std::cout << "Error: please specify gff3_file and gff3_feature together"<<std::endl
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

    std::cout<<"Performing "<<refcount<<"-way coverage projection."<<std::endl;

    std::cout<<"Reference file: "<<ref_filename<<std::endl;
    std::cout<<"Counting... "<<std::flush;
    SkipMerPosition ref(m,n,k,alloc_block);
    ref.create_from_reference_file(ref_filename,true);
    std::cout<< ref.skipmer_positions.size() << " skip-mer positions."<<std::endl;
    std::cout<<"Sorting... "<<std::flush;
    ref.sort();
    if (!gff3_filename.empty()) {
        ref.mark_with_gff3_feature(gff3_filename, gff3_feature, minfeature);
        std::cout << "the feature collection has " << ref.reffeatures.size() << " features" << std::endl;
    }
    std::ofstream prog_cov(output_prefix+"_cov_progression.csv");
    prog_cov<<"Shape,S,Set Count,Total Ref,Total Any,Total All,Feature Ref,Feature Any,Feature All,No Feature Ref,No Feature Any,No Feature All,PTotal Ref,PTotal Any,PTotal All,PFeature Ref,PFeature Any,PFeature All,PNo Feature Ref,PNo Feature Any,PNo Feature All"<<std::endl;
    SkipMerMultiWayCoverageAnalyser mwca(ref,maxcoverage);
    for (auto ref_id=0;ref_id<refcount;++ref_id) {
        std::cout<<"Processing sequence file #"<<ref_id+1<<": "<<seq_filenames[ref_id]<<std::endl;
        SkipMerSpectrum skms(m,n,k,alloc_block);
        std::cout<<"Counting... "<<std::flush;
        skms.count_from_file(seq_filenames[ref_id],true);
        std::cout<<"Sorting and collapsing... "<<std::flush;
        skms.sort_and_collapse();
        std::cout<<"Adding coverage track... "<<std::flush;
        mwca.add_coverage_track(skms,maxcoverage);
        prog_cov<<(int) m<<"-"<<(int) n<<'-'<<(int) k<<","<<skms.S<<","<<ref_id+1;
        auto cs=mwca.covered_by_all_stats();
        for (auto covx:cs) prog_cov<<","<<covx;
        prog_cov<<","<<cs[0]-cs[3]<<","<<cs[1]-cs[4]<<","<<cs[2]-cs[5];
        cs=mwca.covered_by_all_projected_stats();
        for (auto covx:cs) prog_cov<<","<<covx;
        prog_cov<<","<<cs[0]-cs[3]<<","<<cs[1]-cs[4]<<","<<cs[2]-cs[5];
        prog_cov<<std::endl;
        std::cout<<"DONE!"<<std::endl;
    }
    std::cout<<std::endl<<"Dumping feature conservation scores... "<<std::flush;
    if (!gff3_filename.empty()) mwca.dump_features_conservation_scores(output_prefix);
    std::cout<<"DONE!"<<std::endl;
    return 0;

}
