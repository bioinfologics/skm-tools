#include <iostream>
#include <fstream>
#include "deps/cxxopts/include/cxxopts.hpp"
#include "SkipMer.h"
#include "SkipMerSpectrum.h"


int main(int argc, char * argv[]) {
    std::string seq_filename,output_prefix;
    uint16_t min_freq=1,max_freq=UINT16_MAX;
    uint8_t m=1;
    uint8_t n=3;
    uint8_t k=21;
    bool fasta_input=false,stats_only=false;
    uint64_t alloc_block=10000000;

    try
    {
        cxxopts::Options options("skm-count", "Skip-mer counter");

        options.add_options()
                ("help", "Print help")
                ("i,input", "input fasta/fastq to count skip-mers from", cxxopts::value<std::string>(seq_filename))
                ("a,fasta_input", "input is a fasta file (default=0)", cxxopts::value<bool>(fasta_input))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("s,stats_only", "do not dump spectra, only stats and hist (default=0)", cxxopts::value<bool>(stats_only));
        options.add_options("Counting")
                ("f,min_freq","minimum skip-mer frequency", cxxopts::value<uint16_t>(min_freq))
                ("F,max_freq","maximum skip-mer frequency", cxxopts::value<uint16_t>(max_freq));
        options.add_options("Skip-mer shape (m every n, total k)")
                ("m,used_bases", "m", cxxopts::value<uint8_t>(m))
                ("n,skipped_bases", "n", cxxopts::value<uint8_t>(n))
                ("k,total_bases", "k", cxxopts::value<uint8_t>(k));
        options.add_options("Performance")
                ("alloc_block", "allocation block for the skip-mer sets", cxxopts::value<uint64_t>(alloc_block));


        options.parse(argc, argv);

        if (0 != options.count("help"))
        {
            std::cout << options.help({"","Skip-mer counting", "Counting", "Skip-mer shape (m every n, total k)", "Performance"}) << std::endl;
            exit(0);
        }

        if (options.count("o")!=1 or options.count("i")!=1) {
            std::cout << "Error: please specify input file and output prefix"<<std::endl
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


    std::cout<< "Welcome to skm-count"<<std::endl<<std::endl;
    print_skipmer_shape(m,n,k);
    std::cout<<std::endl;

    SkipMerSpectrum skms(m,n,k,alloc_block);
    std::cout<<"Counting..."<<std::endl;
    skms.count_from_file(seq_filename, fasta_input);
    std::cout<<skms.skipmers.size()<<" skip-mers generated"<<std::endl;
    std::cout<<"Aggregating..."<<std::endl;
    if ( min_freq!=1 or max_freq!=UINT16_MAX) skms.sort_and_collapse(min_freq,max_freq);
    else skms.sort_and_collapse();
    skms.generate_stats(output_prefix);
    if (not stats_only) {
        std::cout << "Dumping..." << std::endl;
        skms.dump(output_prefix);
    }
    std::cout<<"DONE!"<<std::endl;
    return 0;
}

