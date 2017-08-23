//
// Created by Bernardo Clavijo (TGAC) on 27/06/2017.
//

#include <iostream>
#include "SkipMer.h"

void print_skipmer_shape(uint8_t _n, uint8_t _m, uint8_t _k){
    int K=_k;
    K=K+((K-1)/(int)_n)*((int)_m-(int)_n);
    std::cout<<"SkipMer("<<(int)_n<<", "<<(int)_m<<", "<<(int)_k<<") spans over "<<K<<"bp"<<std::endl;
    std::cout<<"Shape: ";
    for (auto i=0;i<K;++i) {
        if (i % _m < _n) std::cout << 'X';
        else std::cout << '-';
    }
    std::cout<<std::endl;
}