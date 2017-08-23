#ifndef SKIPMER_H
#define SKIPMER_H


#include <cstdint>

#ifdef USE_LARGE_KMERS
typedef __int128 int128_t;
typedef unsigned __int128 uint128_t;

typedef uint128_t SkipMer;
#else
typedef uint64_t SkipMer;
#endif


void print_skipmer_shape(uint8_t _n, uint8_t _m, uint8_t _k);

#endif //KMATCH_SKIPMER_H
