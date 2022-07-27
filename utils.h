#pragma once

#include "pcg/pcg_basic.h"

inline pcg32_random_t GetRNG()
{
    pcg32_random_t rng;
#if DETERMINISTIC()
    pcg32_srandom_r(&rng, 0x1337FEED, 0);
#else
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution<uint32_t> dist;
    pcg32_srandom_r(&rng, dist(generator), 0);
#endif
    return rng;
}

inline float RandomFloat01(pcg32_random_t& rng)
{
    return float(pcg32_random_r(&rng)) / 4294967295.0f;
}

template <typename T>
inline T Lerp(T A, T B, T t)
{
    return A * (T(1) - t) + B * t;
}

template <typename T>
T Clamp(T x, T themin, T themax)
{
    if (x <= themin)
        return themin;
    if (x >= themax)
        return themax;
    return x;
}
