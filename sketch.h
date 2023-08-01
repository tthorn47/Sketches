#include <inttypes.h>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "hashutil.h"

#pragma once
//xorshift is used to provide quick psuedo-random numbers to hash
struct xorshift64_state {
    uint64_t a;
};

static uint64_t xorshift64(struct xorshift64_state *state)
{
	uint64_t x = state->a;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	return state->a = x;
};

template <typename T = uint64_t>
using hStore = std::pair<uint64_t, T>;

template<typename T = uint64_t>
class Sketch {
    public:
    virtual void add(T x) = 0;
    virtual uint64_t estimate(T x) = 0;
    virtual std::multimap<uint64_t, T, std::greater<uint64_t>> heavyHitters(double phi) = 0;
    virtual uint64_t size() = 0;
    virtual ~Sketch(){};

    protected:
    uint64_t numCounters;
    uint64_t n;
};

template<typename T = uint64_t>
class MisraGries : Sketch<T> {
    public:
    MisraGries(uint64_t size);
    void add(T x);
    uint64_t estimate(T x);
    std::multimap<uint64_t, T, std::greater<uint64_t>> heavyHitters(double phi);
    uint64_t size();
    ~MisraGries();

    private:
    std::unordered_map<T,double>* backingCounts;
};

template<typename T = uint64_t>
class CountSketch : Sketch<T> {
public:
    CountSketch(uint64_t hashes, uint64_t counters, uint64_t seed, uint64_t n, double phi);
    void add(T x);
    uint64_t estimate(T x);
    std::multimap<uint64_t, T, std::greater<uint64_t>> heavyHitters(double phi);
    uint64_t size();
    ~CountSketch();

private:
    uint64_t seed;
    int64_t *counterArray;
    uint64_t *salts;
    uint64_t hashes;
    double phi;

    std::unordered_map<T, double> *counts;
    std::vector<hStore<T>> *hh;
    
    std::multiset<int64_t> iterateCounters(T x, bool inc);
};


template<typename T = uint64_t>
class CountMin : Sketch<T> {
public:
    CountMin(uint64_t hashes, uint64_t counters, uint64_t seed, uint64_t n, double phi);
    void add(T x);
    uint64_t estimate(T x);
    std::multimap<uint64_t, T, std::greater<uint64_t>> heavyHitters(double phi);
    uint64_t size();
    ~CountMin();

private:
    uint64_t seed;
    int64_t *counterArray;
    uint64_t *salts;
    uint64_t hashes;
    double phi;

    std::unordered_map<T, double> *counts;
    std::vector<hStore<T>> *hh;
};