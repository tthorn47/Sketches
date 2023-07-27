#include <inttypes.h>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "hashutil.h"

using heapStore = std::pair<uint64_t, uint64_t>;
using hhStore = std::pair<uint64_t, uint64_t>;

class Sketch {
    public:
    virtual void add(uint64_t x) = 0;
    virtual uint64_t estimate(uint64_t x) = 0;
    virtual std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> heavyHitters(double phi) = 0;
    virtual uint64_t size() = 0;
    virtual ~Sketch(){};

    protected:
    uint64_t numCounters;
    uint64_t n;
};

class MisraGries : Sketch {
    public:
    MisraGries(uint64_t size);
    void add(uint64_t x);
    uint64_t estimate(uint64_t x);
    std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> heavyHitters(double phi);
    uint64_t size();
    ~MisraGries();

    private:
    std::unordered_map<uint64_t,double>* backingCounts;
};

class CountSketch : Sketch {
public:
    CountSketch(uint64_t hashes, uint64_t counters, uint64_t seed, uint64_t n, double phi);
    void add(uint64_t x);
    uint64_t estimate(uint64_t x);
    std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> heavyHitters(double phi);
    uint64_t size();
    ~CountSketch();

private:
    uint64_t seed;
    int64_t *counterArray;
    uint64_t *salts;
    uint64_t hashes;
    double phi;

    std::unordered_map<uint64_t, double> *counts;
    std::vector<heapStore> *hh;
};

class CountMin : Sketch {
public:
    CountMin(uint64_t hashes, uint64_t counters, uint64_t seed, uint64_t n, double phi);
    void add(uint64_t x);
    uint64_t estimate(uint64_t x);
    std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> heavyHitters(double phi);
    uint64_t size();
    ~CountMin();

private:
    uint64_t seed;
    int64_t *counterArray;
    uint64_t *salts;
    uint64_t hashes;
    double phi;

    std::unordered_map<uint64_t, double> *counts;
    std::vector<heapStore> *hh;
};