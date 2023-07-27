#include "sketch.h"
#include <stdio.h>
#include <algorithm>
#include <queue>
#include <limits>
#include <cassert>

//xorshift is used to provide quick psuedo-random numbers to hash
struct xorshift64_state {
    uint64_t a;
};

uint64_t xorshift64(struct xorshift64_state *state)
{
	uint64_t x = state->a;
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	return state->a = x;
}


MisraGries::MisraGries(uint64_t size)
{
    numCounters = size;
    n = 0;
    backingCounts = new std::unordered_map<uint64_t, double>();
}

MisraGries::~MisraGries()
{
}

void MisraGries::add(uint64_t x)
{
    n++;
    if (backingCounts->count(x) != 0)
    {
        // printf("found %lu\n", x);
        (*(backingCounts))[x] += 1;
    }
    else if (backingCounts->size() < numCounters)
    {
        // printf("insert %lu\n", x);
        (*(backingCounts))[x] = 1;
    }
    else
    {
        // printf("missed %lu\n", x);
        std::vector<uint64_t> toDel;
        auto keyIter = backingCounts->begin();
        while (keyIter != backingCounts->end())
        {
            (*(backingCounts))[(*keyIter).first]--;
            if ((*(backingCounts))[(*keyIter).first] == 0)
            {
                toDel.push_back((*keyIter).first);
            }
            keyIter++;
        }
        for (size_t i = 0; i < toDel.size(); i++)
        {
            backingCounts->erase(toDel[i]);
        }
    }
}

uint64_t MisraGries::estimate(uint64_t x)
{
    if (backingCounts->count(x) != 0)
    {
        //printf("est found %lu: %lu / %lu\n", x, (*(backingCounts))[x], n);
        return (*backingCounts)[x];
    }
    return 0;
}

std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> MisraGries::heavyHitters(double phi)
{
    auto keyIter = backingCounts->begin();
    std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> ret;
    while (keyIter != backingCounts->end())
    {
        uint64_t est = (*keyIter).second;
        if ((double)est/(double)n >= phi)
        {
            ret.insert(hhStore(est,(*keyIter).first));
        }
        keyIter++; 
    }
    return ret;
}

uint64_t MisraGries::size()
{
    return sizeof(*backingCounts) + backingCounts->size() * sizeof(uint64_t) * 2 
    + sizeof(*this);
}

//-----------------------------------------------------
// COUNT SKETCH
//
//
//
//
//

CountSketch::CountSketch(uint64_t hashes, uint64_t counters, uint64_t seed, uint64_t n, double phi)
{

    hashes = hashes;
    numCounters = counters;
    this->seed = seed;
    this->phi = phi;
    this->n = n;

    xorshift64_state s;
    s.a = seed;
    salts = (uint64_t *)malloc(hashes * sizeof(uint64_t) * 2);

    for (size_t i = 0; i < hashes * 2; i++)
    {
        salts[i] = xorshift64(&s);
    }

    this->counterArray = (int64_t *)calloc(hashes * counters, sizeof(int64_t));
    hh = new std::vector<heapStore>();
    counts = new std::unordered_map<uint64_t, double>();
}

CountSketch::~CountSketch()
{
}

void CountSketch::add(uint64_t x)
{
    std::multiset<int64_t> ret;
    for (size_t i = 0; i < hashes; i++)
    {
        uint64_t toHash = x;
        uint64_t row = i * numCounters;
        uint64_t column = MurmurHash64A(&toHash, sizeof(uint64_t), seed ^ salts[i]) % numCounters;
        int64_t inc = MurmurHash64A(&toHash, sizeof(uint64_t), (seed) ^ salts[i + hashes]) % 2 == 0 ? 1 : -1;
        counterArray[row + column] += inc;
        int64_t ins =  inc * counterArray[row + column];
        ret.insert(ins);
    }

    auto medianCheck = ret.begin();
    int64_t freq;
    std::advance(medianCheck, ret.size()/2);
    if(hashes % 2 == 0){
        freq = *medianCheck;
        medianCheck--;
        freq += *medianCheck;
        freq /= 2;
    } else {
        freq = *medianCheck;
    }
    if(freq < 0)
        freq = 0;

    if ((double)freq/(double)n > phi)
    {

        if (hh->size() > 0)
        {
            auto top = hh->at(0);
            if (counts->count(x) != 0)
            { // replace item already in heap
                uint64_t oldFreq = (*(counts))[x];
                auto place = std::find(hh->begin(), hh->end(), heapStore(oldFreq, x));
                hh->erase(place);
                (*(counts))[x] = freq;
                hh->push_back((heapStore(freq, x)));
                std::make_heap(hh->begin(), hh->end(), std::greater<heapStore>());
                return;
            }
            else if (freq >= top.first && hh->size() >= hashes)
            { // kick lowest item
                std::pop_heap(hh->begin(), hh->end()); hh->pop_back();
                counts->extract(top.second);
            }
            else if (hh->size() >= hashes)
            { // not added
                return;
            }
        }

        (*(counts))[x] = freq;
        hh->push_back((heapStore(freq, x))); std::push_heap(hh->begin(), hh->end());
    }
}

uint64_t CountSketch::estimate(uint64_t x)
{
    std::multiset<int64_t> ret;
    for (size_t i = 0; i < hashes; i++)
    {
        uint64_t toHash = x;
        uint64_t row = i * numCounters;
        uint64_t column = MurmurHash64A(&toHash, sizeof(uint64_t), seed ^ salts[i]) % numCounters;
        int64_t inc = MurmurHash64A(&toHash, sizeof(uint64_t), (seed) ^ salts[i + hashes]) % 2 == 0 ? 1 : -1;
        counterArray[row + column] += inc;
        int64_t ins =  inc * counterArray[row + column];
        ret.insert(ins);
    }

    auto medianCheck = ret.begin();
    int64_t freq;
    std::advance(medianCheck, ret.size()/2);

    if(hashes % 2 == 0){
        freq = *medianCheck;
        medianCheck--;
        freq += *medianCheck;
        freq /= 2;
    } else {
        //medianCheck++;
        freq = *medianCheck;
    }

    if(freq < 0)
        freq = 0;

    return freq;
}

std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> CountSketch::heavyHitters(double phi)
{
    std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> ret;
    std::sort_heap(hh->begin(), hh->end(),std::greater<heapStore>());
    for (auto kv : *(hh))
    {
        if((double)kv.first/(double)n < phi)
            return ret;
        ret.insert(kv);
    }
    return ret;
}

uint64_t CountSketch::size()
{
    return sizeof(*counts) + (counts->size() * (sizeof(uint64_t) + sizeof(double))) 
                + sizeof(*hh) + (hh->size() * sizeof(heapStore)) 
                + sizeof(*salts) * hashes * 2
                + sizeof(*this);
}





CountMin::CountMin(uint64_t hashes, uint64_t counters, uint64_t seed, uint64_t n, double phi)
{
    hashes = hashes;
    numCounters = counters;
    seed = seed;
    phi = phi;
    this->n = n;
    xorshift64_state s;
    s.a = seed;
    salts = (uint64_t *)malloc(hashes * sizeof(uint64_t));

    for (size_t i = 0; i < hashes; i++)
    {
        salts[i] = xorshift64(&s);
    }

    this->counterArray = (int64_t *)calloc(hashes * counters, sizeof(int64_t));
    hh = new std::vector<heapStore>();
    counts = new std::unordered_map<uint64_t, double>();
}

CountMin::~CountMin()
{
}

void CountMin::add(uint64_t x)
{
    double freq = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < hashes; i++)
    {
        uint64_t toHash = x;
        uint64_t row = i * numCounters;
        uint64_t column = MurmurHash64A(&toHash, sizeof(uint64_t), seed ^ salts[i]) % numCounters;
        counterArray[row + column] += 1;
        uint64_t newFreq = (counterArray[row + column]);
        if (freq > newFreq)
        {
            freq = newFreq;
        }
    }

    if (freq/(double)n > phi)
    {

        if (hh->size() > 0)
        {
            auto top = hh->at(0);
            if (counts->count(x) != 0)
            { // replace item already in heap
                uint64_t oldFreq = (*(counts))[x];
                auto place = std::find(hh->begin(), hh->end(), heapStore(oldFreq, x));
                hh->erase(place);
                std::make_heap(hh->begin(), hh->end(), std::greater<heapStore>());
            }
            else if (freq >= top.first && hh->size() >= hashes)
            { // kick lowest item
                std::pop_heap(hh->begin(), hh->end()); hh->pop_back();
                counts->extract(top.second);
            }
            else if (hh->size() >= hashes)
            { // not added
                return;
            }
        }

        (*(counts))[x] = freq;
        hh->push_back((heapStore(freq, x))); std::push_heap(hh->begin(), hh->end());
    }
}

uint64_t CountMin::estimate(uint64_t x)
{
    double freq = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < hashes; i++)
    {
        uint64_t toHash = x + salts[i];
        uint64_t row = i * numCounters;
        uint64_t column = MurmurHash64A(&toHash, sizeof(uint64_t), seed) % numCounters;
        double newFreq = (counterArray[row + column]);
        if (freq > newFreq)
        {
            freq = newFreq;
        }
    }
    return freq;
}

std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> CountMin::heavyHitters(double phi)
{
    std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> ret;
    std::sort_heap(hh->begin(), hh->end(),std::greater<heapStore>());
    for (auto kv : *hh)
    {
        
        if((double)kv.first/(double)n < phi)
            return ret;
        ret.insert(kv);
    }
    return ret;
}

uint64_t CountMin::size()
{
    return sizeof(*counts) + (counts->size() * (sizeof(uint64_t) + sizeof(double))) 
                + sizeof(*hh) + (hh->size() * sizeof(heapStore)) + sizeof(*this);
}

// int main()
// {
//     CountMin *test = new CountMin(8, 8, 69, 0.1);
//     test->add(69);
//     test->add(35);
//     test->add(69);
//     test->add(69);
//     test->add(69);
//     test->add(69);
//     test->add(69);
//     test->add(69);
//     test->add(69);
//     test->add(69);

//     printf("%lf\n", test->estimate(69));
//     delete test;
// }