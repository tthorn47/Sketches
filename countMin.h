#include "sketch.h"




template<typename T>
CountMin<T>::CountMin(uint64_t hashes, uint64_t counters, uint64_t seed, uint64_t n, double phi)
{
    this->hashes = hashes;
    this->numCounters = counters;
    this->seed = seed;
    this->phi = phi;
    this->n = n;
    
    xorshift64_state s;
    s.a = seed;
    salts = (uint64_t *)malloc(hashes * sizeof(uint64_t));

    for (size_t i = 0; i < hashes; i++)
    {
        salts[i] = xorshift64(&s);
    }

    this->counterArray = (int64_t *)calloc(hashes * counters, sizeof(int64_t));
    hh = new std::vector<hStore<T>>();
    counts = new std::unordered_map<uint64_t, double>();
}

template<typename T>
CountMin<T>::~CountMin()
{
}

template<typename T>
void CountMin<T>::add(T x)
{
    double freq = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < hashes; i++)
    {
        T toHash = x;
        uint64_t row = i * this->numCounters;
        uint64_t column = MurmurHash64A(&toHash, sizeof(T), seed ^ salts[i]) % this->numCounters;
        counterArray[row + column] += 1;
        uint64_t newFreq = (counterArray[row + column]);
        if (freq > newFreq)
        {
            freq = newFreq;
        }
    }

    if (freq/(double)this->n > phi)
    {

        if (hh->size() > 0)
        {
            auto top = hh->at(0);
            if (counts->count(x) != 0)
            { // replace item already in heap
                uint64_t oldFreq = (*(counts))[x];
                auto place = std::find(hh->begin(), hh->end(), hStore(oldFreq, x));
                hh->erase(place);
                std::make_heap(hh->begin(), hh->end(), std::greater<hStore<T>>());
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
        hh->push_back((hStore(freq, x))); std::push_heap(hh->begin(), hh->end());
    }
}

template<typename T>
uint64_t CountMin<T>::estimate(T x)
{
    double freq = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < hashes; i++)
    {
        uint64_t toHash = x + salts[i];
        uint64_t row = i * this->numCounters;
        uint64_t column = MurmurHash64A(&toHash, sizeof(uint64_t), seed) % this->numCounters;
        double newFreq = (counterArray[row + column]);
        if (freq > newFreq)
        {
            freq = newFreq;
        }
    }
    return freq;
}

template<typename T>
std::multimap<uint64_t, T, std::greater<uint64_t>> CountMin<T>::heavyHitters(double phi)
{
    std::multimap<uint64_t, T, std::greater<uint64_t>> ret;
    std::sort_heap(hh->begin(), hh->end(),std::greater<hStore<T>>());
    for (auto kv : *hh)
    {
        
        if((double)kv.first/(double)this->n < phi)
            return ret;
        ret.insert(kv);
    }
    return ret;
}

template<typename T>
uint64_t CountMin<T>::size()
{
    return sizeof(*counts) + (counts->size() * (sizeof(uint64_t) + sizeof(double))) 
                + sizeof(*hh) + (hh->size() * sizeof(hStore<T>)) + sizeof(*this);
}