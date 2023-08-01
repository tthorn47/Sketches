#include "sketch.h"

//-----------------------------------------------------
// COUNT SKETCH
//
//
//
//
//
template <typename T>
CountSketch<T>::CountSketch(uint64_t hashes, uint64_t counters, uint64_t seed, uint64_t n, double phi)
{
    assert(hashes > 0);
    assert(counters > 0);
    assert(n > 0);
    assert(phi > 0);
    
    this->hashes = hashes;
    this->numCounters = counters;
    this->seed = seed;
    this->phi = phi;
    this->n = n;

    xorshift64_state s;
    s.a = seed;

    salts = (uint64_t *)malloc(hashes * sizeof(uint64_t) * 2);
    this->counterArray = (int64_t *)calloc(hashes * counters, sizeof(int64_t));

    if(salts == NULL || this->counterArray == NULL){
        std::cout << "salts/counter malloc failed!\n";
        exit(1); 
    }

    for (size_t i = 0; i < hashes * 2; i++)
    {
        salts[i] = xorshift64(&s);
    }

    hh = new std::vector<hStore<T>>();
    counts = new std::unordered_map<uint64_t, double>();
}

template <typename T>
CountSketch<T>::~CountSketch()
{
}

template <typename T>
std::multiset<int64_t> CountSketch<T>::iterateCounters(T x, bool inc){
    std::multiset<int64_t> ret;
    for (size_t i = 0; i < hashes; i++)
    {
        uint64_t toHash = x;
        uint64_t row = i * this->numCounters;
        uint64_t column = MurmurHash64A(&toHash, sizeof(uint64_t), seed ^ salts[i]) % this->numCounters;
        int64_t inc = MurmurHash64A(&toHash, sizeof(uint64_t), (seed) ^ salts[i + hashes]) % 2 == 0 ? 1 : -1;
        counterArray[row + column] += inc;
        int64_t ins = inc * counterArray[row + column];
        ret.insert(ins);
    }
    return ret;
}

template <typename T>
void CountSketch<T>::add(T x)
{
    std::multiset<int64_t> ret = iterateCounters(x, true);

    auto medianCheck = ret.begin();

    int64_t freq;
    if (hashes != 1)
    {
        std::advance(medianCheck, ret.size() / 2);
        if (hashes % 2 == 0)
        {
            freq = *medianCheck;
            medianCheck--;
            freq += *medianCheck;
            freq /= 2;
        }
        else
        {
            freq = *medianCheck;
        }
    }
    else
    {
        freq = *medianCheck;
    }

    if (freq < 0)
        freq = 0;

    if ((double)freq / (double)this->n > phi)
    {

        if (hh->size() > 0)
        {
            auto top = hh->at(0);
            if (counts->count(x) != 0)
            { // replace item already in heap
                uint64_t oldFreq = (*(counts))[x];
                auto place = std::find(hh->begin(), hh->end(), hStore(oldFreq, x));
                hh->erase(place);
                (*(counts))[x] = freq;
                hh->push_back((hStore(freq, x)));
                std::make_heap(hh->begin(), hh->end(), std::greater<hStore<T>>());
                return;
            }
            else if (freq >= top.first && hh->size() >= hashes)
            { // kick lowest item
                std::pop_heap(hh->begin(), hh->end());
                hh->pop_back();
                counts->extract(top.second);
            }
            else if (hh->size() >= hashes)
            { // not added
                return;
            }
        }

        (*(counts))[x] = freq;
        hh->push_back((hStore(freq, x)));
        std::push_heap(hh->begin(), hh->end());
    }
}

template <typename T>
uint64_t CountSketch<T>::estimate(T x)
{
    std::multiset<int64_t> ret;
    for (size_t i = 0; i < hashes; i++)
    {
        uint64_t toHash = x;
        uint64_t row = i * this->numCounters;
        uint64_t column = MurmurHash64A(&toHash, sizeof(uint64_t), seed ^ salts[i]) % this->numCounters;
        int64_t inc = MurmurHash64A(&toHash, sizeof(uint64_t), (seed) ^ salts[i + hashes]) % 2 == 0 ? 1 : -1;
        counterArray[row + column] += inc;
        int64_t ins = inc * counterArray[row + column];
        ret.insert(ins);
    }

    auto medianCheck = ret.begin();
    int64_t freq;
    std::advance(medianCheck, ret.size() / 2);

    if (hashes % 2 == 0)
    {
        freq = *medianCheck;
        medianCheck--;
        freq += *medianCheck;
        freq /= 2;
    }
    else
    {
        // medianCheck++;
        freq = *medianCheck;
    }

    if (freq < 0)
        freq = 0;

    return freq;
}

template <typename T>
std::multimap<uint64_t, T, std::greater<uint64_t>> CountSketch<T>::heavyHitters(double phi)
{
    std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> ret;
    std::sort_heap(hh->begin(), hh->end(), std::greater<hStore<T>>());
    for (auto kv : *(hh))
    {
        if ((double)kv.first / (double)this->n < phi)
            return ret;
        ret.insert(kv);
    }
    return ret;
}

template <typename T>
uint64_t CountSketch<T>::size()
{
    return sizeof(*counts) + (counts->size() * (sizeof(uint64_t) + sizeof(double))) + sizeof(*hh) + (hh->size() * sizeof(hStore<T>)) + sizeof(*salts) * hashes * 2 + sizeof(*this);
}
