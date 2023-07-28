/**





*/


#include "sketch.h"


template<typename T>
MisraGries<T>::MisraGries(uint64_t size)
{
    this->numCounters = size;
    this->n = 0;
    backingCounts = new std::unordered_map<uint64_t, double>();
}

template<typename T>
MisraGries<T>::~MisraGries()
{
}

template<typename T>
void MisraGries<T>::add(T x)
{
    this->n++;
    if (backingCounts->count(x) != 0)
    {
        // printf("found %lu\n", x);
        (*(backingCounts))[x] += 1;
    }
    else if (backingCounts->size() < this->numCounters)
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

template<typename T>
uint64_t MisraGries<T>::estimate(T x)
{
    if (backingCounts->count(x) != 0)
    {
        //printf("est found %lu: %lu / %lu\n", x, (*(backingCounts))[x], n);
        return (*backingCounts)[x];
    }
    return 0;
}

template<typename T>
std::multimap<uint64_t, T, std::greater<uint64_t>> MisraGries<T>::heavyHitters(double phi)
{
    auto keyIter = backingCounts->begin();
    std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> ret;
    while (keyIter != backingCounts->end())
    {
        uint64_t est = (*keyIter).second;
        if ((double)est/(double)this->n >= phi)
        {
            ret.insert(hStore<T>(est,(*keyIter).first));
        }
        keyIter++; 
    }
    return ret;
}

template<typename T>
uint64_t MisraGries<T>::size()
{
    return sizeof(*backingCounts) + backingCounts->size() * sizeof(uint64_t) * 2 
    + sizeof(*this);
}