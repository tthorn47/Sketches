
#include <iostream>
#include <cassert>
#include <chrono>
#include <openssl/rand.h>
#include <map>
#include <unordered_map>
#include <fstream>
#include "sketch.h"
#include "misra.h"
#include "count.h"
#include "countMin.h"
#include "zipf.h"

using namespace std::chrono;

#define UNIVERSE 1ULL << 30
#define EXP 1.5 

double elapsed(high_resolution_clock::time_point t1, high_resolution_clock::time_point t2) {
	return (duration_cast<duration<double> >(t2 - t1)).count();
}

int main(int argc, char** argv)
{
	
	if (argc < 3) {
		std::cerr << "Specify the number of items N and phi.\n";
		exit(1);
	}
	uint64_t N = atoi(argv[1]);
	double tau = atof(argv[2]);
	uint64_t *numbers = (uint64_t *)malloc(N * sizeof(uint64_t));
	if(!numbers) {
		std::cerr << "Malloc numbers failed.\n";
		exit(0);
	}
	MisraGries mg(16000);
	CountSketch cs(160,160,69,N,tau);
	CountMin cm(160, 160, 69, N, tau);

	high_resolution_clock::time_point t1, t2;
	t1 = high_resolution_clock::now();
	generate_random_keys(numbers, UNIVERSE, N, EXP);
	t2 = high_resolution_clock::now();
	std::cout << "Time to generate " << N << " items: " << elapsed(t1, t2) << " secs\n";

	std::unordered_map<uint64_t, uint64_t> map(N);

	t1 = high_resolution_clock::now();
	for (uint64_t i = 0; i < N; ++i) 
		map[numbers[i]]++;
	t2 = high_resolution_clock::now();
	std::cout << "HashTable Time to count " << N << " items: " << elapsed(t1, t2) << " secs\n";
	double htTime = elapsed(t1, t2);
	t1 = high_resolution_clock::now();
	for (uint64_t i = 0; i < N; ++i) 
		mg.add(numbers[i]);
	t2 = high_resolution_clock::now();
	std::cout << "Misra Time to count " << N << " items: " << elapsed(t1, t2) << " secs\n";
	double mgTime = elapsed(t1, t2);
	t1 = high_resolution_clock::now();
	for (uint64_t i = 0; i < N; ++i) 
		cs.add(numbers[i]);
	t2 = high_resolution_clock::now();
	std::cout << "CountSketch Time to count " << N << " items: " << elapsed(t1, t2) << " secs\n";
	double csTime = elapsed(t1, t2);
	t1 = high_resolution_clock::now();
	for (uint64_t i = 0; i < N; ++i) 
		cm.add(numbers[i]);
	t2 = high_resolution_clock::now();
	std::cout << "CountMin Time to count " << N << " items: " << elapsed(t1, t2) << " secs\n";
	double cmTime = elapsed(t1, t2);

	free(numbers); // free stream

	// Compute heavy hitters
	uint64_t total = 0;
	double threshold = tau * N;
	std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> topK;
	t1 = high_resolution_clock::now();
	for (auto it = map.begin(); it != map.end(); ++it) {
		if (it->second >= threshold) {
			topK.insert(std::make_pair(it->second, it->first));
		}
		total += it->second;
	}
	t2 = high_resolution_clock::now();
	std::cout << "Time to compute phi-heavy hitter items: " << elapsed(t1, t2) << " secs\n";
	assert(total == N);

	std::cout << "Heavy hitter items: \n";
	for (auto it = topK.begin(); it != topK.end(); ++it) {
		std::cerr << "Item: " << it->second << " Count: " << it->first << "\n";
	}
	std::cout << "HT space = " << sizeof(map) + map.size() * sizeof(uint64_t) * 2 << "\n";


	t1 = high_resolution_clock::now();
	std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> res = mg.heavyHitters(tau);
	t2 = high_resolution_clock::now();
	std::set<uint64_t> resSet;
	std::cout << "MG Time to compute phi-heavy hitter items: " << elapsed(t1, t2) << " secs\n";
	for(auto i : res){
		resSet.insert(i.second);
		std::cerr << "Item: " << i.second  << " Count: " << i.first << "\n";
	}
	std::cout << "MG space = " << mg.size() << "\n";

	size_t sum = 0;
	for(auto hh : topK){
		if(resSet.count(hh.second) != 0)
			sum++;
	}
	double mgAcc = (double)sum/(double)res.size();
	double mgRec = (double)sum/(double)topK.size();

	t1 = high_resolution_clock::now();
	std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> csres = cs.heavyHitters(tau);
	t2 = high_resolution_clock::now();
	resSet.clear();
	std::cout << "CS Time to compute phi-heavy hitter items: " << elapsed(t1, t2) << " secs\n";
	for(auto i : csres){
		resSet.insert(i.second);
		std::cerr << "Item: " << i.second  << " Count: " << i.first << "\n";
	}
	std::cout << "CS space = " << cs.size() << "\n";

	sum = 0;
	std::cout << topK.size() << "\n";
	for(auto hh : topK){
		std::cout << hh.second << "\n";
		if(resSet.count(hh.second) != 0)
			sum++;
	}

	double csAcc = (double)sum/(double)csres.size();
	double csRec = (double)sum/(double)topK.size();

	t1 = high_resolution_clock::now();
	std::multimap<uint64_t, uint64_t, std::greater<uint64_t>> cmres = cm.heavyHitters(tau);
	t2 = high_resolution_clock::now();
	resSet.clear();
	std::cout << "CM Time to compute phi-heavy hitter items: " << elapsed(t1, t2) << " secs\n";
	for(auto i : cmres){
		resSet.insert(i.second);
		std::cerr << "Item: " << i.second  << " Count: " << i.first << "\n";
	}
	std::cout << "CM space = " << cm.size() << "\n";

	sum = 0;
	for(auto hh : topK){
		if(resSet.count(hh.second) != 0)
			sum++;
	}

	double cmAcc = (double)((double)sum/(double)cmres.size());
	double cmRec = (double)sum/(double)topK.size();
	std::cout << sum << std::endl;

	std::ofstream out;
	out.open("out.csv", std::ios::app);
	out << 1 << ", " << 1 << ", " << sizeof(map) + map.size() * sizeof(uint64_t) * 2 << ", " << htTime << "\n" ;
	out << mgAcc << ", " << mgRec << ", " << mg.size() << ", " << mgTime <<"\n"; 
	out << csAcc << ", " << csRec << ", " << cs.size() << ", " << csTime <<"\n"; 
	out << cmAcc << ", " << cmRec << ", " << cm.size() << ", " << cmTime <<"\n";
	out.close(); 
	return 0;
}


template class MisraGries<>;
