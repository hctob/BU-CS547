#include <iostream>
#include <vector>
#include <random>
#include <atomic>
#include <stdlib.h>
#include <thread>
#include <pthread.h>
#include <stdio.h>
#include <algorithm>
#include <math.h>
#include <errno.h>

double func1(const double& x) {
	double num = sin(x);
	double denom = x;
	return static_cast<double>((double)num / (double)denom);
}

std::atomic<double> total(0);
//std::atomic<double> total = std::atomic<double>(0);

template<typename T>
T atomic_fetch_add(std::atomic<T> &obj, T arg) {
	T expected = obj.load();
	while(!atomic_compare_exchange_weak(&obj, &expected, expected + arg))
		;
	return expected;
}
int main(int argc, char** argv) {
	if (argc < 5) {
		printf("Error: ./integrate a b samples threads");
		exit(1);
	}
	else {
		double a = static_cast<double>(atof(argv[1]));
		double b = static_cast<double>(atof(argv[2]));
		size_t n = static_cast<size_t>(atoi(argv[3]));
		int num_threads = static_cast<int>(atoi(argv[4]));
		total = total - total;
		size_t workload = (size_t)(n / num_threads);
		std::vector<std::thread> threads;
		for (auto i = 0; i < num_threads; ++i) {
			threads.push_back(std::thread([&] {
				std::default_random_engine rand;
				std::uniform_real_distribution<double> rand_gen(a, b);
				size_t i = 0;
				
				double thread_res = 0.0;
				do {
					i++;
					double x = rand_gen(rand);
					double val = func1(x);
					thread_res += val;

				} while(i < workload);
				//breaks here
				atomic_fetch_add<double>(total, thread_res);
			}));	
		}
		for (auto i = 0; i < num_threads; ++i) {
			threads[i].join();
		}
		double res = total;
		res *= (b - a) / (double)n;
		printf("%.10f\n", static_cast<double>(res));
	}
	return 0;
}
