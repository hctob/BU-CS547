#include <immintrin.h>
#include <iostream>
#include <cmath>
#include <functional>
#include <chrono>
#include <cassert>
#include <random>

const int N = 16*1'000'000; //16 million boy

double time(const std::function<void()> &f) {
	f(); //warmup call
	auto start = std::chrono::system_clock::now();
	f();
	auto stop = std::chrono::system_clock::now();
	//return the runtime of the algorithm in seconds
	return std::chrono::duration<double>(stop - start).count();
}

int main() {
	alignas(32) static float x[N], y[N], z[N], t[N], x1[N], y1[N], z1[N], t1[N];
	
	//generate data for line segments
	std::default_random_engine r; //random engine 
	std::uniform_real_distribution<float> dist(-1, 1);

	/*sequential distance algorithm*/

	for(auto i = 0; i < N; i++) {
		//first line segment in 4D
		x[i] = dist(r);
		y[i] = dist(r);
		z[i] = dist(r);
		t[i] = dist(r);
		//second line segment in 4D
		x1[i] = dist(r);
		y1[i] = dist(r);
		z1[i] = dist(r);
		t1[i] = dist(r);
	}

	static float l_s[N];
	auto seq = [&]() {
		for(auto i = 0; i < N; i++) {
			const float x2 = std::pow(x[i] - x1[i], 2);
			const float y2 = std::pow(y[i] - y1[i], 2);
			const float z2 = std::pow(z[i] - z1[i], 2);
			const float t2 = std::pow(t[i] - t1[i], 2);
			l_s[i] = std::sqrt(x2 + y2 + z2 + t2);
		}
	};

	std::cout << "Line segment (sequential): " << (N/time(seq)) /1000000 << " Mops/s" << std::endl;

	//Data parallelism - SIMD instructions
	//
	alignas(32) static float l_v[N];

	auto parallel = [&]() {
		for(int i = 0; i < N/8; i++) {
			//load values into first line point
			__m256 ymm_x = _mm256_load_ps(x + 8*i);
			__m256 ymm_y = _mm256_load_ps(y + 8*i);
			__m256 ymm_z = _mm256_load_ps(z + 8*i);
			__m256 ymm_t = _mm256_load_ps(t + 8*i);

			//load values into second line point
			__m256 ymm_x1 = _mm256_load_ps(x1 + 8*i);
			__m256 ymm_y1 = _mm256_load_ps(y1 + 8*i);
			__m256 ymm_z1 = _mm256_load_ps(z1 + 8*i);
			__m256 ymm_t1 = _mm256_load_ps(t1 + 8*i);
			
			//compute difference between two points
			__m256 ymm_x2 = _mm256_sub_ps(ymm_x, ymm_x1);
			__m256 ymm_y2 = _mm256_sub_ps(ymm_y, ymm_y1);
			__m256 ymm_z2 = _mm256_sub_ps(ymm_z, ymm_z1);
			__m256 ymm_t2 = _mm256_sub_ps(ymm_t, ymm_t1);
			//sqrt for euclidean distance
			__m256 ymm_sqrt = _mm256_sqrt_ps(_mm256_mul_ps(ymm_x2, ymm_x2) + _mm256_mul_ps(ymm_y2, ymm_y2) + _mm256_mul_ps(ymm_z2, ymm_z2) + _mm256_mul_ps(ymm_t2, ymm_t2));
			_mm256_store_ps(l_v + 8*i, ymm_sqrt);
		}
	};

	std::cout << "Line segment (parallel): " << (N/time(parallel)) / 1000000 << " Mops/s" << std::endl;

	for(int i = 0; i < N; i++) {
		if(std::abs(l_s[i] - l_v[i]) > 0.001) {
			assert(false);
		}
	}
}
