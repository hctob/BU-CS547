#include <iostream>
#include <omp.h>
#include <random>
#include <vector>
#include <array>
#include <chrono>
#include <fstream>
#include <assert.h>

const size_t NUM_POINTS = 100'000;

const short MIN_DIMENSION = 2;

//https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates#generating_random_points
auto better_solution(size_t num_points, size_t dims = 16) {
	std::default_random_engine _rand;
	std::normal_distribution<float> uniform(0, 1);
	//will be used for histogram calculations
	//when magnitude of random point is generated, it will be multipled by 
	//100 and used as the index into the array, which will then increment
	//to reflect the number of points in that slice (i.e. 0.01 = index 1)
	std::vector<std::array<size_t, 100>> distributions;
	
	//from 2 dimensions to our specifed max, 16 by default
	for(size_t d = 2; d <= dims; d += 1) {
		//100 slices - represent distances of 0.00 to 1 from origin 
		std::array<size_t, 100> dim_dist;
		for(size_t i = 0; i < num_points; i++) {
			//the radius, or distance from the origin point
			float radius = 0.0;
			std::vector<float> point;
			//generate a n+2 dimensional point
			for(size_t i = 0; i < d + 2; i++) {
				auto p = uniform(_rand);
				radius += std::pow(p, 2);
				point.push_back(p);
			}
			//sqrt the dimensions^2 to get the magnitude
			radius = std::sqrt(radius);
			for(auto& el : point) {
				//normalize the dimensions with the magnitude
				el /= radius;
			}
			//remove the N+1 and N+2 dimensions from the normalized point
			//this guarantees the randomly generated point will be within the surface area of the d-sphere 
			point.pop_back();
			point.pop_back();

			float sum = 0.0;
			for(auto el : point) {
				sum += std::pow(el, 2);
			}
			sum = std::sqrt(sum);
			if(sum >= 1) {
				//if the normalized point we generated has
				//magnitude >= 1, we reject it
				i--; 
			} 
			else {
				//increase the relevant accumlator
				dim_dist[(int)sum * 100] += 1;
			}
		}
		distributions.push_back(dim_dist);		
	}
	return distributions;
	
}



auto better_solution_parallel(size_t num_points, size_t dims = 16, size_t threads = omp_get_max_threads()) {
	std::default_random_engine _rand;
	std::normal_distribution<double> uniform(0, 1);
	//will be used for histogram calculations
	//when magnitude of random point is generated, it will be multipled by 
	//100 and used as the index into the array, which will then increment
	//to reflect the number of points in that slice (i.e. 0.01 = index 1)
	std::vector<std::vector<int>> distributions;
	
	//from 2 dimensions to our specifed max, 16 by default
	#pragma omp parallel for num_threads(threads) shared(distributions)
	for(size_t d = 2; d <= dims; d += 1) {
		//100 slices - represent distances of 0.00 to 1 from origin 
		std::vector<int> dim_dist(100);
		for(size_t i = 0; i < num_points; i++) {
			//the radius, or distance from the origin point
			double radius = 0.0;
			std::vector<float> point;
			//generate a n+2 dimensional point
			for(size_t j = 0; j < d + 2; j++) {
				auto p = uniform(_rand);
				radius += std::pow(p, 2);
				point.push_back(p);
			}
			//sqrt the dimensions^2 to get the magnitude
			radius = std::sqrt(radius);
			for(auto& el : point) {
				//normalize the dimensions with the magnitude
				el /= radius;
			}
			//remove the N+1 and N+2 dimensions from the normalized point
			//this guarantees the randomly generated point will be within the surface area of the d-sphere 
			point.pop_back();
			point.pop_back();

			double sum = 0.0;
			for(auto el : point) {
				sum += std::pow(el, 2);
			}
			sum = std::sqrt(sum); 
			//sum = std::sqrt(sum);
			//increase the relevant accumlator
			++dim_dist[static_cast<size_t>(std::sqrt(sum) * 100)];
		}
		distributions.push_back(dim_dist);		
	}
	return distributions;
	
}

auto reject_v2(size_t num_points, size_t dims = 16, size_t threads = omp_get_max_threads()) {
	std::random_device _rand;
	std::mt19937 gen(_rand());
	std::uniform_real_distribution<double> uniform(-1, 1);
	std::vector<std::vector<size_t>> distributions;
	
	//#pragma omp parallel for num_threads(threads) shared(distributions)
	for(size_t d = 2; d <= dims; d++) {
		std::vector<size_t> dim_hist(100, 0);
		#pragma omp parallel for num_threads(threads) shared(dim_hist)
		for(size_t i = 0; i < num_points; i++) {
			float radius = 0.0;
			std::vector<double> point;
			for(size_t j = 0; j < d; j++) {
				auto num = uniform(gen);
				point.push_back(num);
				radius += (num*num);
			}
			if(radius > 1) continue;
			++dim_hist[radius];
		}
		distributions.push_back(dim_hist);
	}
	return distributions;
}

auto rejection_method(size_t num_points, size_t dims = 16, size_t threads = omp_get_max_threads()) {	
	std::default_random_engine _rand;
	std::uniform_real_distribution<double> uniform(-1, 1);
	//will be used for histogram calculations
	//when magnitude of random point is generated, it will be multipled by 
	//100 and used as the index into the array, which will then increment
	//to reflect the number of points in that slice (i.e. 0.01 = index 1)
	std::vector<std::array<size_t, 100>> distributions;
	
	//from 2 dimensions to our specifed max, 16 by default
	//#pragma omp parallel for num_threads(threads) shared(distributions)
	for(size_t d = 2; d <= dims; d += 1) {
		//100 slices - represent distances of 0.00 to 1 from origin 
		std::array<size_t, 100> dim_dist;
		#pragma omp parallel for num_threads(threads) shared(dim_dist)
		for(size_t i = 0; i < num_points; i++) {
			//the radius, or distance from the origin point
			float radius = 0.0;
			std::vector<float> point;
			//generate a n+2 dimensional point
			for(size_t i = 0; i < d; i++) {
				auto p = uniform(_rand);
				radius += std::pow(p, 2);
				point.push_back(p);
			}
			radius = std::sqrt(radius);
			if (radius > 1.0) {
				i--;
			}
			else {	
				dim_dist[(int)radius * 100] += 1;
			}
		}
		distributions.push_back(dim_dist);
	}
	return distributions;
}

auto inverse_transform_v2(size_t num_points, size_t dims = 16, size_t threads = omp_get_max_threads()) {
	std::random_device _rand;
	std::mt19937 gen(_rand());
	std::uniform_real_distribution<double> uniform(0, 1.0);
	std::vector<std::vector<double>> gen_points(num_points);

	std::vector<std::vector<int>> distributions;
	//#pragma omp parallel for num_threads(threads) shared(distributions)
	for(size_t dim = 2; dim <= dims; ++dim) {
		std::vector<int> dim_dist(100, 0);
		#pragma omp parallel for num_threads(threads) shared(dim_dist)	
		for(size_t i = 0; i < num_points; ++i) {
			double radius = std::pow(uniform(gen), 1/ dim);
			++dim_dist[radius];
		}
		distributions.push_back(dim_dist);
	}
	return distributions;		
}

auto inverse_transform_method(size_t num_points, size_t dims = 16, size_t threads = omp_get_max_threads()) {	
	//std::default_random_engine _rand;
	//std::uniform_real_distribution<double> uniform(0, 1);

	std::random_device _rand;
	std::mt19937 gen(_rand());
	std::uniform_real_distribution<double> uniform(0, 1.0);
	std::vector<std::vector<double>> gen_points(num_points);
	//will be used for histogram calculations
	//when magnitude of random point is generated, it will be multipled by 
	//100 and used as the index into the array, which will then increment
	//to reflect the number of points in that slice (i.e. 0.01 = index 1)
	std::vector<std::vector<size_t>> distributions;
	
	//from 2 dimensions to our specifed max, 16 by default
	#pragma omp parallel for num_threads(threads) shared(distributions)
	for(size_t d = 2; d <= dims; d += 1) {
		//100 slices - represent distances of 0.00 to 1 from origin 
		std::vector<size_t> dim_dist(100, 0);
		for(size_t i = 0; i < num_points; i++) {
			//the radius, or distance from the origin point
			double radius = 0.0;
			std::vector<float> point;
			//generate a n+2 dimensional point
			for(size_t i = 0; i < d; i++) {
				auto p = uniform(gen);
				radius += std::pow(p, 2);
				point.push_back(p);
			}
			//radius = std::sqrt(radius);
			/*for(auto& e : point) {
				e /= radius;
			}*/
			//auto u = uniform(_rand);
			//uniformly sampled number, d-root of u
			radius = std::pow(radius, 1 / d);
			//float new_radius = 0.0;
			/*for(auto& e : point) {
				new_radius += std::pow(e * u, 2);
			}*/
			//new_radius = std::sqrt(new_radius);
			//double new_radius = 1 - radius;
			std::cout << "index = " << 1 - radius << std::endl;		
			dim_dist[(size_t)(radius / 0.01)] += 1;
		}
		distributions.push_back(dim_dist);
	}
	return distributions;
}
/*auto inverse_transform_method_p(size_t num_points, size_t dims = 16, size_t threads = omp_get_max_threads()) {	
	//std::default_random_engine _rand;
	//std::uniform_real_distribution<double> uniform(0, 1);

	std::random_device _rand;
	std::mt19937 gen(_rand());
	std::uniform_real_distribution<double> uniform(0, 1.0);
	std::vector<std::vector<double>> gen_points(num_points);
	//will be used for histogram calculations
	//when magnitude of random point is generated, it will be multipled by 
	//100 and used as the index into the array, which will then increment
	//to reflect the number of points in that slice (i.e. 0.01 = index 1)
	std::vector<std::vector<size_t>> distributions;
	
	//from 2 dimensions to our specifed max, 16 by default
	#pragma omp parallel for num_threads(threads) shared(distributions)
	for(size_t d = 2; d <= dims; d += 1) {
		//100 slices - represent distances of 0.00 to 1 from origin 
		std::vector<size_t> dim_dist(100, 0);
		for(size_t i = 0; i < num_points; i++) {
			//the radius, or distance from the origin point
			double radius = 0.0;
			std::vector<float> point;
			//generate a n+2 dimensional point
			for(size_t i = 0; i < d; i++) {
				auto p = uniform(gen);
				radius += std::pow(p, 2);
				point.push_back(p);
			}
			//radius = std::sqrt(radius);
			for(auto& e : point) {
				e /= radius;
			}
			//auto u = uniform(_rand);
			//uniformly sampled number, d-root of u
			radius = std::pow(radius, 1 / d);
			//float new_radius = 0.0;
			for(auto& e : point) {
				new_radius += std::pow(e * u, 2);
			}
			//new_radius = std::sqrt(new_radius);
			//double new_radius = 1 - radius;
			std::cout << "index = " << 1 - radius << std::endl;		
			dim_dist[(size_t)(radius / 0.01)] += 1;
		}
		distributions.push_back(dim_dist);
	}
	return distributions;
}*/

auto reject_v3(const size_t& num_points, const size_t& dims = 16) {
	std::random_device _rand;
	std::mt19937 gen(_rand());
	std::uniform_real_distribution<double> dist(-1.0, 1.0);
	std::vector<std::vector<double>> gen_points(num_points);
	size_t count = 0;
	size_t loops = 0;
	while(count < num_points) {
		++loops;
		std::vector<double> point(dims); //point of dims dimensiosn
		double radius = 0.0;
		for(size_t i = 0; i < dims; ++i) {
			double dim = dist(gen);
			point[i] = dim;
			radius += (dim * dim);
		}
		radius = std::sqrt(radius);
		if(radius > 1) continue;
		gen_points[count] = point;
		++count;
	}
	return gen_points;
}

auto generate_histogram(const std::vector<std::vector<double>>& points, const size_t& dims) {
	std::vector<size_t> hist(100, 0);
	for(size_t i = 0; i < points.size(); ++i) {
		double radius = 0.0;
		for(size_t j = 0; j < dims; ++j) {
			radius += (points[i][j] * points[i][j]);
		}
		radius = std::sqrt(radius);
		std::cout << "radius before multiply (distance from surface area): " << radius << std::endl;
		radius *= 100;
		std::cout << "index in histogram: " << radius << std::endl;
		++hist[static_cast<size_t>(radius)];	
	}
	
	return hist;
}

int main(int argc, char** argv) {
	if(argc != 5) {
		std::cout << "./d-sphere <num_points> <max_dimension> <threads> <output file name>\n";
		exit(0);
	}
	else {
		size_t num_points = static_cast<size_t>(atoi(argv[1]));
		size_t dim = static_cast<size_t>(atoi(argv[2]));
		
		size_t threads = static_cast<size_t>(atoi(argv[3]));
		auto output_name = argv[4];
		/*auto start = std::chrono::system_clock::now();
		auto histogram = better_solution(num_points, dim);
		auto stop = std::chrono::system_clock::now();
		std::chrono::duration<double> diff = (stop - start);
		std::cout << "Sequential: " << diff.count() << " seconds.\n";

		auto startp = std::chrono::system_clock::now();
		auto histogramp = better_solution_parallel(num_points, dim, threads);
		auto stopp = std::chrono::system_clock::now();
		std::chrono::duration<double> diffp = (stopp - startp);
		std::cout << "Parallel: " << diffp.count() << " seconds.\n";
		*/

		auto start = std::chrono::system_clock::now();
		auto histogram2 = better_solution(num_points, dim);
		//auto points = reject_v3(num_points, dim);
		//auto hist = generate_histogram(points, dim);
		auto stop = std::chrono::system_clock::now();
		std::chrono::duration<double> diff = (stop - start);
		std::cout << "extra credit method (seq): " << diff.count() << " seconds.\n";
		auto start1 = std::chrono::system_clock::now();
		auto histogram = better_solution_parallel(num_points, dim, threads);
		//auto points = reject_v3(num_points, dim);
		//auto hist = generate_histogram(points, dim);
		auto stop1 = std::chrono::system_clock::now();
		std::chrono::duration<double> diff2 = (stop1 - start1);
		std::cout << "extra credit method (parallel): " << diff2.count() << " seconds.\n";
		/*const auto interval = 0.01;
		for(size_t i = 0; i < hist.size(); ++i) {
			std::cout << interval * i << "," << 1.0 * hist.at(i) / num_points * 100.00 << "%\n";
		}*/
		std::ofstream os(output_name);
		for(size_t i = 0; i < histogram.size(); i++) {
			const auto interval = 0.01;
			std::cout << "dim = " <<  i+2 << "\n";
			os << "dim = " << i + 2 << "\n";
			for(size_t slice = 0; slice < histogram[i].size(); slice += 1) { 
				std::cout << interval * slice << ", " << 1.0 * histogram[i].at(slice) / num_points * 100 << "%\n"; 
				os << interval * slice << ", " << 1.0 * histogram[i].at(slice) / num_points * 100 << "%\n"; 
				//std::cout << interval * slice << ", " << histogram[i][slice] / num_points * 100 << "\n"; 
			}
		}
	}
	return 0;
}

