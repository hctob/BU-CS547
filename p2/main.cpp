#include <random>
#include <cstdlib>
#include <cstdio>
#include <errno.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <iostream>
#include <assert.h>
#include <sys/mman.h>
#include <linux/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iomanip>
#include <sys/time.h>
#include <chrono>
#include <sys/resource.h>
#include <thread>
#include "KDTree.hpp"
#include "ThreadPool.hpp"


#define handle_error_en(en, msg) \
	do { errno = en; perror(msg); exit(EXIT_FAILURE);} while (0)

auto comp_char_arr(char* a, char* b, int n) {
	for(auto i = 0; i < n; i++) {
		if(a[i] != b[i]) return false;
	}
	return true;
}

//https://courses.cs.vt.edu/~cs2604/fall00/binio.html#nonchar
//#pragma pack(32)
struct TrainingFile { //
	using ll = unsigned long long;
	char file_name[8]; //8
	ll file_id; //8
	ll num_points; //8
	ll dims; //8
};
//#pragma pack(264)
struct QueryFile { //size is 8 + 4(64)
	using ll = unsigned long long;
	char file_name[8];
	ll file_id;
	ll num_queries; //number of queries - will also be used for number of threads
	ll dims; //dimensions
	ll k; //number of nearest neighbours to return
};
//#pragma pack(392)
struct Result { //size is 8 + 6(64)
	using ll = unsigned long long;
	char file_name[8];
	ll training_file_id; //training file used to generate this result
	ll query_file_id; //query file used to generate this result
	ll result_file_id; //randomly generated
	ll num_queries; //number of queries ran
	ll dims; //dimensions
	ll neighbours; //neighbours

	friend std::ostream& operator<<(std::ostream& os, const Result& res) {
		os << res.file_name << res.training_file_id << res.query_file_id << res.result_file_id << res.num_queries << res.dims << res.neighbours;
		return os;
	}
};

void print_training(const TrainingFile& tf) {
	for(auto i = 0; i < 8; i++)
		printf("%c", tf.file_name[i]);
	std::cout << ", Training ID: " << tf.file_id << ", num_points = " << tf.num_points << ", dimensions = " << tf.dims << std::endl;
}

int main(int argc, char** argv) {
	if(argc < 5 || argc > 5) {
		std::cout << "Usage: ./k-nn-test num_cores training_file query_file result_file" << std::endl;
	}
	else {
		size_t n_cores = static_cast<size_t>(atoi(argv[1]));
		assert(n_cores > 0);
		const char* training_file = argv[2];
		const char* query_file = argv[3];
		const char* result_file = argv[4];
		
		//set CPU affinity 
		int s, j;
		cpu_set_t cpuset;
		pthread_t thread;

		thread = pthread_self();

		CPU_ZERO(&cpuset);
		for(j = 0; j < 24; j++) {
			CPU_SET(j, &cpuset);
		}
		s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
		if(s != 0) {
			handle_error_en(s, "pthread_setaffinity_np");
		}
		TrainingFile tf;
		
		auto start = std::chrono::steady_clock::now();
		int fd = open(training_file, O_RDONLY);
		if (fd < 0) {
			int en = errno;
			std::fprintf(stderr, "Couldn't open %s: %s\n", training_file, strerror(en));
			exit(2);
		}

		struct stat sb;
		int rv = fstat(fd, &sb); assert(rv == 0);
		// std::cout << sb.st_size << std::endl;
		// Make sure that the size is a multiple of the size of a double.
		assert(sb.st_size%sizeof(float) == 0);
		void *vp = mmap(nullptr, sb.st_size, PROT_READ, MAP_PRIVATE|MAP_POPULATE, fd, 0);
		if (vp == MAP_FAILED) {
			int en = errno;
			fprintf(stderr, "mmap() failed: %s\n", strerror(en));
			exit(3);
		}
		rv = madvise(vp, sb.st_size, MADV_SEQUENTIAL|MADV_WILLNEED); assert(rv == 0);

		rv = close(fd); assert(rv == 0);

		//array = (float *) vp;
		char* training_file_pointer = (char*)vp;
		//unsigned long long  n_floats = (sb.st_size - 32)/sizeof(float); //8 chars + 3 unsigned long longs = 32
		//std::cout << "N-floats from m-map: " << n_floats << std::endl;
		for(int i = 0; i < 8; i++) {
			tf.file_name[i] = training_file_pointer[i];
			//printf("%c", tf.file_name[i]);
		}
		memcpy((char*)&tf.file_id, &training_file_pointer[8], 8);
		memcpy((char*)&tf.num_points, &training_file_pointer[16], 8);
		memcpy((char*)&tf.dims, &training_file_pointer[24], 8);
		if(tf.num_points == 0) {
			std::cout << "ERROR: number of training points must be greater than 0.\n";
			exit(1);
		}
		else if (tf.dims == 0) {
			std::cout << "ERROR: training point dimensions must be greater than 0.\n";
			exit(1);
		}
		print_training(tf);

		//reading in floats
		std::vector<std::vector<float>> points;
		//float** points;
		//std::vector<float*> ponts;
		for(size_t i = 0; i < tf.num_points; i++) {
			//points[i] = new float[tf.dims];
			std::vector<float> dims;
			for(size_t j = 0; j < tf.dims; j++) {
				size_t start = 32;
				float val;
				memcpy((char*)&val, &training_file_pointer[start + (sizeof(float) * tf.dims * i) + (sizeof(float)*(j))], sizeof(float));
				//std::cout << points[i][j] << " ";
				dims.push_back(val);
				//std::cout << i << ": " << dims[j] << " ";
			}
			points.push_back(dims);
			//std::cout << std::endl;
		}
		//std::cout << "Number of points: " << points.size() << std::endl;
		//std::cout << "Dimensions: " << points[0].size() << std::endl;

		//time to make a mohootin tree
		size_t num_threads =  std::min((size_t)36, (2 * n_cores) - 1);
		if(num_threads == 0) num_threads = 36;
		//size_t num_threads = (n_cores == 8) ? n_cores : 8;
		//std::cout << "Num_threads: " << num_threads << std::endl;
		//auto start = std::chrono::steady_clock::now();
		auto tree = new cs547::KDTree<float>(points, tf.dims, num_threads); // = new KDTree<float>(points, tf.dims, num_threads);
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> diff = (end - start);
		std::cout << "Build time: " << diff.count() << " seconds."  << std::endl;
		//tree.print_root();
		if(tf.num_points < 128) {
			tree->print_tree();
		}
		// tree.print_tree();

		/*Query this mohootin tree*/
		std::ifstream query(query_file);
		QueryFile qf;
		std::vector<std::vector<float>> query_points;
		
		auto query_start = std::chrono::steady_clock::now();
		if(query) {
			query.read((char*)&qf, sizeof(QueryFile));
			std::cout << "Query information:" << std::endl;
			std::cout << qf.file_name << ", " << qf.file_id << ", " << qf.num_queries << ", " << qf.dims << ", " << qf.k << std::endl;
			for(size_t i = 0; i < qf.num_queries; i++) {
				std::vector<float> point;
				if(i == 0) printf("Query #%d: ", (int)i + 1);
				for(size_t j = 0; j < qf.dims; j++) {
					float val = 0;
					query.read((char*)&val, sizeof(float));
					point.push_back(val);
					if(i == 0) printf("%f ", point[j]);
				}
				if(i == 0) printf("\n");
				query_points.push_back(point);
			}
			query.close();
		}
		if(qf.num_queries == 0) {
			std::cout << "ERROR: cannot run k-nn with 0 queries.\n";
			exit(1);
		}
		else if(qf.dims == 0) {
			std::cout << "ERROR: query points cannot have zero dimensions.\n";
			exit(1);
		}
		else if(qf.k == 0) {
			std::cout << "ERROR: number of neighbors to query cannot be 0.\n";
			exit(1);
		}
		//testing
		/*for(size_t i = 0; i < qf.num_queries; i++) {
		  printf("Query %d: ", (int)i);
		  for(size_t j = 0; j < qf.dims; j++) {
		  std::cout << query_points[i][j];
		  }
		  std::cout << std::endl;
		  }*/
		//
		// cs547::KDTree<float> tree(points, tf.dims, num_threads);
		//tree.print_tree();
		std::cout << std::endl;
		int count = 0;
		/*for(size_t k = 0; k < query_points.size(); k++) {
		
			std::cout << "Querying " << qf.k << " neighbors to point #" << count << ": ";
			for(auto i : query_points[k]) {
				std::cout << i  << " ";
			}
			std::cout << std::endl;
			auto knn = tree->knn(query_points[k], qf.k);

			for(auto i = (int)knn.size() - 1; i >= 0; i--) {
				std::cout << "\tNeighbor #" << i + 1 << ": ";
				for(auto j : knn[i]) {
					std::cout << j << " ";
				}
				std::cout << std::endl;
			}
			count++;
		}*/
		/*printf("\nExecuting queries in parallel:\n");
		//parallelized queries using a threadpool
		auto query_start = std::chrono::steady_clock::now();
		ThreadPool tp(query_points.size());
		count = 0;
		for(size_t i = 0; i < query_points.size(); i++) {
			auto res = tp.enqueue_work(&cs547::KDTree<float>::knn_query, tree, query_points[i], qf.k);
			auto knn = res.get();
			std::cout << "K-Nearest neighbors for ";
			for(auto k : query_points[i]) {
				printf("%f ", k);
			}
			printf("\n");
			count = 0;
			for(auto k = (int)knn.size() - 1; k >= 0; k--) {
				printf("\tNeighbor #%d: ", count + 1);
				for(auto j : knn[k]) {
					printf("%f ", j);
				}
				printf("\n");
				count++;
			}
			//count++;
		}
		auto query_end = std::chrono::steady_clock::now();
		std::chrono::duration<double> q_diff = (query_end - query_start);

		std::cout << "Time to execute queries in parallel: " << q_diff.count() << " seconds." << std::endl;*/
		//tree->print_tree();
		/*Time to generate a basic output file*/
		Result rf;
		//rf.file_name[0] = 'R';
		strcpy(rf.file_name, "RESULT\0\0");
		rf.training_file_id = tf.file_id;
		rf.query_file_id = qf.file_id;
		//todo: pick output id
		//srand (time(NULL));
		//std::default_random_engine generator;
		//std::uniform_int_distribution<int> distribution(1,100000000);
		//int dice_roll = distribution(generator);
		//rf.result_file_id = static_cast<unsigned long long>(rand() % 10000000 + 1);
		//https://stackoverflow.com/questions/2572366/how-to-use-dev-random-or-urandom-in-c/11990066

		FILE* _rand = fopen("/dev/urandom", "r");
		unsigned long long r_id = 0;
		fread(&r_id, 1, sizeof(unsigned long long), _rand);
		fclose(_rand);
		rf.result_file_id = r_id;
		std::cout << "Result file ID from /dev/urandom/: " << r_id << std::endl;
		/*if (_rand < 0) {
		  printf("Error: could not read result_id from /dev/urandom.\n");
		  exit(1); 
		  }    
		  else {
		  unsigned long long r_id = 0;
		  fread(&r_id, sizeof(unsigned long long), 1)
		  }*/
		rf.num_queries = qf.num_queries;
		rf.dims = tf.dims;
		rf.neighbours = qf.k;
		std::ofstream result;
		result.open(result_file, std::ofstream::out | std::ios_base::binary);
		//std::ofstream result(result_file);
		result << rf.file_name;
		result << '\0';
		result << '\0';
		result.write((char*)&rf.result_file_id, sizeof(unsigned long long));
		result.write((char*)&rf.training_file_id, sizeof(unsigned long long));
		result.write((char*)&rf.query_file_id, sizeof(unsigned long long));
		result.write((char*)&rf.num_queries, sizeof(unsigned long long));
		result.write((char*)&rf.dims, sizeof(unsigned long long));
		result.write((char*)&rf.neighbours, sizeof(unsigned long long));
		//result << rf.file_name << rf.training_file_id << rf.query_file_id << rf.result_file_id << rf.num_queries << rf.dims << rf.neighbours;
		
		//auto query_start = std::chrono::steady_clock::now();
		ThreadPool tp(num_threads);
		//std::cout << "Created ThreadPool of size " << query_points.size() << std::endl;
		count = 0;
		for(size_t i = 0; i < query_points.size(); i++) {
			auto res = tp.enqueue_work(&cs547::KDTree<float>::knn_query, tree, query_points[i], qf.k);
			auto knn = res.get();
			if(i == 0) {
				std::cout << "K-Nearest neighbors for ";
				for(auto k : query_points[i]) {
					printf("%f ", k);
				//result << k;
				//result.write((char*)&k, sizeof(float));
				}
			}
			//printf("\n");
			count = 0;
			for(auto k = (int)knn.size() - 1; k >= 0; k--) {
				if(i == 0) printf("\tNeighbor #%d: ", count + 1);
				for(auto j : knn[k]) {
					if(i == 0) printf("%f ", j);
					//result << j;
					result.write((char*)&j, sizeof(float));
				}
				if(i == 0) printf("\n");
				count++;
			}
			//count++;
		}
		auto query_end = std::chrono::steady_clock::now();
		std::chrono::duration<double> q_diff = (query_end - query_start);

		std::cout << "Time to execute queries in parallel: " << q_diff.count() << " seconds." << std::endl;
		//result << rf.file_name << rf.training_file_id << rf.query_file_id;
		result.close();
		delete tree;
		std::cout << "Result file written to " << result_file << std::endl;
	}
	return 0;
}
