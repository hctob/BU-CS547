#ifndef BOTTCHER_K_D_TREE_HPP
#define BOTTCHER_K_D_TREE_HPP

#include <atomic>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>
#include <vector>
#include <initializer_list>
#include <queue>
#include <thread>
#include <future>
#include <map>
#include <memory>
#include <limits>
#include <pthread.h>

namespace cs547 {
	size_t COUNT = 20;
	template<typename T>
		struct KDTree {
			struct Node {
				int axis;
				int depth;
				//std::vector<Point> points;
				//Point point;
				//std::array<T, k> point;
				std::vector<T> point;
				Node* parent;
				Node* left_child;
				Node* right_child;
				Node() = delete;
				explicit Node(std::vector<T>& point, Node* left_child, Node* right_child) : point(std::move(point)), left_child(left_child), right_child(right_child) {}
				/*explicit Node(std::vector<T> point, Node* left_child, Node* right_child) {
				  this->point = point;
				  this->left_child = left_child;
				  this->right_child = right_child;
				  }*/

				explicit Node(const int& axis, const int& depth, std::vector<T>& point) : axis(axis), depth(depth), point(point) {
					parent = nullptr;
					left_child = nullptr;
					right_child = nullptr;
				}
				Node(const int& axis, const int& depth, std::vector<T> point, Node* parent, Node* left, Node* right) {
					this->axis = axis;
					this->depth = depth;
					this->point = point;
					this->parent = parent;
					this->left_child = left;
					this->right_child = right;
				}
				Node(const int& axis, const int& depth, std::vector<T> point, Node* parent) : axis(axis), depth(depth), point(point), parent(parent), left_child(nullptr), right_child(nullptr) {
				}
				//Node(const int& axis, const std::vector<Point> points) : axis(axis), points(points), left_child(nullptr), right_child(nullptr) {}
				//Node(const int& axis, const std::array<T, k>& point): axis(axis), point(point), left_child(nullptr), right_child(nullptr) {}
				//Node(const int& axis, const std::vector<T>& point): axis(axis), point(point), left_child(nullptr), right_child(nullptr) {}
				friend std::ostream& operator<<(std::ostream& os, const Node& p)
				{
					os << '(';
					for (size_t i = 0; i < p.point.size(); ++i)
					{
						if (i > 0)
							os << ", ";
						os << p.point[i];
					}
					os << ')';
					os << std::endl;
					return os;
				}
				//returns the euclidean distance from this point to another point
				constexpr float distance(const Node& n) {
					assert(point.size() == n.point.size());
					float dist = 0.0f;
					for(size_t i = 0; i < point.size(); i++) {
						dist += std::pow((point[i] - n.point[i]), 2);
					}
					return std::sqrt(dist);
				}
				constexpr double point_distance(std::vector<float> f) { //returns distance from point f
					assert(point.size() == f.size());
					double dist = 0.0f;
					for(size_t i = 0; i < f.size(); i++) {
						dist += std::pow((point[i] - f[i]), 2);
						//todo: add overflow exception
						/*if(dist > std::numeric_limits<double>::max()) {
							std::cout << "Overflow error: ";
							exit(1);
						}*/
					}
					//return std::sqrt(dist);
					return dist;
				}

				constexpr double pole_distance(std::vector<float> f, int axis) {
					assert(point.size() == f.size());
					double dist = 0.0f;
					dist += std::pow(f[axis] - point[axis], 2);
					return dist;
				}
			};

			using node_dist = std::pair<Node*, double>;

			Node* cmp_node(std::vector<std::vector<float>> point, Node* x, Node* y) {
				if(x == nullptr) return y;
				else if(y == nullptr) return x;
				else {
					auto dist1 = x->point_distance(point);
					auto dist2 = y->point_distance(point);

					if(dist1 < dist2) {
						return x;
					}
					else {
						return y;
					}
				}

			}
			//public data of the KDTree
			Node* root; //root of the tree - median of points split along first axisi
			Node* closest; //closest possible node
			//std::vector<Node*> node_pointers; //vector of node pointers for eas(ier) sorting, hopefully
			size_t dimensions;
			std::vector<std::vector<T>> _points;
			std::vector<Node*> nodes;
			size_t max_threads;
			std::atomic<size_t> num_threads;
			std::atomic<size_t> num_nodes;
			std::map<float, std::vector<float>> distance_map;
			//using nodes = std::vector<Node>;
			//std::unique_ptr<Node> best;
			KDTree() = delete;
			KDTree(std::vector<std::vector<T>> points, size_t dimensions, size_t max_threads) : root(nullptr), dimensions(dimensions), max_threads(max_threads) {
				//printf("Constructor: # of points = %llu\n", points.size());
				_points = points;
				nodes.reserve(points.size());
				//max_threads = max_threads;
				num_threads = 0;
				//root = thread_ctor(points, 0, root);
				//root = _ctor(points, 0, root); //working non-parallel build
				std::promise<Node*>* root_promise = new std::promise<Node*>;
				//distance_map = new map<float, std::vector<float>>();
				//size_t num_threads = 0;
				//root = ctor(points,0, num_threads, root_promise); //currently buggy parallel build
				//std::atomic<size_t> num_threads(0);
				num_threads = 0;
				num_nodes = 0;
				root = thread_ctor(points, 0, root_promise);
				//root->parent = root;
				//	delete root_promise;
				std::cout << "Number of nodes constructed: " << num_nodes << std::endl;
				delete root_promise;
			}
			~KDTree() {
				destroy(root);
			}
			//returns pointer to the Node at which the point currently resides
			Node* search(Node* n, const std::vector<float>& point) const {
				if(n == nullptr || n->point == point) {
					return n;
				}
				else {
					std::vector<float> curr(point);
					int size = (int)point.size();
					int depth = n->depth;
					if(point[depth % size] < curr[depth % size]) {
						return (n->left_child == nullptr) ? curr : search(n->left_child, point);
					}
					else {
						return (n->right_child == nullptr) ? curr : search(n->right_child, point);
					}
				}
			}

			bool exists(const std::vector<float>& point) const {
				return (search(root, point) == nullptr) ? false : true;
			}

			constexpr double minkowski_distance(const std::vector<float>& a, const std::vector<float>& b) {
				assert(a.size() == b.size());
				double dist = 0.0;
				for(size_t i = 0; i < a.size(); i++) {
					dist += std::pow(a[i] - b[i], a.size());
				}

				return std::pow(dist, 1 / a.size());

			}

			std::vector<T> nearest(std::vector<T> query) {
				std::vector<T> mins(query.size());
				for(auto i = 0; i < query.size(); i++) {
					mins[i] = 99999999;
				}
				for (auto& p : _points) {
					mins = std::min(mins, p);
				}

				return mins;
			}
			//inserts point into tree, partition based on dimensions of point
			Node* insert(Node* root, const std::vector<float>& point) {
				return _insert(root, point, 0);
			}
			void print_root() {
				this->print_node(this->root);
			}
			friend std::ostream& operator<<(std::ostream& os, const KDTree& kd) {
				print_tree();
				return os;
			}
			// Wrapper over print2DUtil()
			void print_tree()
			{
				// Pass initial space count as 0
				print2DUtil(this->root, 0);
			}
			static std::vector<std::vector<float>> knn_query(KDTree* kd, std::vector<T> query, size_t k) {
				return kd->knn(query, k);
			}
			std::vector<std::vector<float>> knn(const std::vector<T>& query, size_t k) {
				//todo: priority queue
				//if(root == nullptr) return nullptr;
				//
				//distanceMap[distance] = point
				auto dist_queue = new std::priority_queue<node_dist, std::vector<node_dist>, node_dist_cmp>(); //should be std::less by default
				//auto dist_queue = std::priority_queue<float, std::vector<float>>
				std::vector<std::vector<float>> ret;
				ret.reserve(k);

				knn_logic(root, query, 0, k, dist_queue);
				for(size_t j = 0; j < k && dist_queue->empty() == false; j++) {
					//pushing minimum point values into return vector
					ret.push_back(dist_queue->top().first->point); //return point from within Node* of smallest distance
					//std::cout << "New k-n value, j = " << j << std::endl;
					dist_queue->pop();
				}
				delete dist_queue;
				return ret; //return set of k nearest neighbors to query point
			}
			std::vector<std::vector<float>> parallel_knn(const std::vector<T>& query, size_t k) {
				//todo: priority queue
				//if(root == nullptr) return nullptr;
				//
				//distanceMap[distance] = point
				auto dist_queue = new std::priority_queue<node_dist, std::vector<node_dist>, node_dist_cmp>(); //should be std::less by default
				//auto dist_queue = std::priority_queue<float, std::vector<float>>
				std::vector<std::vector<float>> ret;
				ret.reserve(k);

				parallel_query(root, query, 0, k, dist_queue);
				for(size_t j = 0; j < k && dist_queue->empty() == false; j++) {
					//pushing minimum point values into return vector
					ret.push_back(dist_queue->top().first->point); //return point from within Node* of smallest distance
					//std::cout << "New k-n value, j = " << j << std::endl;
					dist_queue->pop();
				}
				delete dist_queue;
				return ret; //return set of k nearest neighbors to query point
			}

private:
			class node_dist_cmp {
				public:
				bool operator() (const node_dist& a, const node_dist& b) const {
					return a.second < b.second;
				}
			};

		void parallel_query(Node* root, const std::vector<float>& query,
			size_t depth, size_t k, std::priority_queue<node_dist, std::vector<node_dist>, node_dist_cmp>* dist_queue) {
				if(root == nullptr) { return; };
				auto d = depth % dimensions;
				dist_queue->push(std::make_pair(root, root->point_distance(query)));
		}

		void knn_logic(Node* root, const std::vector<float>& query,
				size_t depth, size_t k,  std::priority_queue<node_dist, std::vector<node_dist>, node_dist_cmp>* dist_queue){
			//todo: k-nn checking other inverted branches
			if (root == nullptr) { return; }
			size_t d = depth % dimensions;
			dist_queue->push(std::make_pair(root, root->point_distance(query)));  //using euclidean distance

			//dist_queue->push(std::make_pair(root, minkowski_distance(root->point, query)));
			//std::cout << "Queue size: " << dist_queue->size() << std::endl;
			//std::cout << "Queue top: " << dist_queue->top().first->point[0] << std::endl;
			//printf("Distance from %f to query [%f]: %f\n", root->point[0], query[0], root->point_distance(query));
			if(dist_queue->size() > k) {
				dist_queue->pop(); //pop element with greatest distance from query point to maintain k neighbors
			}
			//todo: logic to check in other branches for lesser values
			if(root->point.at(d) > query.at(d)) {
				//go right
				knn_logic(root->right_child, query, depth + 1,k, dist_queue);
				//todo: prune away tree results
				node_dist furthest_k = dist_queue->top();
				auto root_perp_dist = root->pole_distance(query, d);
				//printf("node: %f, furthest k in queue: %f\n", root->point.at(d), furthest_k.first->point.at(d));
				if(fabs(root_perp_dist) <= fabs(furthest_k.second)) {
					//point on other side of parition wall is closer than the current furthest k, must prune away subtree by going down other branch
					//printf("Going towards %f", root->left_child->point[0]);
					knn_logic(root->left_child, query, depth + 1,k, dist_queue);
				}
			}
			else {
				knn_logic(root->left_child, query, depth + 1,k, dist_queue);
				//todo: prune away tree results
				node_dist furthest_k = dist_queue->top();
				auto root_perp_dist = root->pole_distance(query, d); //distance from query of point perpendicular to the root node
				if(fabs(root_perp_dist) <= fabs(furthest_k.second)) {
					//go down other branch of the tree in order to check for other viable neighbors
					knn_logic(root->right_child, query, depth + 1, k, dist_queue);
				}
			}
		}

		/*void* ctor_wrapper(void* n) {
		  BuildInfo* bi = (BuildInfo*)n;
		  auto node = ken_thread_ctor(bi->points, bi->depth, bi->parent);
		  delete bi;
		  return node;
		  }*/
			public:
		//median of medians stuff
		float median(std::vector<float> point) {
			return point.at(point.size() / 2);
		}

		float median_of_medians(std::vector<std::vector<float>> points) {
			std::vector<float> medians;
			for(size_t i = 0; i < points.size(); i++) {
				auto med = median(points[i]);
				medians.push_back(med);
			}
			return median(medians);
		}

		/*float select_mom() {
		}*/
		Node* thread_ctor(std::vector<std::vector<T>>& points, int depth, std::promise<Node*>* promise) {
			//size_t num_threads_ = num_threads;
			//std::atomic<int> num_threads = new std::atomic<int>(0);
			if(points.size() == 1) {
				int d = depth % dimensions;
				Node* n = new Node{d, depth, points[0]};
				//nodes.push_back(n);
				num_nodes++;
				if(promise) promise->set_value_at_thread_exit(n);
				return n;
			}
			else if(points.size() == 0) {
				if(promise) promise->set_value_at_thread_exit(nullptr);
				return nullptr;
			}

			int axis = depth % dimensions;
			size_t median = points.size() / 2;
			/*std::sort(points.begin(), points.end(), [&](std::vector<T> a, std::vector<T> b) {
					return a.at(axis) < b.at(axis);
					});*/
			std::nth_element(points.begin(), points.begin() + median, points.end()); //lambda comparison doubles build time, same accuracy without, might be more accurate though?

			//begin parallelism
			std::vector<std::thread> handles;
			//handles.reserve(max_threads);
			//std::vector<std::thread>* handles = new std::vector<std::thread>();
			std::promise<Node*> left_p;
			std::promise<Node*> right_p;

			auto left_future = left_p.get_future();
			auto right_future = right_p.get_future();
			bool left_processing = false;
			bool right_processing = false;

			std::vector<std::vector<T>> l_slice{points.begin(), points.begin() + median};
			std::vector<std::vector<T>> r_slice{points.begin() + median + 1, points.end()};

			//create threads
			if(num_threads < max_threads) {
				num_threads++;
				//std::cout << "Thread #" << num_threads << " created for left\n";
				//printf("Created thread #%d\n", num_threads);
				//printf("Processing left in parallel\n");
				//std::vector<std::vector<T>> l_slice{points.begin(), points.begin() + median};
				handles.emplace_back(&KDTree::thread_ctor, this, std::ref(l_slice), depth + 1, &left_p);
				left_processing = true;
			}
			if(num_threads < max_threads) {
				num_threads++;
				//std::cout << "Thread #" << num_threads << " created for right\n";
				//printf("Processing right in parallel\n");
				//std::vector<std::vector<T>> r_slice{points.begin() + median + 1, points.end()};
				handles.emplace_back(&KDTree::thread_ctor, this, std::ref(r_slice), depth + 1, &right_p);
				right_processing = true;
			}
			Node* left_child = nullptr;
			Node* right_child = nullptr;

			if(!left_processing) {
				//no thread is currently processing left subtree. do it normally
				//std::cout << "Processing left normally\n";
				//std::vector<std::vector<T>> l_slice{points.begin(), points.begin() + median};
				left_child = thread_ctor(l_slice, depth + 1, nullptr);
			}
			if(!right_processing) {
				//std::cout << "Processing right normally\n";
				//std::vector<std::vector<T>> r_slice{points.begin() + median + 1, points.end()};
				right_child = thread_ctor(r_slice, depth + 1, nullptr);
			}
			for(auto& i : handles) {
				if (i.joinable()) {
					i.detach();
				}
			}
			handles.clear();

			if(left_processing) {
				left_child = left_future.get(); //either return new left subtree or block until it returns
				//num_threads--;
			}
			if(right_processing) {
				right_child = right_future.get(); //either return new right subtree or block until it returns
				//num_threads--;
			}
			auto node = new Node{axis, depth, points[median]};
			//node->axis = axis;
			//node->depth = depth;
			//node->point = points[median];
			node->left_child = left_child;
			node->right_child = right_child;
			//nodes.push_back(node);
			num_nodes++;
			//printf("new node created\n");
			if(promise) {
				promise->set_value_at_thread_exit(node);
				return node;
			}
			return node;

		}
			private:
		Node* _ctor(std::vector<std::vector<T>>& points, int depth, Node* parent) {
			//std::cout << "Number of points passed into tree: " << points.size() << std::endl;
			//https://en.wikipedia.org/wiki/K-d_tree#Construction
			//https://gopalcdas.com/2017/05/24/construction-of-k-d-tree-and-using-it-for-nearest-neighbour-search/
			//if(points.size() == 0) return nullptr;
			if(points.size() == 1) {
				int d = depth % dimensions;
				Node* n = new Node{d, depth, points[0], parent};
				nodes.push_back(n);
				return n;
			}
			else if (points.size() == 0) {
				return nullptr;
			}

			int axis = depth % dimensions; //axis along which to split partition - chosen in round robin fashion
			size_t median = points.size() / 2;
			//https://stackoverflow.com/questions/1719070/what-is-the-right-approach-when-using-stl-container-for-median-calculation
			//std::nth_element(points.begin(), points.begin() + median, points.end());
			std::sort(points.begin(), points.end(), [&](std::vector<T> a, std::vector<T> b) {
					return a.at(axis) < b.at(axis);
					});
			//if(depth == 0) {
			//	std::sort(points.begin(), points.end());
			//}
			Node* node = new Node{axis, depth, points.at(median), parent};
			nodes.push_back(node);
			//node->axis = axis;
			std::vector<std::vector<T>> left_slice{points.begin(), points.begin() + median}; //O(N)
			std::vector<std::vector<T>> right_slice{points.begin() + median + 1, points.end()}; //O(N)

			node->left_child = _ctor(left_slice, depth + 1, node); //possibly divide amount of threads, likely will have to
			node->right_child = _ctor(right_slice, depth + 1, node);
			return node;
		}
		void destroy(Node* node) {
			if (node)
			{
				destroy(node->left_child);
				destroy(node->right_child);
				delete node;
			}
		}

		Node* _insert(Node* root, const std::vector<float>& point, int depth) {
			int axis = depth % dimensions;
			if(root == nullptr) {
				return new Node(axis, depth, point, root);
			}
			else {
				if(point[axis] < root->point[axis]) {
					root->left_child =  insert(root->left_child, point, depth + 1);
				}
				else {
					root->right_child =  insert(root->right_child, point, depth + 1);
				}
			}
			return root;
		}

		//nearest neighbour garbage
		void recursive_nearest(Node* n, const std::vector<float>& query, std::priority_queue<T>& queue) {
			double distance = n->distance(query);
			//queue.enqueue(distance, );
		}
		void print_node(Node* n) {
			if(n == nullptr) {
				std::cout << "Error: empty node";
			}
			else
				std::cout << *n;
		}
		//vertical printing of tree
		void print2DUtil(Node *root, int space)
		{
			// Base case
			if (root == NULL)
				return;

			// Increase distance between levels
			space += COUNT;

			// Process right child first
			print2DUtil(root->right_child, space);

			// Print current node after space
			// count
			std::cout<<std::endl;
			for (int i = COUNT; i < space; i++)
				std::cout<<" ";
			std::cout<<*root<<"\n";

			// Process left child
			print2DUtil(root->left_child, space);
		}
		};
}
#endif
