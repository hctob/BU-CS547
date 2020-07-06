#ifndef BOTTCHER_THREAD_POOL_HPP
#define BOTTCHER_THREAD_POOL_HPP
#include <thread>
#include <future>
#include <pthread.h>
#include <mutex>
#include <vector>
#include <queue>
#include <condition_variable>
#include <functional>
#include <memory>
#include <algorithm>

struct ThreadPool {
	//ThreadPool() : num_threads(std::thread::hardware_concurrency()), handles(num_threads) {}
	inline ThreadPool(const size_t& num_threads) {
		stop = false;
		//lock = new std::mutex;
		for(size_t i = 0; i < num_threads; i++) {
			handles.emplace_back([this] {
					for(;;) {
						std::function<void()> task;
						//printf("Locking...\n");
						//acquire lock for this scope block
						{
							std::unique_lock<std::mutex>lock(this->lock); 
							cond.wait(lock, [this] {
								return (this->stop || !this->task_queue.empty());
							});
							if(stop && task_queue.empty()) return;
						//printf("Popping queue task into thread\n");
							if(!task_queue.empty()) task = std::move(task_queue.front());
							task_queue.pop();
						}
							
						task();

					}
			});
		}
	}

	template<typename F, typename... Args>
	auto enqueue_work(F&& f, Args&& ...args) -> std::future<typename std::result_of<F(Args...)>::type> {
		using ret = typename std::result_of<F(Args...)>::type;
		//https://en.cppreference.com/w/cpp/thread/packaged_task
		//wrap function handle in std::packaged_task so that it can be invoked asynchronously by the thread pool
		auto task = std::make_shared<std::packaged_task<ret()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
		std::future<ret> res = task->get_future(); //assign future in order to receive return value of worker thread

		std::unique_lock<std::mutex> lock(this->lock); //acquire the mutex for this scope block

		//printf("enqueing work\n");
		task_queue.emplace([task] { (*task)(); });

		lock.unlock(); //explictly releasing lock so that I can see where my deadlock issue is arising
		/*task_queue.push([=] {
			f(args...);
		});*/
		cond.notify_one(); //notify the condition variable that
		return res;
	}

	inline ~ThreadPool() {
		//cond.notify_all();
		{
			std::unique_lock<std::mutex>lock(this->lock);
			stop = true;
		}
		cond.notify_all();
		for (auto& t : handles) {
			t.join();
		}
	}
private:
	using Task = std::function<void()>;
	std::vector<std::thread> handles;
	std::queue<Task> task_queue;
	//std::atomic<size_t> num_tasks;
	std::mutex lock;
	std::condition_variable cond;
	bool stop;
};

#endif
