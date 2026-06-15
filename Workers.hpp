#pragma once

// derived from
/*
 * Developer: Wayne Moore <wmoore@stanford.edu> 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * License: BSD 3 clause
 */
// abstract class representing a unit of work to be done
// virtual functions let subclasses specialize tasks
#include <queue>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <future>
#include <functional>
#include <type_traits>

class Work
{
public:

    // many threads can execute this in parallel
    virtual void parallel() noexcept = 0;

    // then only one thread at a time can run
    virtual void serial() noexcept = 0;

    virtual ~Work();

protected:
    explicit Work() noexcept {};
};

// a generic worker thread. looks for work, does it, deletes it
// virtual functions in the work object do all the real work
class Worker
{
private:
    static std::queue<Work *> work_list;
    static std::atomic<bool> kiss_of_death;
    std::thread th;

protected:
    static std::mutex serialize;
    static std::mutex mutex;
    static std::condition_variable work_available;

    void work(bool blocking = true) noexcept;
    Work *dequeue(bool blocking) noexcept;
    bool idle();

public:
    static void enqueue(Work *work) noexcept;
    static void kiss() noexcept;
    static void revive() noexcept;

    Worker(bool threaded = true) noexcept;
    ~Worker();

    Worker(Worker&& other) noexcept;
    Worker& operator=(Worker&& other) noexcept;
};

template <typename T>
class ThreadSafeQueue {
private:
    std::queue<T> queue_;
    std::mutex mutex_;
    std::condition_variable cv_;
    bool stopped_ = false;

public:
    void push(T item) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            queue_.push(std::move(item));
        }
        cv_.notify_one();
    }

    // Returns false if the queue is stopped and empty
    bool pop(T& item) {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this]() { return !queue_.empty() || stopped_; });

        if (queue_.empty() && stopped_) {
            return false;
        }

        item = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    void stop() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopped_ = true;
        }
        cv_.notify_all();
    }
};

class ThreadPool {
private:
    std::vector<std::thread> workers;
    ThreadSafeQueue<std::function<void()>> tasks;

public:
    explicit ThreadPool(size_t threads) {
        for(size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this] {
                std::function<void()> task;
                // Threads sleep here until a task is pushed or the queue stops
                while(tasks.pop(task)) {
                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        tasks.stop();
        for(std::thread &worker: workers) {
            if(worker.joinable()) worker.join();
        }
    }

    // Enqueue a function and return a std::future containing its return value
    // Modern C++: use lambda captures instead of std::bind for passing arguments
    template<class F>
    auto enqueue(F&& f) -> std::future<typename std::invoke_result<F>::type> {
        using return_type = typename std::invoke_result<F>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(std::forward<F>(f));
        std::future<return_type> res = task->get_future();
        tasks.push([task]() { (*task)(); });
        return res;
    }
};
