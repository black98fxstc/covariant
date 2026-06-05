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
class Worker : public std::thread
{
private:
    static std::queue<Work *> work_list;
    volatile static bool kiss_of_death;

protected:
    static std::mutex serialize;
    static std::mutex mutex;
    static std::condition_variable work_available;

    void work() noexcept;
    Work *dequeue() noexcept;
    bool idle();
    void wait();

public:
    static void enqueue(Work *work) noexcept;
    static void kiss() noexcept;
    static void revive() noexcept;

    Worker(bool threaded = true) noexcept;
};
