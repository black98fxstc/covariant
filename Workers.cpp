#include "Workers.hpp"

// Explicit static member instantiations
std::queue<Work *> Worker::work_list;
volatile bool Worker::kiss_of_death = false;
std::mutex Worker::serialize;
std::mutex Worker::mutex;
std::condition_variable Worker::work_available;

Work::~Work() = default;

void Worker::work() noexcept
{
    Work *work = dequeue();
    if (!work) // spurious wakeup
        return;

    work->parallel();
    {
        std::lock_guard<std::mutex> lock(serialize);
        work->serial();
    }

    delete work;
}

Work *Worker::dequeue() noexcept
{
    std::lock_guard<std::mutex> lock(mutex);
    if (work_list.empty())
        return nullptr;
    Work *work = work_list.front();
    work_list.pop();
    return work;
}

bool Worker::idle()
{
    std::lock_guard<std::mutex> lock(mutex);
    return work_list.empty();
}

void Worker::wait()
{
    std::unique_lock<std::mutex> lock(mutex);
    while (work_list.empty() && !kiss_of_death)
        work_available.wait(lock);
}

void Worker::enqueue(Work *work) noexcept
{
    {
        std::lock_guard<std::mutex> lock(mutex);
        work_list.push(work);
    }
    work_available.notify_one();
}

void Worker::kiss() noexcept
{
    {
        std::lock_guard<std::mutex> lock(mutex);
        kiss_of_death = true;
    }
    work_available.notify_all();
}

void Worker::revive() noexcept
{
    std::lock_guard<std::mutex> lock(mutex);
    while (!work_list.empty())
    {
        delete work_list.front();
        work_list.pop();
    }
    kiss_of_death = false;
}

Worker::Worker(bool threaded) noexcept
{
    if (threaded)
        while (!kiss_of_death)
            if (idle())
                wait();
            else
                work();
    else
        while (!idle())
            work();
}