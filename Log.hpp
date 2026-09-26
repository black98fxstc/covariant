#pragma once

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>

// say and log format into a private buffer. The temporary's destructor appends
// that finished text under one mutex, so worker threads may format at the same
// time and a message still arrives whole in the session file and, when enabled,
// on the console.

enum class LogLevel
{
    Say,
    Log
};

class SessionLog;

class LogLine
{
public:
    LogLine(SessionLog *session, LogLevel level);
    LogLine(const LogLine &) = delete;
    LogLine &operator=(const LogLine &) = delete;
    LogLine(LogLine &&other) noexcept;
    LogLine &operator=(LogLine &&) = delete;
    ~LogLine();

    template <typename T>
    LogLine &operator<<(T &&value)
    {
        buf_ << std::forward<T>(value);
        return *this;
    }

    LogLine &operator<<(std::ostream &(*manip)(std::ostream &));
    LogLine &operator<<(std::ios_base &(*manip)(std::ios_base &));

private:
    SessionLog *session_;
    LogLevel level_;
    std::ostringstream buf_;
    bool active_ = true;
};

class LogChannel
{
public:
    LogChannel(SessionLog *session, LogLevel level) : session_(session), level_(level) {}

    template <typename T>
    LogLine operator<<(T &&value) const
    {
        LogLine line(session_, level_);
        line << std::forward<T>(value);
        return line;
    }

    LogLine operator<<(std::ostream &(*manip)(std::ostream &)) const;
    LogLine operator<<(std::ios_base &(*manip)(std::ios_base &)) const;

private:
    SessionLog *session_;
    LogLevel level_;
};

class SessionLog
{
    friend class LogLine;

public:
    // Creates directory and a session file named with the local ISO-8601 time.
    // Colons in the clock time are written as hyphens so the name is legal on Windows.
    bool open(const std::filesystem::path &directory);

    // Always appended to the session file and written to stderr.
    void error(const std::string &text);

    void set_quiet(bool quiet);
    void set_verbose(bool verbose);

    bool is_open() const;
    std::filesystem::path path() const;

private:
    void write(LogLevel level, const std::string &text);

    mutable std::mutex mutex_;
    std::ofstream file_;
    std::filesystem::path path_;
    bool quiet_ = false;
    bool verbose_ = false;
};

inline LogLine::LogLine(SessionLog *session, LogLevel level)
    : session_(session), level_(level)
{
}

inline LogLine::LogLine(LogLine &&other) noexcept
    : session_(other.session_), level_(other.level_), buf_(std::move(other.buf_)), active_(other.active_)
{
    other.active_ = false;
}

inline LogLine::~LogLine()
{
    if (!active_ || session_ == nullptr)
        return;
    session_->write(level_, buf_.str());
}

inline LogLine &LogLine::operator<<(std::ostream &(*manip)(std::ostream &))
{
    buf_ << manip;
    return *this;
}

inline LogLine &LogLine::operator<<(std::ios_base &(*manip)(std::ios_base &))
{
    buf_ << manip;
    return *this;
}

inline LogLine LogChannel::operator<<(std::ostream &(*manip)(std::ostream &)) const
{
    LogLine line(session_, level_);
    line << manip;
    return line;
}

inline LogLine LogChannel::operator<<(std::ios_base &(*manip)(std::ios_base &)) const
{
    LogLine line(session_, level_);
    line << manip;
    return line;
}

inline bool SessionLog::open(const std::filesystem::path &directory)
{
    std::lock_guard<std::mutex> lock(mutex_);

    std::error_code ec;
    std::filesystem::create_directories(directory, ec);
    if (ec)
        return false;

    const auto now = std::chrono::system_clock::now();
    const std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    if (localtime_s(&tm, &tt) != 0)
        return false;
#else
    if (localtime_r(&tt, &tm) == nullptr)
        return false;
#endif

    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    std::ostringstream stamp;
    stamp << std::put_time(&tm, "%Y-%m-%dT%H-%M-%S")
          << '.' << std::setw(3) << std::setfill('0') << static_cast<int>(ms.count());

    if (file_.is_open())
        file_.close();
    file_.clear();

    path_ = directory / (stamp.str() + ".log");
    file_.open(path_, std::ios::out | std::ios::app);
    if (!file_)
    {
        path_.clear();
        return false;
    }

    std::ostringstream opened;
    opened << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S")
           << '.' << std::setw(3) << std::setfill('0') << static_cast<int>(ms.count());
    file_ << opened.str() << " Leonard\n";
    file_.flush();
    return true;
}

inline void SessionLog::error(const std::string &text)
{
    std::lock_guard<std::mutex> lock(mutex_);
    if (file_.is_open())
    {
        file_ << text;
        file_.flush();
    }
    std::cerr << text;
    std::cerr.flush();
}

inline void SessionLog::set_quiet(bool quiet)
{
    std::lock_guard<std::mutex> lock(mutex_);
    quiet_ = quiet;
}

inline void SessionLog::set_verbose(bool verbose)
{
    std::lock_guard<std::mutex> lock(mutex_);
    verbose_ = verbose;
}

inline bool SessionLog::is_open() const
{
    std::lock_guard<std::mutex> lock(mutex_);
    return file_.is_open();
}

inline std::filesystem::path SessionLog::path() const
{
    std::lock_guard<std::mutex> lock(mutex_);
    return path_;
}

inline void SessionLog::write(LogLevel level, const std::string &text)
{
    if (text.empty())
        return;

    std::lock_guard<std::mutex> lock(mutex_);
    if (file_.is_open())
    {
        file_ << text;
        file_.flush();
    }

    const bool echo = level == LogLevel::Say ? !quiet_ : verbose_;
    if (echo)
    {
        std::cout << text;
        std::cout.flush();
    }
}
