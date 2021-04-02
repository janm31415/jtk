/*
   Do this:
      #define JTK_LOG_IMPLEMENTATION
   before you include this file in *one* C++ file to create the implementation.
   // i.e. it should look like this:
   #include ...
   #include ...
   #include ...
   #define JTK_LOG_IMPLEMENTATION
   #include "jtk/log.h"
 */

#ifndef JTK_LOG_H
#define JTK_LOG_H

#include <sstream>
#include <string>
#include <stdio.h>
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>

#ifndef JTK_LOG_API
#if defined(JTK_DLL_EXPORT)
#define JTK_LOG_API __declspec (dllexport)
#elif defined(JTK_DLL_IMPORT)
#define JTK_LOG_API __declspec (dllimport)
#else
#define JTK_LOG_API
#endif
#endif

#else
#include <sys/time.h>
#ifndef JTK_LOG_API
#define JTK_LOG_API
#endif
#endif

#ifndef JTKLDEF
#define JTKLDEF extern
#endif

namespace jtk
  {

  enum struct log_level
    {
    error,
    warning,
    info,
    debug,
    debug1,
    debug2,
    debug3,
    debug4
    };

  JTK_LOG_API JTKLDEF void init_log_stream(std::ostream* p_stream, bool print_time = true, log_level verbose_level = log_level::info);
  JTK_LOG_API JTKLDEF void release_log_stream();

  JTKLDEF std::string now_time();


  class log
    {
    public:
      log();
      ~log();
      std::ostringstream& get(log_level level = log_level::info);

      static std::string to_string(log_level level);
      static int to_int(log_level level);

    protected:
      std::ostringstream os;
      std::ostringstream dummy_os;

    private:
      log(const log&);
      log& operator =(const log&);
    };


  } // namespace jtk


#ifdef JTK_LOG_IMPLEMENTATION

namespace jtk
  {
  std::ostream* log_stream_ptr = nullptr;
  log_level verbose_level = log_level::debug4;
  bool log_with_time_stamp = true;

  JTKLDEF void init_log_stream(std::ostream* p_stream, bool with_time_stamp, log_level vlevel)
    {
    log_stream_ptr = p_stream;
    verbose_level = vlevel;
    log_with_time_stamp = with_time_stamp;
    }

  JTKLDEF void release_log_stream()
    {
    log_stream_ptr = nullptr;
    }

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
  JTKLDEF std::string now_time()
    {
    const int MAX_LEN = 200;
    char buffer[MAX_LEN];
    if (GetTimeFormatA(LOCALE_USER_DEFAULT, 0, 0,
      "HH':'mm':'ss", buffer, MAX_LEN) == 0)
      return "Error in now_time()";

    char result[100] = { 0 };
    static DWORD first = GetTickCount();
    std::sprintf(result, "%s.%03ld", buffer, (long)(GetTickCount() - first) % 1000);
    return result;
    }
#else
  JTKLDEF std::string now_time()
    {
    char buffer[11];
    time_t t;
    time(&t);
    tm r = { 0 };
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&t, &r));
    struct timeval tv;
    gettimeofday(&tv, 0);
    char result[100] = { 0 };
    std::sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000);
    return result;
    }
#endif

  log::log()
    {
    }

  log::~log()
    {
    if (!os.str().empty())
      {
      os << std::endl;
      if (log_stream_ptr)
        {
        *log_stream_ptr << os.str();
        log_stream_ptr->flush();
        }
      }
    }

  int log::to_int(log_level level)
    {
    switch (level)
      {
      case log_level::error: return 0;
      case log_level::warning: return 1;
      case log_level::info: return 2;
      case log_level::debug: return 3;
      case log_level::debug1: return 4;
      case log_level::debug2: return 5;
      case log_level::debug3: return 6;
      case log_level::debug4: return 7;
      }
    return -1;
    }

  std::string log::to_string(log_level level)
    {
    static const char* const buffer[] = { "ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4" };
    return buffer[to_int(level)];
    }

  std::ostringstream& log::get(log_level level)
    {
    if (to_int(level) > to_int(verbose_level))
      return dummy_os;
    if (log_with_time_stamp)
      {
      os << "- " << now_time();
      }
    os << " " << to_string(level) << ": ";
    switch (level)
      {
      default: break;
      case log_level::debug1: os << std::string(1, '\t'); break;
      case log_level::debug2: os << std::string(2, '\t'); break;
      case log_level::debug3: os << std::string(3, '\t'); break;
      case log_level::debug4: os << std::string(4, '\t'); break;
      }
    return os;
    }

  } // namespace jtk

#endif

#endif

