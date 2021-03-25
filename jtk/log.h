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

   If you want to log functionality in a shared library/dll, and you want to set the 
   logging pointer from your main executable, then do the following:
    - make sure JTK_DLL_EXPORT is defined in your shared library so that the function 
      init_log_stream is exported.
    - define JTK_DLL_IMPORT in your main executable where you want to call init_log_stream.
    - call init_log_stream from your main executable with the requested ostream pointer.
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
  JTK_LOG_API JTKLDEF void init_log_stream(std::ostream* p_stream);
  JTK_LOG_API JTKLDEF void release_log_stream();

  JTKLDEF std::string now_time();

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


  class log
    {
    public:
      log();
      ~log();
      std::ostringstream& get(enum struct log_level level = log_level::info);

      static std::string to_string(enum struct log_level level);

    protected:
      std::ostringstream os;

    private:
      log(const log&);
      log& operator =(const log&);
    };


  } // namespace jtk


#ifdef JTK_LOG_IMPLEMENTATION

namespace jtk
  {
  std::ostream* log_stream_ptr = nullptr;

  JTKLDEF void init_log_stream(std::ostream* p_stream)
    {
    log_stream_ptr = p_stream;
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
    os << std::endl;     
    if (log_stream_ptr)
      {
      *log_stream_ptr << os.str();
      log_stream_ptr->flush();
      }
    }

  std::string log::to_string(enum struct log_level level)
    {
    static const char* const buffer[] = { "ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4" };
    switch (level)
      {
      case log_level::error: return buffer[0];
      case log_level::warning: return buffer[1];
      case log_level::info: return buffer[2];
      case log_level::debug: return buffer[3];
      case log_level::debug1: return buffer[4];
      case log_level::debug2: return buffer[5];
      case log_level::debug3: return buffer[6];
      case log_level::debug4: return buffer[7];
      }
    return std::string();
    }

  std::ostringstream& log::get(enum struct log_level level)
    {
    os << "- " << now_time();
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