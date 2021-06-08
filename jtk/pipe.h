/*
   Do this:
      #define JTK_PIPE_IMPLEMENTATION
   before you include this file in *one* C++ file to create the implementation.
   // i.e. it should look like this:
   #include ...
   #include ...
   #include ...
   #define JTK_PIPE_IMPLEMENTATION
   #include "jtk/pipe.h"
 */

#ifndef JTK_PIPE_H
#define JTK_PIPE_H

#ifndef JTKPIPEDEF
#ifdef JTK_PIPE_STATIC
#define JTKPIPEDEF static
#define JTK_FILE_UTILS_STATIC
#else
#define JTKPIPEDEF extern
#endif
#endif

#ifndef JTKPIPEINLINE
#define JTKPIPEINLINE inline
#endif

#ifdef JTK_FILE_UTILS_IMPLEMENTATION
  #undef JTK_FILE_UTILS_IMPLEMENTATION
  #include "file_utils.h"
  #define JTK_FILE_UTILS_IMPLEMENTATION
#else
  #include "file_utils.h"
#endif


#if defined(_WIN32)
#include <windows.h>
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(UNIX)
#include <linux/limits.h>
#include <unistd.h>
#elif defined(__APPLE__)
#include <limits.h>
#include <unistd.h>
#endif

namespace jtk
  {

  class active_folder;
#ifdef _WIN32
  JTKPIPEDEF int create_pipe(const char *path, char * const * argv, const char* current_dir, void** pr);
  JTKPIPEDEF void close_pipe(void* pr);
  JTKPIPEDEF void destroy_pipe(void* pr, int signal);
  JTKPIPEDEF int send_to_pipe(void* process, const char* message);
  JTKPIPEDEF std::string read_from_pipe(void* process, int time_out);
  JTKPIPEDEF std::string read_std_input(int time_out);
  JTKPIPEDEF int run_process(const char *path, char * const * argv, const char* current_dir, void** pr);
  JTKPIPEDEF void destroy_process(void* pr, int signal);
#else
  JTKPIPEDEF int create_pipe(const char *path, char* const* argv, const char* current_dir, int* pipefd);
  JTKPIPEDEF void close_pipe(int* pipefd);
  JTKPIPEDEF void destroy_pipe(int* pipefd, int);
  JTKPIPEDEF int send_to_pipe(int* pipefd, const char* message);
  JTKPIPEDEF std::string read_from_pipe(int* pipefd, int time_out);
  JTKPIPEDEF std::string read_std_input(int time_out);
  JTKPIPEDEF int run_process(const char *path, char * const * argv, const char* current_dir, pid_t* process);
  JTKPIPEDEF void destroy_process(pid_t process, int signal);
#endif

#if defined(_WIN32)
  #define JTK_MAX_PATH MAX_PATH
#else
  #define JTK_MAX_PATH PATH_MAX
#endif

  class active_folder
    {
    public:
      active_folder(const char* folder);
      ~active_folder();

    private:
    #ifdef _WIN32
      wchar_t buf[JTK_MAX_PATH];
    #else
      char buf[JTK_MAX_PATH];
    #endif
    };
    
#ifdef _WIN32

  JTKPIPEINLINE active_folder::active_folder(const char* folder)
    {
    GetCurrentDirectoryW(JTK_MAX_PATH, buf);
    if (folder)
      {
      std::wstring wdir(convert_string_to_wstring(std::string(folder)));
      SetCurrentDirectoryW(wdir.c_str());
      }
    }

  JTKPIPEINLINE active_folder::~active_folder()
    {
    SetCurrentDirectoryW(buf);
    }
#else

  JTKPIPEINLINE active_folder::active_folder(const char* folder)
    {
    ::getcwd(buf, sizeof(buf));
    if (folder)
      chdir(folder);
    }

  JTKPIPEINLINE active_folder::~active_folder()
    {
    ::chdir(buf);
    }
    
#endif

  } // namespace jtk
  
#endif // #ifndef JTK_PIPE_H


#ifdef JTK_PIPE_IMPLEMENTATION

#ifdef _WIN32
//#include <windows.h>
#include <string>
#include <winsock.h>
#else
#include <sstream>
#include <iostream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#if defined(unix) || defined(__unix) || defined(__unix__) || defined(UNIX)
#include <sys/prctl.h>
//#include <linux/limits.h>
#endif
#include <signal.h>
#include <fcntl.h>
#endif

#include <chrono>
#include <thread>

namespace jtk
  {
    
#define MAX_PIPE_BUFFER_SIZE 4096

#ifdef _WIN32

  struct pipe_process
    {
    HANDLE hProcess;
    DWORD pid;
    HANDLE hTo;
    HANDLE hFrom;
    };

  JTKPIPEDEF int create_pipe(const char *path, char * const * argv, const char* current_dir, void** pr)
    {
    HANDLE hChildStdinRd, hChildStdinWr, hChildStdoutRd, hChildStdoutWr;
    HANDLE hChildStdinWrDup, hChildStdoutRdDup;
    SECURITY_ATTRIBUTES saAttr;
    BOOL fSuccess;
    PROCESS_INFORMATION piProcInfo;
    STARTUPINFOW siStartInfo;
    pipe_process *cp;

    DWORD err;

    *pr = nullptr;

    /* Set the bInheritHandle flag so pipe handles are inherited. */
    saAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
    saAttr.bInheritHandle = TRUE;
    saAttr.lpSecurityDescriptor = NULL;

    /*
    * The steps for redirecting child's STDOUT:
    *     1. Create anonymous pipe to be STDOUT for child.
    *     2. Create a noninheritable duplicate of read handle,
    *         and close the inheritable read handle.
    */

    /* Create a pipe for the child's STDOUT. */
    if (!CreatePipe(&hChildStdoutRd, &hChildStdoutWr, &saAttr, 0))
      {
      return GetLastError();
      }

    /* Duplicate the read handle to the pipe, so it is not inherited. */
    fSuccess = DuplicateHandle(GetCurrentProcess(), hChildStdoutRd,
      GetCurrentProcess(), &hChildStdoutRdDup, 0,
      FALSE,	/* not inherited */
      DUPLICATE_SAME_ACCESS);
    if (!fSuccess)
      {
      return GetLastError();
      }
    CloseHandle(hChildStdoutRd);

    /*
    * The steps for redirecting child's STDIN:
    *     1. Create anonymous pipe to be STDIN for child.
    *     2. Create a noninheritable duplicate of write handle,
    *         and close the inheritable write handle.
    */

    /* Create a pipe for the child's STDIN. */
    if (!CreatePipe(&hChildStdinRd, &hChildStdinWr, &saAttr, 0)) {
      return GetLastError();
      }

    /* Duplicate the write handle to the pipe, so it is not inherited. */
    fSuccess = DuplicateHandle(GetCurrentProcess(), hChildStdinWr,
      GetCurrentProcess(), &hChildStdinWrDup, 0,
      FALSE,	/* not inherited */
      DUPLICATE_SAME_ACCESS);
    if (!fSuccess) {
      return GetLastError();
      }
    CloseHandle(hChildStdinWr);

    /* Arrange to (1) look in dir for the child .exe file, and
    * (2) have dir be the child's working directory.  Interpret
    * dir relative to the directory WinBoard loaded from. */

    active_folder af(current_dir);

    std::wstring wcmdLine;
    wcmdLine = convert_string_to_wstring(std::string(path));
    std::replace(wcmdLine.begin(), wcmdLine.end(), L'/', L'\\');
    // i = 0 equals path, is for linux
    size_t i = 1;
    const char* arg = argv[i];
    while (arg)
      {
      std::wstring warg = convert_string_to_wstring(std::string(arg));
      std::replace(warg.begin(), warg.end(), L'/', L'\\');
      wcmdLine.append(L" " + warg);
      arg = argv[++i];
      }

    /* Now create the child process. */

    siStartInfo.cb = sizeof(STARTUPINFOW);
    siStartInfo.lpReserved = NULL;
    siStartInfo.lpDesktop = NULL;
    siStartInfo.lpTitle = NULL;
    siStartInfo.dwFlags = STARTF_USESTDHANDLES;
    siStartInfo.cbReserved2 = 0;
    siStartInfo.lpReserved2 = NULL;
    siStartInfo.hStdInput = hChildStdinRd;
    siStartInfo.hStdOutput = hChildStdoutWr;
    siStartInfo.hStdError = hChildStdoutWr;

    fSuccess = CreateProcessW(NULL,
      (LPWSTR)(wcmdLine.c_str()),	   /* command line */
      NULL,	   /* process security attributes */
      NULL,	   /* primary thread security attrs */
      TRUE,	   /* handles are inherited */
      DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
      NULL,	   /* use parent's environment */
      NULL,
      &siStartInfo, /* STARTUPINFO pointer */
      &piProcInfo); /* receives PROCESS_INFORMATION */

    err = GetLastError();

    if (!fSuccess)
      {
      return err;
      }

    SetPriorityClass(piProcInfo.hProcess, 0x00000080);

    /* Close the handles we don't need in the parent */
    CloseHandle(piProcInfo.hThread);
    CloseHandle(hChildStdinRd);
    CloseHandle(hChildStdoutWr);

    /* Prepare return value */
    cp = (pipe_process *)calloc(1, sizeof(pipe_process));
    cp->hProcess = piProcInfo.hProcess;
    cp->pid = piProcInfo.dwProcessId;
    cp->hFrom = hChildStdoutRdDup;
    cp->hTo = hChildStdinWrDup;

    *pr = (void *)cp;

    return NO_ERROR;
    }

  JTKPIPEDEF void close_pipe(void* pr)
    {
    pipe_process *cp;

    cp = (pipe_process *)pr;
    if (cp == nullptr)
      return;


    /* TerminateProcess is considered harmful, so... */
    CloseHandle(cp->hTo); /* Closing this will give the child an EOF and hopefully kill it */
    if (cp->hFrom) CloseHandle(cp->hFrom);  /* if NULL, InputThread will close it */
                                            /* The following doesn't work because the chess program
                                            doesn't "have the same console" as WinBoard.  Maybe
                                            we could arrange for this even though neither WinBoard
                                            nor the chess program uses a console for stdio? */
                                            /*!!if (signal) GenerateConsoleCtrlEvent(CTRL_BREAK_EVENT, cp->pid);*/

                                            /* [AS] Special termination modes for misbehaving programs... */


    CloseHandle(cp->hProcess);

    free(cp);
    }

  JTKPIPEDEF void destroy_pipe(void* pr, int signal)
    {
    pipe_process *cp;
    int result;

    cp = (pipe_process *)pr;
    if (cp == nullptr)
      return;


    /* TerminateProcess is considered harmful, so... */
    CloseHandle(cp->hTo); /* Closing this will give the child an EOF and hopefully kill it */
    if (cp->hFrom) CloseHandle(cp->hFrom);  /* if NULL, InputThread will close it */
                                            /* The following doesn't work because the chess program
                                            doesn't "have the same console" as WinBoard.  Maybe
                                            we could arrange for this even though neither WinBoard
                                            nor the chess program uses a console for stdio? */
                                            /*!!if (signal) GenerateConsoleCtrlEvent(CTRL_BREAK_EVENT, cp->pid);*/

                                            /* [AS] Special termination modes for misbehaving programs... */
    if (signal == 9)
      {
      result = TerminateProcess(cp->hProcess, 0);
      }
    else if (signal == 10)
      {
      DWORD dw = WaitForSingleObject(cp->hProcess, 3 * 1000); // Wait 3 seconds at most

      if (dw != WAIT_OBJECT_0)
        {
        result = TerminateProcess(cp->hProcess, 0);
        }
      }

    CloseHandle(cp->hProcess);

    free(cp);
    }

  JTKPIPEDEF int send_to_pipe(void* process, const char* message)
    {
    int count;
    DWORD dOutCount;
    if (process == nullptr)
      return ERROR_INVALID_HANDLE;
    count = (int)strlen(message);
    if (WriteFile(((pipe_process *)process)->hTo, message, count, &dOutCount, NULL))
      {
      if (count == (int)dOutCount)
        return NO_ERROR;
      else
        return (int)GetLastError();
      }
    return SOCKET_ERROR;
    }

  JTKPIPEDEF std::string read_from_pipe(void* process, int time_out)
    {
    std::string input;

    if (process == nullptr)
      return "";
    pipe_process *cp = (pipe_process *)process;

    DWORD bytes_left = 0;

    auto tic = std::chrono::steady_clock::now();
    auto toc = std::chrono::steady_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count();

    bool check_at_least_once = true;

    while (time_elapsed < time_out || check_at_least_once)
      {
      check_at_least_once = false;
      if (PeekNamedPipe(cp->hFrom, NULL, 0, NULL, &bytes_left, NULL) && bytes_left)
        {
        check_at_least_once = true;
        bool bytes_left_to_read = bytes_left > 0;

        char line[MAX_PIPE_BUFFER_SIZE];

        while (bytes_left_to_read)
          {
          int n = bytes_left;
          if (n >= MAX_PIPE_BUFFER_SIZE)
            n = MAX_PIPE_BUFFER_SIZE - 1;
          memset(line, 0, MAX_PIPE_BUFFER_SIZE);
          DWORD count;
          if (!ReadFile(cp->hFrom, line, n, &count, nullptr))
            return input;
          std::string str(line);
          input.append(str);

          bytes_left -= (DWORD)str.size();
          bytes_left_to_read = bytes_left > 0;
          }
        }
      std::this_thread::sleep_for(std::chrono::milliseconds(time_out / 2 > 10 ? 10 : time_out / 2));
      toc = std::chrono::steady_clock::now();
      time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count();
      }
    return input;
    }

  JTKPIPEDEF std::string read_std_input(int time_out)
    {
    pipe_process pr;
    pr.hFrom = GetStdHandle(STD_INPUT_HANDLE);
    return read_from_pipe(&pr, time_out);
    }


  struct process_info
    {
    HANDLE hProcess;
    DWORD pid;
    };

  JTKPIPEDEF int run_process(const char *path, char * const * argv, const char* current_dir, void** pr)
    {
    STARTUPINFOW siStartInfo;
    process_info *cp;
    DWORD err;
    BOOL fSuccess;
    PROCESS_INFORMATION piProcInfo;

    *pr = nullptr;

    active_folder af(current_dir);

    std::wstring wcmdLine;
    wcmdLine = convert_string_to_wstring(std::string(path));
    std::replace(wcmdLine.begin(), wcmdLine.end(), L'/', L'\\');

    // i = 0 equals path, is for linux
    size_t i = 1;
    const char* arg = argv[i];
    while (arg)
      {
      std::wstring warg = convert_string_to_wstring(std::string(arg));
      std::replace(warg.begin(), warg.end(), L'/', L'\\');
      wcmdLine.append(L" " + warg);
      arg = argv[++i];
      }
    /* Now create the child process. */

    siStartInfo.cb = sizeof(STARTUPINFOW);
    siStartInfo.lpReserved = NULL;
    siStartInfo.lpDesktop = NULL;
    siStartInfo.lpTitle = NULL;
    siStartInfo.dwFlags = 0;
    siStartInfo.cbReserved2 = 0;
    siStartInfo.lpReserved2 = NULL;
    siStartInfo.hStdInput = NULL;
    siStartInfo.hStdOutput = NULL;
    siStartInfo.hStdError = NULL;

    fSuccess = CreateProcessW(NULL,
      (LPWSTR)(wcmdLine.c_str()),	   /* command line */
      NULL,	   /* process security attributes */
      NULL,	   /* primary thread security attrs */
      TRUE,	   /* handles are inherited */
      CREATE_NEW_CONSOLE,
      NULL,	   /* use parent's environment */
      NULL,
      &siStartInfo, /* STARTUPINFO pointer */
      &piProcInfo); /* receives PROCESS_INFORMATION */

    err = GetLastError();

    if (!fSuccess)
      {
      return err;
      }

    SetPriorityClass(piProcInfo.hProcess, 0x00000080);

    /* Close the handles we don't need in the parent */
    CloseHandle(piProcInfo.hThread);

    /* Prepare return value */
    cp = (process_info *)calloc(1, sizeof(process_info));
    cp->hProcess = piProcInfo.hProcess;
    cp->pid = piProcInfo.dwProcessId;

    *pr = (void *)cp;

    return NO_ERROR;
    }

  JTKPIPEDEF void destroy_process(void* pr, int signal)
    {
    process_info *cp;
    int result;

    cp = (process_info *)pr;
    if (cp == nullptr)
      return;

    if (signal == 9)
      {
      result = TerminateProcess(cp->hProcess, 0);
      }
    else if (signal == 10)
      {
      DWORD dw = WaitForSingleObject(cp->hProcess, 3 * 1000); // Wait 3 seconds at most

      if (dw != WAIT_OBJECT_0)
        {
        result = TerminateProcess(cp->hProcess, 0);
        }
      }

    CloseHandle(cp->hProcess);

    free(cp);
    }

#else

  JTKPIPEDEF int create_pipe(const char *path, char* const* argv, const char* current_dir, int* pipefd)
    {
    pid_t pid = 0;
    int inpipefd[2];
    int outpipefd[2];
    if (pipe(inpipefd) != 0)
      {
      throw std::runtime_error("failed to pipe");
      }
    if (pipe(outpipefd) != 0)
      {
      throw std::runtime_error("failed to pipe");
      }
    pid = fork();
    if (pid < 0)
      throw std::runtime_error("failed to fork");
    if (pid == 0)
      {
      dup2(outpipefd[0], STDIN_FILENO);
      dup2(inpipefd[1], STDOUT_FILENO);

      close(inpipefd[0]);
      close(inpipefd[1]);
      close(outpipefd[0]);
      close(outpipefd[1]);
#if defined(unix) || defined(__unix) || defined(__unix__) || defined(UNIX)
      prctl(PR_SET_PDEATHSIG, SIGTERM);
#endif

      execv(path, argv);
      //if (execv(path, argv) == -1)
      //  throw std::runtime_error("failed to pipe (execl failed)"); // remove throwing because it causes error messages on the os when a pipe fails
      exit(1);
      }

    close(outpipefd[0]);
    close(inpipefd[1]);
    pipefd[0] = outpipefd[1];
    pipefd[1] = inpipefd[0];

    auto flags = fcntl(pipefd[0], F_GETFL, 0);
    flags |= O_NONBLOCK;
    fcntl(pipefd[0], F_SETFL, flags);

    flags = fcntl(pipefd[1], F_GETFL, 0);
    flags |= O_NONBLOCK;
    fcntl(pipefd[1], F_SETFL, flags);

    pipefd[2] = pid;

    return 0;
    }

  JTKPIPEDEF void close_pipe(int* pipefd)
    {
    if (pipefd[2] >= 0)
      {
      close(pipefd[0]);
      close(pipefd[1]);
      }
    }

  JTKPIPEDEF void destroy_pipe(int* pipefd, int)
    {
    if (pipefd[2] >= 0)
      {
      kill(pipefd[2], SIGKILL);
      int status;
      waitpid(pipefd[2], &status, 0);
      close(pipefd[0]);
      close(pipefd[1]);
      }
    }
    
  JTKPIPEDEF void handler(int s) {
    std::cerr << "Caught SIGPIPE with signal " << s << std::endl;
    }

  JTKPIPEDEF int send_to_pipe(int* pipefd, const char* message)
    {
    signal(SIGPIPE, handler);
    auto v = write(pipefd[0], message, strlen(message));
    return v == strlen(message) ? 0 : 1;
    }

  JTKPIPEDEF std::string read_from_pipe(int* pipefd, int time_out)
    {
    std::stringstream ss;
    auto tic = std::chrono::steady_clock::now();
    auto toc = std::chrono::steady_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count();

    bool check_at_least_once = true;
    char buffer[MAX_PIPE_BUFFER_SIZE];
    while (time_elapsed < time_out || check_at_least_once)
      {
      check_at_least_once = false;
      memset(buffer, 0, MAX_PIPE_BUFFER_SIZE);
      int num_read = read(pipefd[1], buffer, MAX_PIPE_BUFFER_SIZE - 1);
      if (num_read > 0)
        {
        ss << buffer;
        check_at_least_once = true;
        }
      toc = std::chrono::steady_clock::now();
      time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count();
      std::this_thread::sleep_for(std::chrono::milliseconds(time_out / 2 > 10 ? 10 : time_out / 2));
      }
    return ss.str();
    }



  JTKPIPEDEF std::string read_std_input(int time_out)
    {

    std::stringstream ss;


    auto flags = fcntl(STDIN_FILENO, F_GETFL, 0);
    flags |= O_NONBLOCK;
    fcntl(STDIN_FILENO, F_SETFL, flags);


    auto tic = std::chrono::steady_clock::now();
    auto toc = std::chrono::steady_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count();

    bool check_at_least_once = true;
    char buffer[MAX_PIPE_BUFFER_SIZE];
    while (time_elapsed < time_out || check_at_least_once)
      {
      check_at_least_once = false;
      memset(buffer, 0, MAX_PIPE_BUFFER_SIZE);
      int num_read = read(STDIN_FILENO, buffer, MAX_PIPE_BUFFER_SIZE - 1);
      if (num_read > 0)
        {
        ss << buffer;
        check_at_least_once = true;
        }
      toc = std::chrono::steady_clock::now();
      time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count();
      std::this_thread::sleep_for(std::chrono::milliseconds(time_out / 2 > 10 ? 10 : time_out / 2));
      }
    return ss.str();

    }

  JTKPIPEDEF int run_process(const char *path, char * const * argv, const char* current_dir, pid_t* process)
    {
    pid_t pid = fork();

    if (pid < 0)
      {
      printf("\nfailed to fork child\n");
      return -1;
      }
    else if (pid == 0)
      {
      if (execv(path, argv) < 0)
        printf("\nCould not run process (execv failed)\n");
      exit(0);
      }
    *process = pid;
    return 0;
    }

  void destroy_process(pid_t process, int signal)
    {
    if (signal == 9)
      kill(process, SIGKILL);

    if (signal == 10)
      {
      int status;
      waitpid(process, &status, 0);
      }
    }

#endif

  } // namespace jtk

#endif // JTK_PIPE_IMPLEMENTATION
