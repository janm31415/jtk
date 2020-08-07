#pragma once

#include <string>
#include <vector>
#include <iostream>

namespace jtk
  {
  void csv_read(std::vector<std::vector<std::string>>& data, FILE* stream, const char* separator = ",");
  bool csv_read(std::vector<std::vector<std::string>>& data, const char* filename, const char* separator = ",");

  void csv_write(const std::vector<std::vector<std::string>>& data, FILE* stream, const char* separator = ",");
  bool csv_write(const std::vector<std::vector<std::string>>& data, const char* filename, const char* separator = ",");

  namespace details
    {
    inline std::string clean_str(const std::string& str)
      {
      std::string cleaned;
      for (auto ch : str)
        {
        if (ch > 31 && ch < 127)
          cleaned.push_back(ch);
        }
      return cleaned;
      }

    inline std::string remove_brackets(const std::string& str)
      {
      if (str.empty())
        return str;
      if (str.front() == '"' && str.back() == '"')
        {
        return std::string(str.begin() + 1, str.end() - 1);
        }
      return str;
      }

    inline std::string add_brackets_iff_separator(const std::string& str, const char* separator = ",")
      {
      std::string w(str);
      if (w.find_first_of(separator) != std::string::npos)
        {
        w.insert(0, 1, '"');
        w.push_back('"');
        }
      return w;
      }

    inline void csv_print(const std::vector<std::vector<std::string>>& data, const char* separator = " | ")
      {
      for (const auto& line : data)
        {
        for (size_t i = 0; i < line.size() - 1; ++i)
          {
          std::cout << line[i] << separator;
          }
        std::cout << line.back() << "\n";
        }
      }

    inline const char* strpbrk_brackets(const char* str1, const char* str2)
      {
      const char* targ = strpbrk(str1, str2);
      if (!targ)
        return nullptr;
      const char* brackets1 = strpbrk(str1, "\"");
      if (brackets1)
        {
        if (brackets1 < targ)
          {
          const char* brackets2 = strpbrk(brackets1 + 1, "\"");
          if (!brackets2)
            return nullptr;
          return strpbrk_brackets(brackets2 + 1, str2);
          }
        else
          return targ;
        }
      else
        return targ;
      }
    }

  inline void csv_read(std::vector<std::vector<std::string>>& data, FILE* stream, const char* separator)
    {
    using namespace details;
    char line[16384];
    while (fgets(line, 16383, stream))
      {
      std::vector<std::string> dataline;
      const char* first = line;
      const char* last = strpbrk_brackets(line, separator);
      if (last != nullptr)
        dataline.push_back(remove_brackets(clean_str(std::string(first, last))));
      else
        dataline.push_back(remove_brackets(clean_str(std::string(first))));
      while (last != nullptr)
        {
        first = last + 1;
        last = strpbrk_brackets(last + 1, separator);
        if (last != nullptr)
          dataline.push_back(remove_brackets(clean_str(std::string(first, last))));
        else
          dataline.push_back(remove_brackets(clean_str(std::string(first))));
        }
      data.push_back(dataline);
      }
    }

  inline bool csv_read(std::vector<std::vector<std::string>>& data, const char* filename, const char* separator)
    {
    FILE *f = fopen(filename, "r");
    if (f == nullptr)
      return false;
    csv_read(data, f, separator);
    fclose(f);
    return true;
    }

  inline void csv_write(const std::vector<std::vector<std::string>>& data, FILE* stream, const char* separator)
    {
    using namespace details;
    for (const auto& line : data)
      {
      for (size_t i = 0; i < line.size() - 1; ++i)
        {
        std::string w = add_brackets_iff_separator(clean_str(line[i]), separator);
        fprintf(stream, "%s%s", w.c_str(), separator);
        }
      std::string w = add_brackets_iff_separator(clean_str(line.back()));
      fprintf(stream, "%s\n", w.c_str());
      }
    }

  bool csv_write(const std::vector<std::vector<std::string>>& data, const char* filename, const char* separator)
    {
    FILE *f = fopen(filename, "w");
    if (f == nullptr)
      return false;
    csv_write(data, f, separator);
    fclose(f);
    return true;
    }

  } // namespace jtk