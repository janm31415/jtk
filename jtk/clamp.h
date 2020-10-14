#pragma once

#undef min
#undef max

#include <limits.h>

namespace jtk
  {

  namespace clamp_details
    {
    template <class T2>
    struct _Clamp
      {
      template <class T1>
      T1 clmp(T2 val)
        {
        if (val > T2(std::numeric_limits<T1>::max()))
          return std::numeric_limits<T1>::max();
        if (val < T2(std::numeric_limits<T1>::min()))
          return std::numeric_limits<T1>::min();
        return T1(val);
        }

      template <>
      T2 clmp<T2>(T2 val)
        {
        return val;
        }

      };
    }

  template <class T1, class T2>
  T1 clamp(T2 val)
    {
    return clamp_details::_Clamp<T2>().template clmp<T1>(val);
    }

  }
