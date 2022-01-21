/*
 Do this:
 #define JTK_HALFFLOAT_IMPLEMENTATION
 before you include this file in *one* C++ file to create the implementation.
 // i.e. it should look like this:
 #include ...
 #include ...
 #include ...
 #define JTK_HALFFLOAT_IMPLEMENTATION
 #include "jtk/halffloat.h"
 */

#ifndef JTK_HALFFLOAT_H
#define JTK_HALFFLOAT_H

#ifndef JTKHALFFLOATDEF
#ifdef JTK_HALFFLOAT_STATIC
#define JTKHALFFLOATDEF static
#else
#define JTKHALFFLOATDEF extern
#endif
#endif

#include <array>
#include <stdint.h>

namespace jtk
{

JTKHALFFLOATDEF float halffloat_to_float(uint16_t h);
JTKHALFFLOATDEF uint16_t float_to_halffloat(float f);

} // namespace jtk

#endif //JTK_HALFFLOAT_H

#ifdef JTK_HALFFLOAT_IMPLEMENTATION

namespace jtk
{

namespace halffloat_details
{

inline uint32_t convertmantissa(uint32_t i)
{
  auto m = i << 13; // Zero pad mantissa bits
  uint32_t e = 0;   // Zero exponent
  while (!(m & 0x00800000))
  {
    // While not normalized
    e -= 0x00800000; // Decrement exponent (1<<23)
    m <<= 1;         // Shift mantissa
  }
  m &= ~0x00800000; // Clear leading 1 bit
  e += 0x38800000;  // Adjust bias ((127-14)<<23)
  return m | e;     // Return combined number
}

std::array<uint32_t, 2048> generate_halffloat_mantissatable()
{
  std::array<uint32_t, 2048> out{};
  out[0] = 0;
  for (uint32_t i = 1; i < 1024; ++i)
    out[i] = convertmantissa(i);
  for (uint32_t i = 1024; i < 2048; ++i)
    out[i] = 0x38000000 + ((i - 1024) << 13);
  return out;
}

std::array<uint32_t, 64> generate_halffloat_exponenttable()
{
  std::array<uint32_t, 64> out{};
  out[0] = 0;
  out[32] = 0x80000000;
  for (uint32_t i = 1; i < 31; ++i)
    out[i] = i << 23;
  for (uint32_t i = 33; i < 63; ++i)
    out[i] = 0x80000000 + ((i - 32) << 23);
  out[31] = 0x47800000;
  out[63] = 0xC7800000;
  return out;
}

std::array<uint32_t, 64> generate_halffloat_offsettable()
{
  std::array<uint32_t, 64> out{};
  for (uint32_t i = 1; i < 64; ++i)
    out[i] = 1024;
  out[0] = 0;
  out[32] = 0;
  return out;
}

std::array<uint8_t, 512> generate_halffloat_shifttable()
{
  std::array<uint8_t, 512> shifttable;
  unsigned int i;
  int e;
  for (i = 0; i < 256; ++i) {
    e = i - 127;
    if (e < -24) { // Very small numbers map to zero
      shifttable[i | 0x000] = 24;
      shifttable[i | 0x100] = 24;
    }
    else if (e < -14) { // Small numbers map to denorms
      shifttable[i | 0x000] = -e - 1;
      shifttable[i | 0x100] = -e - 1;
    }
    else if (e <= 15) { // Normal numbers just lose precision
      shifttable[i | 0x000] = 13;
      shifttable[i | 0x100] = 13;
    }
    else if (e < 128) { // Large numbers map to Infinity
      shifttable[i | 0x000] = 24;
      shifttable[i | 0x100] = 24;
    }
    else { // Infinity and NaN's stay Infinity and NaN's
      shifttable[i | 0x000] = 13;
      shifttable[i | 0x100] = 13;
    }
  }
  return shifttable;
}

std::array<uint16_t, 512> generate_halffloat_basetable()
{
  std::array<uint16_t, 512> basetable;
  unsigned int i;
  int e;
  for (i = 0; i < 256; ++i) {
    e = i - 127;
    if (e < -24) { // Very small numbers map to zero
      basetable[i | 0x000] = 0x0000;
      basetable[i | 0x100] = 0x8000;
    }
    else if (e < -14) { // Small numbers map to denorms
      basetable[i | 0x000] = (0x0400 >> (-e - 14));
      basetable[i | 0x100] = (0x0400 >> (-e - 14)) | 0x8000;
    }
    else if (e <= 15) { // Normal numbers just lose precision
      basetable[i | 0x000] = ((e + 15) << 10);
      basetable[i | 0x100] = ((e + 15) << 10) | 0x8000;
    }
    else if (e < 128) { // Large numbers map to Infinity
      basetable[i | 0x000] = 0x7C00;
      basetable[i | 0x100] = 0xFC00;
    }
    else { // Infinity and NaN's stay Infinity and NaN's
      basetable[i | 0x000] = 0x7C00;
      basetable[i | 0x100] = 0xFC00;
    }
  }
  return basetable;
}

}

JTKHALFFLOATDEF float halffloat_to_float(uint16_t h)
{
  static const auto halffloat_mantissatable = halffloat_details::generate_halffloat_mantissatable();
  static const auto halffloat_exponenttable = halffloat_details::generate_halffloat_exponenttable();
  static const auto halffloat_offsettable = halffloat_details::generate_halffloat_offsettable();
  auto repr = halffloat_mantissatable[static_cast<size_t>(halffloat_offsettable[h >> 10]) + (h & 0x3ff)] + halffloat_exponenttable[h >> 10];
  return *reinterpret_cast<float*>(&repr);
}

JTKHALFFLOATDEF uint16_t float_to_halffloat(float f)
{
  static const auto halffloat_basetable = halffloat_details::generate_halffloat_basetable();
  static const auto halffloat_shifttable = halffloat_details::generate_halffloat_shifttable();
  uint32_t f32 = *reinterpret_cast<uint32_t*>(&f);
  uint16_t h = halffloat_basetable[(f32 >> 23) & 0x1ff] + (uint16_t)((f32 & 0x007fffff) >> halffloat_shifttable[(f32 >> 23) & 0x1ff]);
  return h;
  
}


} // namespace jtk
#endif // JTK_HALFFLOAT_IMPLEMENTATION
