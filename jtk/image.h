#pragma once

#ifdef _JTK_FOR_ARM
#include "sse2neon.h"
#else
#include <immintrin.h>
#include <emmintrin.h>
#endif

#include <array>
#include <cassert>
#include <fstream>
#include <stdint.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cmath>
#include <vector>

namespace jtk
  {

  template <class T>
  class image;

  inline bool load_pgm(image<uint16_t>& im, const std::string& filename);

  template <class T>
  inline image<T> span_to_image(uint32_t w, uint32_t h, uint32_t stride, const T* p_image);

  template <class T>
  void draw_line(T* raw, double x1, double y1, double x2, double y2, const uint32_t w, const uint32_t h, const uint32_t stride, const T& color);

  template <class T>
  void draw_box(T* raw, int x, int y, int box_width, int box_height, const uint32_t w, const uint32_t h, const uint32_t stride, const T& color);

  template <class T>
  image<T> rotate90(const image<T>& im);

  template <class T>
  image<T> rotate180(const image<T>& im);

  template <class T>
  image<T> rotate270(const image<T>& im);

  template <class T>
  void fill_image(image<T>& im, T value);

  image<uint8_t> image_to_gray(const image<uint32_t>& im);
  image<uint8_t> image_to_gray(const image<float>& im);

  void minmax(float& min, float& max, const image<float>& im, bool skip_infinity);

  template <class T, class TValid>
  void convolve_col_14641_div_16(image<T>& out, const image<T>& im, TValid valid);

  template <class T, class TValid>
  void convolve_row_14641_div_16(image<T>& out, const image<T>& im, TValid valid);

  template <class T, class TValid>
  void gauss(image<T>& im, TValid valid);

  typedef std::vector<uint32_t> histogram_t;

  histogram_t make_histogram(const image<uint8_t>& im);

  histogram_t make_histogram(const image<uint8_t>& im, int x, int y, int w, int h);

  int64_t get_histogram_weight(histogram_t::const_iterator begin, histogram_t::const_iterator end);

  int otsu_threshold(const histogram_t& histogram);

  void threshold_image_region(image<uint8_t>& im, int x, int y, int w, int h, uint8_t threshold);

  image<uint8_t> make_binary_image(const image<uint8_t>& im, int steps_w, int steps_h);

  image<uint8_t> make_binary_image(const image<uint8_t>& im, uint8_t threshold);

  image<uint32_t> three_gray_to_uint32_t(const image<uint8_t>& r, const image<uint8_t>& g, const image<uint8_t>& b);

  image<uint32_t> census_transform(const image<uint8_t>& im);

  void census_transform(image<uint32_t>& census_red, image<uint32_t>& census_green, image<uint32_t>& census_blue, const image<uint32_t>& im);

  template <class T>
  class image
    {
    public:

      typedef T* iterator;
      typedef const T* const_iterator;
      typedef T value_type;

      image() : _w(0), _h(0), _stride(0), _data(nullptr), _access(nullptr) {}

      image(uint32_t w, uint32_t h, bool init_to_zero = true) : _w(w), _h(h), _stride(w), _data(nullptr), _access(nullptr)
        {
        if ((_w * sizeof(T)) & 15)
          {
          switch (sizeof(T))
            {
            case 1:
            {
            _stride += 16 - (_stride & 15);
            break;
            }
            case 2:
            {
            _stride += 8 - (_stride & 7);
            break;
            }
            case 4:
            {
            _stride += 4 - (_stride & 3);
            break;
            }
            default:
            {
            ++_stride;
            while ((_stride * sizeof(T)) & 15)
              ++_stride;
            break;
            }
            }
          }
        _data = (T*)(_mm_malloc(_stride*_h * sizeof(T), 16));
        _access = (T**)(malloc(_h * sizeof(T*)));
        for (uint32_t i = 0; i < h; ++i)
          _access[i] = _data + (i * _stride);
        if (init_to_zero)
          memset(_data, 0, _stride * _h * sizeof(T));
        }

      image(image&& other) : _w(0), _h(0), _stride(0), _data(nullptr), _access(nullptr)
        {
        _data = other._data;
        _access = other._access;
        _w = other._w;
        _h = other._h;
        _stride = other._stride;
        other._data = nullptr;
        other._access = nullptr;
        other._w = 0;
        other._h = 0;
        other._stride = 0;
        }

      image& operator = (image&& other)
        {
        std::swap(_data, other._data);
        std::swap(_access, other._access);
        std::swap(_w, other._w);
        std::swap(_h, other._h);
        std::swap(_stride, other._stride);
        return *this;
        }

      image(const image& other) : _w(other._w), _h(other._h), _stride(other._stride), _data(nullptr), _access(nullptr)
        {
        _data = (T*)(_mm_malloc(_stride*_h * sizeof(T), 16));
        _access = (T**)(malloc(_h * sizeof(T*)));
        for (uint32_t i = 0; i < _h; ++i)
          _access[i] = _data + (i * _stride);
        memcpy(_data, other._data, _stride * _h * sizeof(T));
        }

      void swap(image& other)
        {
        std::swap(_w, other._w);
        std::swap(_h, other._h);
        std::swap(_stride, other._stride);
        std::swap(_data, other._data);
        std::swap(_access, other._access);
        }

      image& operator = (const image& other)
        {
        image temp(other);
        swap(temp);
        return *this;
        }

      ~image()
        {
        _mm_free(_data);
        free(_access);
        }

      const uint32_t get_index(uint32_t x, uint32_t y) const
        {
        return y * _stride + x;
        }

      const uint32_t width() const
        {
        return _w;
        }

      const uint32_t height() const
        {
        return _h;
        }

      const uint32_t stride() const
        {
        return _stride;
        }

      T& operator [] (uint32_t ind)
        {
        return _data[ind];
        }

      const T& operator [] (uint32_t ind) const
        {
        return _data[ind];
        }

      T& operator () (uint32_t x, uint32_t y)
        {
        return _access[y][x];
        }

      const T& operator () (uint32_t x, uint32_t y) const
        {
        return _access[y][x];
        }

      iterator begin()
        {
        return _data;
        }

      iterator end()
        {
        return _data + _stride * _h;
        }

      const_iterator begin() const
        {
        return _data;
        }

      const_iterator end() const
        {
        return _data + _stride * _h;
        }

      T* data()
        {
        return _data;
        }

      const T* data() const
        {
        return _data;
        }

      T* row(uint32_t y)
        {
        return _access[y];
        }

      const T* row(uint32_t y) const
        {
        return _access[y];
        }

    private:
      uint32_t _w, _h, _stride;
      T* _data;
      T** _access;
    };

  namespace details
    {
    const int census_samples[] = {
      -3, -2,
      -3, 0,
      -3, 2,
      -2, -3,
      -2, -1,
      -2, 1,
      -2, 3,
      -1, -2,
      -1, 0,
      -1, 2,
      0, -3,
      0, -1,
      0, 1,
      0, 3,
      1, -2,
      1, 0,
      1, 2,
      2, -3,
      2, -1,
      2, 1,
      2, 3,
      3, -2,
      3, 0,
      3, 2
      };

    inline float horizontal_min_ps(__m128 x)
      {
      __m128 max1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
      __m128 max2 = _mm_min_ps(x, max1);
      __m128 max3 = _mm_shuffle_ps(max2, max2, _MM_SHUFFLE(0, 0, 0, 1));
      __m128 max4 = _mm_min_ps(max2, max3);
      float result = _mm_cvtss_f32(max4);
      return result;
      }

    inline float horizontal_max_ps(__m128 x)
      {
      __m128 max1 = _mm_shuffle_ps(x, x, _MM_SHUFFLE(0, 0, 3, 2));
      __m128 max2 = _mm_max_ps(x, max1);
      __m128 max3 = _mm_shuffle_ps(max2, max2, _MM_SHUFFLE(0, 0, 0, 1));
      __m128 max4 = _mm_max_ps(max2, max3);
      float result = _mm_cvtss_f32(max4);
      return result;
      }

    inline __m128 set_data_128(float value)
      {
      return _mm_set1_ps(value);
      }

    inline __m128 set_data_128(uint8_t value)
      {
      return _mm_castsi128_ps(_mm_set1_epi8(value));
      }

    inline __m128 set_data_128(uint16_t value)
      {
      return _mm_castsi128_ps(_mm_set1_epi16(value));
      }

    inline __m128 set_data_128(uint32_t value)
      {
      return _mm_castsi128_ps(_mm_set1_epi32(value));
      }

    inline std::string pnm_read_line(std::ifstream& file)
      {
      std::string line;
      std::getline(file, line);
      while (!line.empty() && line[0] == '#')
        std::getline(file, line);
      return line;
      }

    template <class T>
    inline void put_pixel(T* raw, const int x, const int y, const uint32_t stride, const T& clr)
      {
      raw[x + y * stride] = clr;
      }

    template <class T>
    inline void put_pixel_checked(T* raw, const int x, const int y, const uint32_t w, const uint32_t h, const uint32_t stride, const T& clr)
      {
      if (x >= 0 && y >= 0 && x < (int)w && y < (int)h)
        put_pixel(raw, x, y, stride, clr);
      }
    } // namespace details


  template <class T>
  inline image<T> span_to_image(uint32_t w, uint32_t h, uint32_t stride, const T* p_image)
    {
    assert(sizeof(T) <= 16);
    assert(16 % sizeof(T) == 0);
    image<T> out(w, h, false);
    uint32_t offset = 16 / sizeof(T);
    for (uint32_t y = 0; y < h; ++y)
      {
      const T* p_im = p_image + y * stride;
      T* p_out = out.data() + y * out.stride();
      for (uint32_t x = 0; x < w; x += offset, p_im += offset, p_out += offset)
        {
        __m128 data = _mm_loadu_ps((const float*)p_im);
        _mm_store_ps((float*)p_out, data);
        }
      }
    return out;
    }


  template <class T>
  void draw_line(T* raw, double x1, double y1, double x2, double y2, const uint32_t w, const uint32_t h, const uint32_t stride, const T& color)
    {
    using namespace details;
    // Bresenham's line algorith
    const bool steep = (std::abs(y2 - y1) > std::abs(x2 - x1));
    if (steep)
      {
      std::swap(x1, y1);
      std::swap(x2, y2);
      }

    if (x1 > x2)
      {
      std::swap(x1, x2);
      std::swap(y1, y2);
      }

    double dx = x2 - x1;
    double dy = std::abs(y2 - y1);

    double error = dx / 2.0;
    const int ystep = (y1 < y2) ? 1 : -1;
    int y = (int)y1;

    const int maxX = (int)x2;

    for (int x = (int)x1; x < maxX; ++x)
      {
      if (steep)
        {
        put_pixel_checked(raw, y, x, w, h, stride, color);
        }
      else
        {
        put_pixel_checked(raw, x, y, w, h, stride, color);
        }
      error -= dy;
      if (error < 0)
        {
        y += ystep;
        error += dx;
        }
      }
    }

  template <class T>
  void draw_box(T* raw, int x, int y, int box_width, int box_height, const uint32_t w, const uint32_t h, const uint32_t stride, const T& color)
    {
    using namespace details;

    int x0 = x;
    int y0 = y;
    int x1 = x + box_width - 1;
    int y1 = y + box_height - 1;

    if (x0 < 0)
      x0 = 0;
    if (x1 >= w)
      x1 = w - 1;
    if (y0 < 0)
      y0 = 0;
    if (y1 >= h)
      y1 = h - 1;

    for (int xx = x0; xx <= x1; ++xx)
      {
      put_pixel(raw, xx, y0, stride, color);
      put_pixel(raw, xx, y1, stride, color);
      }
    for (int yy = y0 + 1; yy < y1; ++yy)
      {
      put_pixel(raw, x0, yy, stride, color);
      put_pixel(raw, x1, yy, stride, color);
      }
    }


  template <class T>
  image<T> rotate90(const image<T>& im)
    {
    const auto w = im.width();
    const auto h = im.height();
    image<T> out(h, w, false);
    for (uint32_t x = 0; x < w; ++x)
      {
      for (uint32_t y = 0; y < h; ++y)
        {
        out(y, w - x - 1) = im(x, y);
        }
      }
    return out;
    }

  template <class T>
  image<T> rotate180(const image<T>& im)
    {
    const auto w = im.width();
    const auto h = im.height();
    image<T> out(w, h, false);
    for (uint32_t x = 0; x < w; ++x)
      {
      for (uint32_t y = 0; y < h; ++y)
        {
        out(w - x - 1, h - y - 1) = im(x, y);
        }
      }
    return out;
    }

  template <class T>
  image<T> rotate270(const image<T>& im)
    {
    const auto w = im.width();
    const auto h = im.height();
    image<T> out(h, w, false);
    for (uint32_t x = 0; x < w; ++x)
      {
      for (uint32_t y = 0; y < h; ++y)
        {
        out(h - y - 1, x) = im(x, y);
        }
      }
    return out;
    }

  inline bool load_pgm(image<uint16_t>& im, const std::string& filename)
    {
    using namespace details;
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file.is_open())
      return false;
    std::stringstream str;
    str << pnm_read_line(file);
    int width, height, max_val;
    width = height = max_val = -1;
    std::string P5;
    str >> P5;
    if (P5 != "P5")
      return false;
    str >> width;
    if (width == -1)
      {
      str.clear(); str.str("");
      str << pnm_read_line(file);
      str >> width;
      }
    str >> height;
    if (height == -1)
      {
      str.clear(); str.str("");
      str << pnm_read_line(file);
      str >> height;
      }
    str >> max_val;
    if (max_val == -1)
      {
      str.clear(); str.str("");
      str << pnm_read_line(file);
      str >> max_val;
      }
    if (max_val > 65535)
      return false;
    if (max_val <= 255)
      return false;
    im = image<uint16_t>(width, height, false);
    file.read((char *)im.data(), width * height * sizeof(uint16_t));
    file.close();
    return true;
    }

  template <class T>
  inline void fill_image(image<T>& im, T value)
    {
    using namespace details;
    const int w = (int)im.width();
    const int s = (int)im.stride();
    const int h = (int)im.height();
    assert(sizeof(T) <= 16);
    assert(16 % sizeof(T) == 0);
    uint32_t offset = 16 / sizeof(T);
    __m128 data = set_data_128(value);
    for (int y = 0; y < h; ++y)
      {
      T* p_im = im.data() + y * s;
      for (int x = 0; x < w; x += offset, p_im += offset)
        {
        _mm_store_ps((float*)p_im, data);
        }
      }
    }

  inline image<uint8_t> image_to_gray(const image<uint32_t>& im)
    {
    image<uint8_t> out(im.width(), im.height(), false);
    for (int y = 0; y < (int)im.height(); ++y)
      {
      const uint32_t* p_f = im.data() + im.stride()*y;
      uint8_t* p_g = out.data() + out.stride() * y;
      for (int x = 0; x < (int)im.width(); ++x, ++p_f, ++p_g)
        {
        *p_g = (*p_f >> 8) & 0xff;
        }
      }
    return out;
    }

  inline void minmax(float& min, float& max, const image<float>& im, bool skip_infinity)
    {
    using namespace details;
    const int w = (int)im.width();
    const int s = (int)im.stride();
    const int h = (int)im.height();

    float single_min = std::numeric_limits<float>::max();
    float single_max = -std::numeric_limits<float>::max();

    __m128 current_min = _mm_set1_ps(std::numeric_limits<float>::max());
    __m128 current_max = _mm_set1_ps(-std::numeric_limits<float>::max());

    __m128 inf = _mm_set1_ps(std::numeric_limits<float>::infinity());

    for (int y = 0; y < h; ++y)
      {
      __m128* p_im = (__m128*)(im.data() + y * s);
      int x = 0;
      for (; x < w - 4; x += 4)
        {
        __m128 val_min = *p_im;
        __m128 val_max = val_min;
        if (skip_infinity)
          {
          __m128 mask = _mm_cmpeq_ps(val_min, inf);
          val_min = _mm_blendv_ps(val_min, current_min, mask);
          val_max = _mm_blendv_ps(val_max, current_max, mask);
          }
        current_min = _mm_min_ps(current_min, val_min);
        current_max = _mm_max_ps(current_max, val_max);
        ++p_im;
        }
      const float* p_imf = im.data() + y * s + x;
      for (; x < w; ++x)
        {
        if (!skip_infinity || (*p_imf != std::numeric_limits<float>::infinity()))
          {
          single_min = std::min<float>(single_min, *p_imf);
          single_max = std::max<float>(single_max, *p_imf);
          }
        ++p_imf;
        }
      }

    min = std::min<float>(single_min, horizontal_min_ps(current_min));
    max = std::max<float>(single_max, horizontal_max_ps(current_max));
    }

  inline image<uint8_t> image_to_gray(const image<float>& im)
    {
    image<uint8_t> out(im.width(), im.height(), false);
    float fmin, fmax;
    minmax(fmin, fmax, im, true);
    float s = 1.f;
    if (fmin != fmax)
      s = 255.f / (fmax - fmin);
    for (int y = 0; y < (int)im.height(); ++y)
      {
      const float* p_f = im.data() + im.stride()*y;
      uint8_t* p_g = out.data() + out.stride() * y;
      for (int x = 0; x < (int)im.width(); ++x, ++p_f, ++p_g)
        {
        *p_g = uint8_t((*p_f - fmin)*s);
        }
      }
    return out;
    }


  template <class T, class TValid>
  void convolve_col_14641_div_16(image<T>& out, const image<T>& im, TValid valid)
    {
    int w = im.width();
    int h = im.height();
    out = image<T>(w, h);
    std::array<T, 5> vals;
    std::array<T, 5> coeff = { { T(1), T(4), T(6), T(4), T(1) } };
    for (int u = 0; u < w; ++u)
      {
      for (int v = 0; v < h; ++v)
        {
        if (!valid(im(u, v), u, v))
          continue;
        T denom = T(0);
        for (int i = 0; i < 5; ++i)
          {
          int vv = v + i - 2;
          if (vv >= 0 && vv < h && valid(im(u, vv), u, vv))
            {
            vals[i] = im(u, vv);
            denom += coeff[i];
            }
          else
            vals[i] = T(0);
          }
        out(u, v) = (vals[0] * coeff[0] + vals[1] * coeff[1] + vals[2] * coeff[2] + vals[3] * coeff[3] + vals[4] * coeff[4]) / denom;
        }
      }
    }

  template <class T, class TValid>
  void convolve_row_14641_div_16(image<T>& out, const image<T>& im, TValid valid)
    {
    int w = im.width();
    int h = im.height();
    out = image<T>(w, h);
    std::array<T, 5> vals;
    std::array<T, 5> coeff = { { T(1), T(4), T(6), T(4), T(1) } };
    for (int v = 0; v < h; ++v)
      {
      for (int u = 0; u < w; ++u)
        {
        if (!valid(im(u, v), u, v))
          continue;
        T denom = T(0);
        for (int i = 0; i < 5; ++i)
          {
          int uu = u + i - 2;
          if (uu >= 0 && uu < w && valid(im(uu, v), uu, v))
            {
            vals[i] = im(uu, v);
            denom += coeff[i];
            }
          else
            vals[i] = T(0);
          }
        out(u, v) = (vals[0] * coeff[0] + vals[1] * coeff[1] + vals[2] * coeff[2] + vals[3] * coeff[3] + vals[4] * coeff[4]) / denom;
        }
      }
    }

  template <class T, class TValid>
  void gauss(image<T>& im, TValid valid)
    {
    image<T> temp;
    convolve_col_14641_div_16(temp, im, valid);
    convolve_row_14641_div_16(im, temp, valid);
    }

  inline histogram_t make_histogram(const image<uint8_t>& im)
    {
    histogram_t hist(256, 0);
    for (auto val : im)
      ++hist[size_t(val)];
    return hist;
    }

  inline histogram_t make_histogram(const image<uint8_t>& im, int x, int y, int w, int h)
    {
    histogram_t hist(256, 0);
    if (x + w > (int)im.width())
      w = (int)im.width() - x;
    if (y + h > (int)im.height())
      h = (int)im.height() - y;
    for (int yy = 0; yy < h; ++yy)
      {
      const uint8_t* p_val = im.data() + (y + yy)*im.stride() + x;
      for (int xx = 0; xx < w; ++xx, ++p_val)
        {
        ++hist[size_t(*p_val)];
        }
      }
    return hist;
    }

  inline int64_t get_histogram_weight(histogram_t::const_iterator begin, histogram_t::const_iterator end)
    {
    int64_t w = 0;
    for (auto it = begin; it != end; ++it)
      w += *it;
    return w;
    }

  inline int otsu_threshold(const histogram_t& histogram)
    {
    int level = 125;
    int64_t total = get_histogram_weight(histogram.begin(), histogram.end());
    int64_t sumB = 0;
    int64_t wB = 0;
    double maximum = 0.0;
    int64_t sum1 = 0;
    for (int i = 0; i < histogram.size(); ++i)
      sum1 += i * histogram[i];
    for (int i = 0; i < histogram.size(); ++i)
      {
      wB += histogram[i];
      if (wB == 0)
        continue;
      int64_t wF = total - wB;
      if (wF == 0)
        break;
      sumB += i * histogram[i];
      double mB = double(sumB) / double(wB);
      double mF = double(sum1 - sumB) / double(wF);
      double between = wB * wF*(mB - mF)*(mB - mF);
      if (between >= maximum)
        {
        level = i;
        maximum = between;
        }
      }
    return level;
    }

  inline void threshold_image_region(image<uint8_t>& im, int x, int y, int w, int h, uint8_t threshold)
    {
    if (x + w > (int)im.width())
      w = (int)im.width() - x;
    if (y + h > (int)im.height())
      h = (int)im.height() - y;
    for (int yy = 0; yy < h; ++yy)
      {
      uint8_t* p_val = im.data() + (y + yy)*im.stride() + x;
      for (int xx = 0; xx < w; ++xx, ++p_val)
        {
        if (*p_val > threshold)
          *p_val = 255;
        else
          *p_val = 0;
        }
      }
    }

  inline image<uint8_t> make_binary_image(const image<uint8_t>& im, int steps_w, int steps_h)
    {
    image<uint8_t> out = im;

    const int w = (int)im.width();
    const int h = (int)im.height();

    const auto rect_w = w / steps_w;
    const auto rect_h = h / steps_h;

    const auto steps_x = w / rect_w + 1;
    const auto steps_y = h / rect_h + 1;

    for (int step_y = 0; step_y < steps_y; ++step_y)
      {
      for (int step_x = 0; step_x < steps_x; ++step_x)
        {
        const auto x = step_x * rect_w;
        const auto y = step_y * rect_h;
        auto hist = make_histogram(im, x, y, rect_w, rect_h);
        int ot = otsu_threshold(hist);
        uint8_t threshold = ot < 0 ? 0 : (ot > 255 ? 255 : (uint8_t)ot);
        threshold_image_region(out, x, y, rect_w, rect_h, threshold);
        }
      }
    return out;
    }

  inline image<uint8_t> make_binary_image(const image<uint8_t>& im, uint8_t threshold)
    {
    image<uint8_t> out(im.width(), im.height());
    const uint8_t* it = im.begin();
    const uint8_t* it_end = im.end();
    uint8_t* o = out.begin();
    for (; it != it_end; ++it, ++o)
      {
      if (*it > threshold)
        *o = 255;
      }
    return out;
    }

  inline image<uint32_t> three_gray_to_uint32_t(const image<uint8_t>& r, const image<uint8_t>& g, const image<uint8_t>& b)
    {
    assert(r.width() == g.width());
    assert(r.width() == b.width());
    assert(r.height() == g.height());
    assert(r.height() == b.height());
    image<uint32_t> out(r.width(), r.height());
    uint32_t* c = out.data();
    const uint8_t* pr = r.data();
    const uint8_t* pg = g.data();
    const uint8_t* pb = b.data();
    const uint32_t* const c_end = c + r.width()*r.height();
    for (; c != c_end; ++c, ++pr, ++pg, ++pb)
      {
      uint32_t clr = 0xff000000 | ((uint32_t)*pb << 16) | ((uint32_t)*pg << 8) | ((uint32_t)*pr);
      *c = clr;
      }
    return out;
    }


  inline image<uint32_t> census_transform(const image<uint8_t>& im)
    {
    /*
    The census transform is an image operator that associates to each pixel of a grayscale image a binary string,
    encoding whether the pixel has smaller intensity than each of its neighbours, one for each bit.
    We use a 9x7 window.
    */
    using namespace details;
    const int w = (int)im.width();
    const int s = (int)im.stride();
    const int h = (int)im.height();
    const int ns = (int)(sizeof(census_samples) / sizeof(int)) / 2;

    image<uint32_t> ct(w, h, true);

    const int cs = (int)ct.stride();

    for (int y = 3; y < h - 3; ++y)
      {
      __m128i descriptor1;
      __m128i descriptor2 = _mm_setzero_si128();
      const uint8_t* p_im = im.data() + s * y;
      for (int x = 4; x < w - 4; ++x, ++p_im)
        {
        uint32_t* ct_pos = ct.data() + (y*cs) + x;
        uint8_t* descriptor_it = (uint8_t*)&descriptor1;
        uint8_t value = *p_im;
        const __m128i value128 = _mm_set1_epi8(value);
        const int pos = y * s + x;
        for (int p = 0; p < 16; ++p)
          {
          *descriptor_it = *(p_im + census_samples[2 * p] + census_samples[2 * p + 1] * s);
          ++descriptor_it;
          }
        const __m128i r1 = _mm_cmpgt_epi8(descriptor1, value128);
        descriptor_it = (uint8_t*)&descriptor2;
        for (int p = 16; p < ns; ++p)
          {
          *descriptor_it = *(p_im + census_samples[2 * p] + census_samples[2 * p + 1] * s);
          ++descriptor_it;
          }
        const __m128i r2 = _mm_cmpgt_epi8(descriptor2, value128);

        *ct_pos = (_mm_movemask_epi8(r1)) | (_mm_movemask_epi8(r2) << 16);
        }
      }
    return ct;
    }

  inline void census_transform(image<uint32_t>& census_red, image<uint32_t>& census_green, image<uint32_t>& census_blue, const image<uint32_t>& im)
    {
    /*
    The census transform is an image operator that associates to each pixel of a grayscale image a binary string,
    encoding whether the pixel has smaller intensity than each of its neighbours, one for each bit.
    We use a 9x7 window.
    */
    using namespace details;
    const int w = (int)im.width();
    const int s = (int)im.stride();
    const int h = (int)im.height();
    const int ns = (int)(sizeof(census_samples) / sizeof(int)) / 2;

    census_red = image<uint32_t>(w, h, true);
    census_green = image<uint32_t>(w, h, true);
    census_blue = image<uint32_t>(w, h, true);

    const int cs = (int)census_red.stride();

    for (int y = 3; y < h - 3; ++y)
      {
      __m128i descriptor1_red;
      __m128i descriptor1_green;
      __m128i descriptor1_blue;
      __m128i descriptor2_red = _mm_setzero_si128();
      __m128i descriptor2_green = _mm_setzero_si128();
      __m128i descriptor2_blue = _mm_setzero_si128();
      const uint32_t* p_im = im.data() + s * y;
      for (int x = 4; x < w - 4; ++x, ++p_im)
        {
        uint32_t* ct_pos_red = census_red.data() + (y*cs) + x;
        uint32_t* ct_pos_green = census_green.data() + (y*cs) + x;
        uint32_t* ct_pos_blue = census_blue.data() + (y*cs) + x;
        uint8_t* descriptor_red_it = (uint8_t*)&descriptor1_red;
        uint8_t* descriptor_green_it = (uint8_t*)&descriptor1_green;
        uint8_t* descriptor_blue_it = (uint8_t*)&descriptor1_blue;
        uint32_t value = *p_im;
        const __m128i red128 = _mm_set1_epi8(value & 0xff);
        const __m128i green128 = _mm_set1_epi8((value >> 8) & 0xff);
        const __m128i blue128 = _mm_set1_epi8((value >> 16) & 0xff);
        const int pos = y * s + x;
        for (int p = 0; p < 16; ++p)
          {
          const uint32_t color = *(p_im + census_samples[2 * p] + census_samples[2 * p + 1] * s);
          *descriptor_red_it = color & 0xff;
          *descriptor_green_it = (color >> 8) & 0xff;
          *descriptor_blue_it = (color >> 16) & 0xff;
          ++descriptor_red_it;
          ++descriptor_green_it;
          ++descriptor_blue_it;
          }
        const __m128i r1 = _mm_cmpgt_epi8(descriptor1_red, red128);
        const __m128i g1 = _mm_cmpgt_epi8(descriptor1_green, green128);
        const __m128i b1 = _mm_cmpgt_epi8(descriptor1_blue, blue128);
        descriptor_red_it = (uint8_t*)&descriptor2_red;
        descriptor_green_it = (uint8_t*)&descriptor2_green;
        descriptor_blue_it = (uint8_t*)&descriptor2_blue;
        for (int p = 16; p < ns; ++p)
          {
          const uint32_t color = *(p_im + census_samples[2 * p] + census_samples[2 * p + 1] * s);
          *descriptor_red_it = color & 0xff;
          *descriptor_green_it = (color >> 8) & 0xff;
          *descriptor_blue_it = (color >> 16) & 0xff;
          ++descriptor_red_it;
          ++descriptor_green_it;
          ++descriptor_blue_it;
          }
        const __m128i r2 = _mm_cmpgt_epi8(descriptor2_red, red128);
        const __m128i g2 = _mm_cmpgt_epi8(descriptor2_green, green128);
        const __m128i b2 = _mm_cmpgt_epi8(descriptor2_blue, blue128);

        *ct_pos_red = _mm_movemask_epi8(r1) | (_mm_movemask_epi8(r2) << 16);
        *ct_pos_green = _mm_movemask_epi8(g1) | (_mm_movemask_epi8(g2) << 16);
        *ct_pos_blue = _mm_movemask_epi8(b1) | (_mm_movemask_epi8(b2) << 16);
        }
      }
    }

  } // namespace jtk
