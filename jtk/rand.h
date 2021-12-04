#pragma once

#include <stdint.h>

namespace jtk
  {

  class xorshift32
    {
    public:
      xorshift32() : _state(0x74382381) {}
      ~xorshift32() {}

      uint32_t operator()()
        {
        _state ^= (_state << 13);
        _state ^= (_state >> 17);
        _state ^= (_state << 5);
        return _state;
        }

      void seed(uint32_t s)
        {
        _state = s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        _state ^= s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        _state ^= s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        _state ^= s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        }

    private:
      uint32_t _state;
    };

  class xorshift64
    {
    public:
      xorshift64() : _state(7967805708226634297) {}
      ~xorshift64() {}

      uint64_t operator()()
        {
        _state ^= (_state << 13);
        _state ^= (_state >> 7);
        _state ^= (_state << 17);
        return _state;
        }

      void seed(uint64_t s)
        {
        _state = s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        _state ^= s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        _state ^= s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        _state ^= s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        }

    private:
      uint64_t _state;
    };

  class xorshift64star
    {
    public:
      xorshift64star() : _state(7967805708226634297) {}
      ~xorshift64star() {}

      uint64_t operator()()
        {
        _state ^= (_state >> 12);
        _state ^= (_state << 25);
        _state ^= (_state >> 27);
        return _state * 2685821657736338717LL;
        }

      void seed(uint64_t s)
        {
        _state = s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        _state ^= s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        _state ^= s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        _state ^= s + s * 17 + s * 121 + (s * 121 / 17);
        this->operator()();
        }

    private:
      uint64_t _state;
    };

  inline uint32_t hash32(int32_t position, uint32_t seed = 0)
    {
    constexpr uint32_t BIT_NOISE1 = 0xB5297A4D;
    constexpr uint32_t BIT_NOISE2 = 0x68E31DA4;
    constexpr uint32_t BIT_NOISE3 = 0x1B56C4E9;

    uint32_t mangled = (uint32_t)position;
    mangled *= BIT_NOISE1;
    mangled += seed;
    mangled ^= (mangled >> 8);
    mangled += BIT_NOISE2;
    mangled ^= (mangled << 8);
    mangled *= BIT_NOISE3;
    mangled ^= (mangled >> 8);
    return mangled;
    }

  inline uint64_t hash64(int64_t position, uint64_t seed = 0)
    {
    constexpr uint64_t BIT_NOISE1 = 0x3687C7A034A83D45;
    constexpr uint64_t BIT_NOISE2 = 0x275B2FC053370EA5;
    constexpr uint64_t BIT_NOISE3 = 0x74AF33B0E4606585;

    uint64_t mangled = (uint64_t)position;
    mangled *= BIT_NOISE1;
    mangled += seed;
    mangled ^= (mangled >> 16);
    mangled += BIT_NOISE2;
    mangled ^= (mangled << 16);
    mangled *= BIT_NOISE3;
    mangled ^= (mangled >> 16);
    return mangled;
    }

  inline uint32_t hash32_2d(int32_t pos_x, int32_t pos_y, uint32_t seed = 0)
    {
    constexpr uint32_t prime = 198491317;
    return hash32(pos_x + (prime * pos_y), seed);
    }

  inline uint32_t hash32_3d(int32_t pos_x, int32_t pos_y, int32_t pos_z, uint32_t seed = 0)
    {
    constexpr uint32_t prime1 = 198491317;
    constexpr uint32_t prime2 = 6542989;
    return hash32(pos_x + (prime1 * pos_y) + (prime2 * pos_z), seed);
    }

  inline uint64_t hash64_2d(int64_t pos_x, int64_t pos_y, uint64_t seed = 0)
    {
    constexpr uint64_t prime = 2563122999052471259LL;
    return hash64(pos_x + (prime * pos_y), seed);
    }

  inline uint64_t hash64_3d(int64_t pos_x, int64_t pos_y, int64_t pos_z, uint64_t seed = 0)
    {
    constexpr uint64_t prime1 = 2563122999052471259LL;
    constexpr uint64_t prime2 = 7156811723171476213LL;
    return hash64(pos_x + (prime1 * pos_y) + (prime2 * pos_z), seed);
    }

  class rndhash32
    {
    public:
      rndhash32() : _state(0), _seed(0) {}
      ~rndhash32() {}

      uint32_t operator()()
        {
        return hash32(_state++, _seed);
        }

      void seed(uint32_t s)
        {
        _seed = s;
        }

    private:
      uint32_t _state;
      uint32_t _seed;
    };

  class rndhash64
    {
    public:
      rndhash64() : _state(0), _seed(0) {}
      ~rndhash64() {}

      uint64_t operator()()
        {
        return hash64(_state++, _seed);
        }

      void seed(uint64_t s)
        {
        _seed = s;
        }

    private:
      uint64_t _state;
      uint64_t _seed;
    };

  } // namespace jtk
