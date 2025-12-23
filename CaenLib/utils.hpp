#ifndef CaenDataReader1725_UTILS_HPP
#define CaenDataReader1725_UTILS_HPP

#include "../LibCo/libCo.hpp"

namespace CaenDataReader1725
{
  //////////////////////
  // Helper Functions //
  //////////////////////

  inline constexpr uint64_t mask(uint bit_up) noexcept {return  (1 << (bit_up + 1)) - 1;}
  inline constexpr uint64_t mask(uint bit_up, uint bit_low) noexcept 
  {
    if (bit_low == 0) return mask(bit_up);
    else              return mask(bit_up) - mask(bit_low-1);
  }

  /// @brief Extracts [bit_up:0] from data (e.g. getBitField(data, 16) exctracts [16,0] for a bit field from bit 16 to bit 0)
  template<typename T> inline constexpr T getBitField(T const & data, uint const & bit_up) noexcept 
  {
  #ifdef DEBUG
    if (bit_up > sizeof(T)*8) Colib::throw_error("getBitField overflow");
  #endif // DEBUG
    return data & mask(bit_up);
  }

  /// @brief Extracts [bit_up:bit_low] from data (e.g. [31,16] for a bit field from bit 31 to bit 16)
  template<typename T> inline constexpr T getBitField(T const & data, uint const & bit_up, uint const & bit_low) noexcept 
  {
  #ifdef DEBUG
    if (bit_up > sizeof(T)*8) Colib::throw_error("getBitField overflow");
  #endif // DEBUG
    return (data & mask(bit_up, bit_low)) >> (bit_low);
  }

  template<typename T> inline constexpr bool getBit(T & data, uint const & bit) 
  {
  #ifdef DEBUG
    if (bit > sizeof(T)*8) Colib::throw_error("getBit overflow");
  #endif // DEBUG
    return bool_cast(getBitField(data, bit, bit));
  }

  ////////////////////
  // Helper Classes //
  ////////////////////

  template<int __size__ = 8>
  class ChannelMaskHelper
  {
  public:

    ChannelMaskHelper& operator=(std::bitset<__size__> const & _mask)
    {
      m_cursor = 0;
      m_mask = _mask;
      return *this;
    }

    int getID()
    {
      while (m_cursor < __size__) 
      {
        if (m_mask.test(m_cursor)) 
        {
          return m_cursor++;
        }
        else ++m_cursor;
      }
      error("ChannelMaskHelper::getID() : pas normal !!");
      return __size__; // We are at the end of the bitset, normally never reaching this part
    }

    auto nbBits() const noexcept {return m_mask.count();}

    friend std::ostream& operator<<(std::ostream& out, ChannelMaskHelper const & mask)
    {
      out << "mask : " << mask.m_mask << "; cursor currently at" << mask.m_cursor;
      return out;
    }

  private:
    int m_cursor = 0;
    std::bitset<__size__> m_mask;
  };

  thread_local uint8_t  tmp_u8  = 0 ;
  thread_local uint16_t tmp_u16 = 0 ;
  thread_local uint32_t tmp_u32 = 0 ;

  template<class Type>
  inline std::istream& read_data(std::istream& data, Type * buff) 
  {
  // #ifdef Cpp17
  //   static_assert(std::is_trivially_copyable_v<Type>);
  // #endif //Cpp17
    return data.read(reinterpret_cast<char*>(buff), sizeof(Type));
  }

  template<class Type>
  inline std::istream& read_data(std::istream& data, Type * buff, size_t & size_read) 
  {
  // #ifdef Cpp17
  //   static_assert(std::is_trivially_copyable_v<Type>);
  // #endif //Cpp17
    auto & ret = read_data(data, buff);
    if (ret) size_read += sizeof(Type);
    return ret;
  }

  inline void skip(std::istream& data, size_t size_to_skip)
  {
    data.seekg(static_cast<std::streamoff>(size_to_skip), std::ios::cur);
  }
  inline void skip(std::istream& data, size_t size_to_skip, size_t & size_read) {
    skip(data, size_to_skip);
    size_read += size_to_skip;
  }

  inline constexpr double ticks_to_ns = 4; // ns
  inline constexpr double ticks_to_ps = ticks_to_ns * 1000.; // ps
}


#endif //CaenDataReader1725_UTILS_HPP
