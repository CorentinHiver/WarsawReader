#ifndef BOARDAGGREGATE_HPP
#define BOARDAGGREGATE_HPP

#include "utils.hpp"
#include "ChannelAggregate.hpp"

namespace CaenDataReader1725
{  
  struct CheckBinMissed
  {} static checkBinMissed;

  class BoardAggregate
  {
  public:
    std::vector<uint32_t> buffer;

    const uint8_t check_bin_ref = 0b1010;
    bool handle_traces = true;

    size_t read_size = 0;
    static constexpr size_t header_size = 4*sizeof(tmp_u32);

    uint32_t  size              = 0; // Size in long words (32 bits)
    uint8_t   check_bin         = 0;
    uint8_t   DUAL_CHANNEL_MASK = 0;
    uint16_t  PATTERN           = 0;
    bool      BF                = false;
    uint8_t   BOARD_ID          = 0;
    uint32_t  COUNTER           = 0;
    uint32_t  TIME_TAG          = 0;

    ChannelMaskHelper<8> maskHelper;
    std::vector<ChannelAggregate> channels;

    int nbChannels() const noexcept {return maskHelper.nbBits();}

    bool readHeader(std::istream& data)
    {
      // 1. Size and Check bin :
      read_data(data, &tmp_u32); // Read the whole word
      debug("size, check_bin:", std::bitset<32>(tmp_u32));

      if (data.eof()) return false;

      if (getBitField (tmp_u32, 31, 28) != check_bin_ref) // Check whether the check bin is correct
      {
        error("Check bin missed : ", std::bitset<32>(tmp_u32)); 
        throw checkBinMissed;
      }

      this->clear();

      read_size += sizeof(tmp_u32); // To take into account the first word read before the clear

      size      = getBitField (tmp_u32, 27    ) * sizeof(uint32_t);
      check_bin = getBitField (tmp_u32, 31, 28);
      
      // 2. Miscellaneous :
      read_data(data, &tmp_u32, read_size); // Read the whole word
      debug("word n°2", std::bitset<32>(tmp_u32));
      DUAL_CHANNEL_MASK = getBitField (tmp_u32, 7     );
      PATTERN           = getBitField (tmp_u32, 22,  8);
      BF                = getBit      (tmp_u32, 26    );
      BOARD_ID          = getBitField (tmp_u32, 31, 27);

      // Handle the channel mask that allows to know the channel ID
      maskHelper = std::bitset<8>(DUAL_CHANNEL_MASK);

      // 3. Counter :
      read_data(data, &COUNTER, read_size);
      debug("COUNTER", std::bitset<32>(COUNTER));
      COUNTER &= mask(22);

      // 4. Timestamp : 
      read_data(data, &TIME_TAG, read_size);
      debug("TIME_TAG", std::bitset<32>(TIME_TAG));

      return true;
    }

    bool hasMoreChannels() const {return read_size < size;}

    void readChannelHeader(ChannelAggregate & channel, std::istream & data)
    {
      channel.handle_traces = handle_traces;
        
      channel.ID = maskHelper.getID();

      read_size += channel.readHeader(data);
    }

    bool newChannel(std::istream& data)
    {
      if (this -> hasMoreChannels())
      {
        channels.emplace_back();
        readChannelHeader(channels.back(), data);
        return true;
      }
      else return false;
    }

    bool readEvent(std::istream& data)
    {
      return (channels.back().readEvent(data, read_size));
    }

    bool skipEvent(std::istream& data)
    {
      return (channels.back().skipEvent(data, read_size));
    }

    bool newChannelAggregate(std::istream& data)
    {
      if (!this->newChannel(data)) return false;

      while(readEvent(data)) continue;

      return true;
    }

    void clear() noexcept 
    {
      read_size         = 0;
      size              = 0;
      check_bin         = 0;
      DUAL_CHANNEL_MASK = 0;
      PATTERN           = 0;
      BF                = false;
      BOARD_ID          = 0;
      COUNTER           = 0;
      TIME_TAG          = 0;
      maskHelper        = 0;
      channels.clear();
    }

    int nbEvents() const 
    {
      int ret = 0; 
      for (auto const & channel : channels) ret+=channel.events.size();
      return ret;
    }

    BoardAggregate& operator=(BoardAggregate const & other)
    {
      size              = other.size;
      read_size         = other.read_size;
      handle_traces     = other.handle_traces;
      check_bin         = other.check_bin;
      DUAL_CHANNEL_MASK = other.DUAL_CHANNEL_MASK;
      PATTERN           = other.PATTERN;
      BF                = other.BF;
      BOARD_ID          = other.BOARD_ID;
      COUNTER           = other.COUNTER;
      TIME_TAG          = other.TIME_TAG;
      maskHelper        = other.maskHelper;
      channels          = other.channels;
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& out, BoardAggregate const & board)
    {
      out << 
        " check_bin " << ((board.check_bin == board.check_bin_ref) 
                            ? (Colib::Color::BLUE + std::string("true") + Colib::Color::RESET) 
                            : (Colib::Color::RED  + std::string("false") + Colib::Color::RESET) ) <<
        " size "              << board.size                               <<
        " DUAL_CHANNEL_MASK " << std::bitset<8>(board.DUAL_CHANNEL_MASK)  <<
        " PATTERN "           << board.PATTERN                            <<
        " BF "                << nicer_bool(board.BF)                     <<
        " BOARD_ID "          << board.BOARD_ID                           <<
        " COUNTER "           << board.COUNTER                            <<
        " TIME_TAG "          << board.TIME_TAG                           <<
        " nb channels "       << board.nbChannels()                       <<
        " nb events "         << board.nbEvents();
      return out;
    }
  };
};
#endif //BOARDAGGREGATE_HPP
