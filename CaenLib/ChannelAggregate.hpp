#ifndef CHANNELAGGREGATE_HPP
#define CHANNELAGGREGATE_HPP

#include "utils.hpp"
#include "CaenEvent.hpp"

namespace Caen1725
{
  struct ChannelAggregate
  {
    uint32_t  size        = 0;
    bool      FI          = false;
    uint16_t  NUM_SAMPLES = 0;
    uint8_t   DP          = 0;
    uint8_t   AP2         = 0;
    uint8_t   AP1         = 0;
    uint8_t   EX          = 0;
    bool      ES          = false;
    bool      E2          = false;
    bool      ET          = false;
    bool      EE          = false;
    bool      DT          = false;
    
    bool    handle_traces = true;
    size_t  read_size     = 0;
    uint8_t ID            = -1;

    std::vector<CaenEvent> events;

    ChannelAggregate() noexcept = default;

    bool hasMoreEvents() {return read_size < size;}

    // Header size + nb_evnts*size_event
    // size_event = timestamp + nb_samples * 2 octets + EXTRAS2 + EXTRA
    double calculated_size() {return 8+events.size()*(4+NUM_SAMPLES*2+8);}

    void clear() {
      size        = 0;
      FI          = false;
      NUM_SAMPLES = 0;
      DP          = 0;
      AP2         = 0;
      AP1         = 0;
      EX          = 0;
      ES          = false;
      E2          = false;
      ET          = false;
      EE          = false;
      DT          = false;
      ID          = -1;
      read_size   = 0;
      events.clear();
    }

    size_t readHeader(std::istream& data)
    {
      // 1. Size and FI
      read_data(data, &tmp_u32, read_size);
      debug("size, FI:", std::bitset<32>(tmp_u32));
      size = getBitField(tmp_u32, 30) * sizeof(uint32_t);
      FI   = getBit     (tmp_u32, 31);

      if (FI)
      {
        // 2. Format
        read_data(data, &tmp_u32, read_size);
        debug("Format:", std::bitset<32>(tmp_u32));
        NUM_SAMPLES = getBitField (tmp_u32, 15) * 8; // The data stored is NUM_SAMPLES/8 (see doc)
        DP          = getBitField (tmp_u32, 19, 16);
        AP2         = getBitField (tmp_u32, 21, 20);
        AP1         = getBitField (tmp_u32, 23, 22);
        EX          = getBitField (tmp_u32, 26, 24);

        ES = getBit (tmp_u32, 27);
        E2 = getBit (tmp_u32, 28);
        ET = getBit (tmp_u32, 29);
        EE = getBit (tmp_u32, 30);
        DT = getBit (tmp_u32, 31);
      }
      return read_size;
    }
    
    bool readEvent(std::istream& data, size_t & board_size)
    {
      if (hasMoreEvents())
      {
        auto const pos_before = data.tellg();

        events.emplace_back(EX);
        events.back().read(data, NUM_SAMPLES, handle_traces);
        
        auto const & read_length = int_cast(data.tellg() - pos_before);
        read_size  += read_length;
        board_size += read_length;
  
        return true;
      }
      else return false;
    }

    bool skipEvent(std::istream& data, size_t & board_size)
    {
      if (hasMoreEvents())
      {
        auto const pos_before = data.tellg();

        events.emplace_back(EX);
        events.back().skip(data, NUM_SAMPLES);
        
        auto const & read_length = int_cast(data.tellg() - pos_before);
        read_size  += read_length;
        board_size += read_length;
  
        return true;
      }
      else return false;
    }

    friend std::ostream& operator<<(std::ostream& out, ChannelAggregate const & channel) {
      out << 
        " ID "           << channel.ID                  <<
        " size "         << channel.size                <<
        " FI "           << nicer_bool(channel.FI)      <<
        " NUM_SAMPLES "  << channel.NUM_SAMPLES         <<
        " DP "           << std::bitset<4>(channel.DP ) <<
        " AP2 "          << std::bitset<2>(channel.AP2) <<
        " AP1 "          << std::bitset<2>(channel.AP1) <<
        " EX "           << std::bitset<3>(channel.EX ) <<
        " ES "           << nicer_bool(channel.ES)      <<
        " E2 "           << nicer_bool(channel.E2)      <<
        " ET "           << nicer_bool(channel.ET)      <<
        " EE "           << nicer_bool(channel.EE)      <<
        " DT "           << nicer_bool(channel.DT)      <<
        " nb events : "  << channel.events.size();
      return out;
    }
  };
};
#endif //CHANNELAGGREGATE_HPP



// ChannelAggregate(ChannelAggregate const & other):
//   size           (other.size),
//   FI             (other.FI),
//   NUM_SAMPLES    (other.NUM_SAMPLES),
//   DP             (other.DP),
//   AP2            (other.AP2),
//   AP1            (other.AP1),
//   EX             (other.EX),
//   ES             (other.ES),
//   E2             (other.E2),
//   ET             (other.ET),
//   EE             (other.EE),
//   DT             (other.DT),
//   handle_traces  (other.handle_traces),
//   ID             (other.ID),
//   events         (other.events)
// {
// }

// ChannelAggregate& operator=(ChannelAggregate const & other)
// {
//   size           = other.size;
//   FI             = other.FI;
//   NUM_SAMPLES    = other.NUM_SAMPLES;
//   DP             = other.DP;
//   AP2            = other.AP2;
//   AP1            = other.AP1;
//   EX             = other.EX;
//   ES             = other.ES;
//   E2             = other.E2;
//   ET             = other.ET;
//   EE             = other.EE;
//   DT             = other.DT;
//   handle_traces  = other.handle_traces;
//   ID             = other.ID;
//   events         = other.events;
//   return *this;
// }