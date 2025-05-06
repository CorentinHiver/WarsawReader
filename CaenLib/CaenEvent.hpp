#ifndef CAENEVENT_HPP
#define CAENEVENT_HPP

#include "utils.hpp"
#include "Sample.hpp"
#include "Extra2.hpp"

namespace CaenDataReader
{
  struct CaenEvent
  {
    uint32_t  TRIGGER_TIME_TAG  = 0;
    bool      CH                = false;
    uint32_t  EXTRAS2           = 0;
    uint16_t  ENERGY            = 0;
    bool      PU                = false;
    uint16_t  EXTRAS            = 0;
    uint8_t   EX                = 0; // ChannelAggregate information, duplicate for convenience
  
    std::vector<Sample> samples; // Sampled waveform
    Extra2 extra; // Extra informations like extended timestamp, trigger count or zero crossing

    CaenEvent() noexcept = default;

    CaenEvent (uint8_t const & _EX) : EX(_EX) {}

    void read(std::istream& data, int const & nb_samples, bool const & handle_traces)
    {
      read_buff(&tmp_u32, data);
      debug("TRIGGER_TIME_TAG, CH:", std::bitset<32>(tmp_u32));
      
      TRIGGER_TIME_TAG = getBitField(tmp_u32, 30);
      CH               = getBit     (tmp_u32, 31);
      
      // Looping through the NUM_SAMPLES samples of the waveform of each event :
      if (handle_traces) {
        debug("Trace with :", nb_samples, "samples");
        samples.resize(nb_samples);
        for (auto & sample : samples)
        {
          read_buff(&tmp_u16, data);
          sample.sample  = getBitField(tmp_u16, 13);
          sample.DP1     = getBit     (tmp_u16, 14);
          sample.T       = getBit     (tmp_u16, 15);
        }
      }
      else {
        auto const & size_to_skip = nb_samples * sizeof(nb_samples);
        skip(data, size_to_skip);
      }
  
      // if (channel.E2)  // TODO: do we need to read EXTRAS2 if it is disabled, i.e. if channel.E2=false ?
        read_buff(&EXTRAS2, data);
      debug("EXTRAS2", std::bitset<32>(EXTRAS2));
  
      extra = Extra2(EXTRAS2, TRIGGER_TIME_TAG, EX);
      
      read_buff(&tmp_u32, data);
  
      debug("ENERGY, PU, EXTRAS:", std::bitset<32>(tmp_u32));
  
      ENERGY = getBitField (tmp_u32, 14    );
      PU     = getBit      (tmp_u32, 15    );
      EXTRAS = getBitField (tmp_u32, 25, 16);
    }
  
    void read(std::istream& data, uint8_t const & _EX, int const & nb_samples, bool const & handle_traces)
    {
      EX = _EX;
      return read(data, nb_samples, handle_traces);
    }

    void clear() 
    {
      TRIGGER_TIME_TAG  = 0;
      CH                = false;
      EXTRAS2           = 0;
      ENERGY            = 0;
      PU                = false;
      EXTRAS            = 0;
      EX                = 0;
      extra             = Extra2();
      samples.clear();
    }
    
    friend std::ostream& operator<<(std::ostream& out, CaenEvent const & event) 
    {
      out << 
        " TRIGGER_TIME_TAG " << event.TRIGGER_TIME_TAG  <<
        " CH "               << nicer_bool(event.CH)    <<
        " EXTRAS2 "          << event.EXTRAS2           <<
        " ENERGY "           << event.ENERGY            <<
        " PU "               << nicer_bool(event.PU)    <<
        " EXTRAS "           << std::bitset<10>(event.EXTRAS) <<
        " extra "            << event.extra;
      return out;
    }
  };
};
#endif //CAENEVENT_HPP

// CaenEvent(CaenEvent const & other) :
//   TRIGGER_TIME_TAG  (other.TRIGGER_TIME_TAG),
//   CH                (other.CH),
//   EXTRAS2           (other.EXTRAS2),
//   ENERGY            (other.ENERGY),
//   PU                (other.PU),
//   EXTRAS            (other.EXTRAS),
//   EX                (other.EX),
//   samples           (other.samples),
//   extra             (other.extra)
// {
// }

// CaenEvent& operator=(CaenEvent const & other)
// {
//   TRIGGER_TIME_TAG  = other.TRIGGER_TIME_TAG;
//   CH                = other.CH;
//   EXTRAS2           = other.EXTRAS2;
//   ENERGY            = other.ENERGY;
//   PU                = other.PU;
//   EXTRAS            = other.EXTRAS;
//   EX                = other.EX;
//   samples           = other.samples;
//   extra             = other.extra;

//   return *this;
// }
