#ifndef EXTRA2_HPP
#define EXTRA2_HPP

#include "utils.hpp"

namespace CaenDataReader
{
  struct Extra2
  {
    constexpr static uint8_t FineTimestamp_flag     = 0b010;
    constexpr static uint8_t ExtendedTimestamp_flag = 0b000;
    constexpr static uint8_t TriggerCount_flag      = 0b100;
    constexpr static uint8_t ZeroCrossing_flag      = 0b101;
    uint8_t flag = 0;
    
    Extra2() = default;
    
    constexpr bool hasExtendedTimestamp() const {return bool((flag & 0b100) >> 2);}

    Extra2(uint32_t const & word, uint32_t const & timestamp, uint8_t const & _flag) : 
      flag(_flag)
    {                  
      switch (flag)
      {
        case FineTimestamp_flag:
          extended_timestamp  = getBitField(word, 31, 16);
          fine_timestamp      = getBitField(word, 9, 0) * ticks_to_ns / 1.024;
          
          // Convert to ps :
          extended_timestamp  = ((extended_timestamp << 32) | timestamp) * ticks_to_ps; // TODO : check if the shift is correct
          precise_timestamp   = extended_timestamp + fine_timestamp;
          break;
          
        case ExtendedTimestamp_flag:
          extended_timestamp  = getBitField(word, 31, 16);
          trapezoid_baseline = getBitField(word, 15, 0) / 4;

          // Convert to ps :
          extended_timestamp  = ((extended_timestamp << 32) | timestamp) * ticks_to_ps; // TODO : check if the shift is correct
          break;

        case TriggerCount_flag:
          lost_trigger_counter  = getBitField(word, 31, 16);
          total_trigger_counter = getBitField(word, 15,  0);
          break;
          
        case ZeroCrossing_flag:
          event_before_zero_crossing = getBitField(word, 31, 16);
          event_after_zero_crossing  = getBitField(word, 15,  0);
          break;
        
      }
    }

    uint64_t extended_timestamp         = 0; // in ps
    uint16_t fine_timestamp             = 0; // in ps
    uint64_t precise_timestamp          = 0; // in ps
    uint16_t trapezoid_baseline         = 0;
    uint16_t lost_trigger_counter       = 0;
    uint16_t total_trigger_counter      = 0;
    uint16_t event_before_zero_crossing = 0;
    uint16_t event_after_zero_crossing  = 0;

    friend std::ostream& operator<< (std::ostream& out, Extra2 const & extra)
    {
      switch (extra.flag)
      {
        case Extra2::FineTimestamp_flag:{
          out << " extended_timestamp " << extra.extended_timestamp
              << " fine_timestamp "     << extra.fine_timestamp
              << " precise_timestamp "  << extra.precise_timestamp;
          break;
        }
        case Extra2::ExtendedTimestamp_flag:{
          out << " extended_timestamp " << extra.extended_timestamp
              << " trapezoid_baseline " << extra.trapezoid_baseline;
          break;
        }
        case Extra2::TriggerCount_flag:{
          out << " lost_trigger_counter "  << extra.lost_trigger_counter
              << " total_trigger_counter " << extra.total_trigger_counter;
          break;
        }
        case Extra2::ZeroCrossing_flag:{
          out << " event_before_zero_crossing "  << extra.event_before_zero_crossing
              << " event_after_zero_crossing "   << extra.event_after_zero_crossing;
          break;
        }
        default:
          out << "EXTRA2 " << std::bitset<3>(extra.flag) << " not handled";
      }
      return out;
    }
  };
  
  // struct FineTimestamp : public Extra2
  // {
    
  
  //   FineTimestamp(uint32_t const & word, uint32_t const & timestamp) : 
  //   Extra2(),
  //   extended_timestamp (getBitField(word, 31, 16)), 
  //   fine_timestamp     (getBitField(word, 9 ,  0) * ticks_to_ns / 1.024)
  //   {                  
  //     extended_timestamp  = ((extended_timestamp << 32) | timestamp) * ticks_to_ps; // TODO : check if the shift is correct
  //     precise_timestamp   = extended_timestamp + fine_timestamp;
  //   }

  //   FineTimestamp* clone() const override {return new FineTimestamp(*this);}
  
  //   uint64_t extended_timestamp   = 0; // in ps
  //   uint16_t const fine_timestamp = 0; // in ps
  //   uint64_t precise_timestamp    = 0; // in ps
  // };
  
  // struct ExtendedTimestamp : public Extra2
  // {
  //   constexpr static uint8_t flag = 0b000;

  //   ExtendedTimestamp(uint32_t const & word, uint32_t const & timestamp) : 
  //     Extra2(),
  //     extended_timestamp (getBitField(word, 31, 16)),
  //     trapezoid_baseline (getBitField(word, 15,  0) / 4)
  //   {
  //     extended_timestamp  = ((extended_timestamp << 32) | timestamp) * ticks_to_ps; // TODO : check if the shift is correct
  //   }

  //   ExtendedTimestamp* clone() const override {return new ExtendedTimestamp(*this);}

  //   uint64_t extended_timestamp = 0;
  //   uint16_t const trapezoid_baseline  = 0;
  // };
  
  // struct TriggerCount : public Extra2
  // {
  //   constexpr static uint8_t flag = 0b100;

  //   TriggerCount(uint32_t const & word) : Extra2(),
  //     lost_trigger_counter  (getBitField(word, 31, 16)),
  //     total_trigger_counter (getBitField(word, 15,  0))
  //   {}

  //   TriggerCount* clone() const override {return new TriggerCount(*this);}

  //   uint16_t const lost_trigger_counter  = 0;
  //   uint16_t const total_trigger_counter = 0;
  // };
  
  // struct ZeroCrossing : public Extra2
  // {
  //   constexpr static uint8_t flag = 0b101;

  //   ZeroCrossing(uint32_t const & word) : Extra2(),
  //     event_before_zero_crossing  (getBitField(word, 31, 16)),
  //     event_after_zero_crossing   (getBitField(word, 15,  0))
  //   {}

  //   ZeroCrossing* clone() const override {return new ZeroCrossing(*this);};

  //   uint16_t const event_before_zero_crossing = 0;
  //   uint16_t const event_after_zero_crossing  = 0;
  // };


};


#endif //EXTRA2_HPP