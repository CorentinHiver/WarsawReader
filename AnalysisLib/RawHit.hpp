#ifndef RAWHIT_HPP
#define RAWHIT_HPP

#include "../CaenLib/BoardAggregate.hpp"

using ADC = uint32_t;
using Energy = double;
using Timestamp = uint64_t; // ps
using Label = uint8_t;


enum DetectorTypes{ EAGLE, BGO, NEDA, DSSDRing, DSSDSector, LaBr3, EMPTY};
std::vector<std::string> DetectorName = {"EAGLE", "BGO", "NEDA", "DSSDRing", "DSSDSector", "LaBr3", "EMPTY"};
static constexpr std::array<int, 10> Boards_map = {EAGLE, EAGLE, EMPTY, EMPTY, EMPTY, NEDA, DSSDRing, DSSDSector, DSSDSector, LaBr3};

/// @brief Detector type
constexpr inline auto getDetectorType(CaenDataReader::BoardAggregate const & board, CaenDataReader::CaenEvent const & event){
  auto detectorType = Boards_map[board.BOARD_ID];
  if (detectorType == EAGLE && event.CH) detectorType = BGO;
  return detectorType;
}

/// @brief detector label
constexpr inline auto labels(CaenDataReader::BoardAggregate const & board, CaenDataReader::ChannelAggregate const & channel, CaenDataReader::CaenEvent const & event) {
  switch (Boards_map[board.BOARD_ID]) {
  case EAGLE:                                   // ID is 0 to 16 : board.BOARD_ID [0;1], channel.ID [0,7]
    return board.BOARD_ID * 8 + channel.ID;

  case NEDA: case DSSDRing : case LaBr3:        // ID is 0 to 16 : channel.ID [0,7], event.CH [0;1]
    return channel.ID * 2 + int_cast(event.CH);
  
  case DSSDSector:                              // ID is 0 to 32 : board.BOARD_ID [7;8], channel.ID [0,7], event.CH [0,1]
    return (board.BOARD_ID - 7) * 16 + channel.ID * 2 + int_cast(event.CH);
  
  default:
    return -1;
  }
}

/// @brief detector global label
constexpr inline auto gLabels (CaenDataReader::BoardAggregate const & board, CaenDataReader::ChannelAggregate const & channel, CaenDataReader::CaenEvent const & event) {
  return 16 * getDetectorType(board, event) + labels(board, channel, event);
}


/// @brief Hit that extracts information from the CaenReader interface
/// @attention Works only for the implementation that was in place at least from 11/2024 to 05/2025
class RawHit
{
public:
  RawHit(){}
  RawHit(CaenDataReader::BoardAggregate const & board, CaenDataReader::ChannelAggregate const & channel, CaenDataReader::CaenEvent const & event) :
  label        (labels(board, channel, event)),
  glabel       (gLabels(board, channel, event)),
  adc          (event.ENERGY),
  detectorType (Boards_map[board.BOARD_ID]),
  samples      (event.samples)
  {
    if (event.extra.hasExtendedTimestamp()) 
    {
      extended_timestamp = event.extra.extended_timestamp;
      if (event.extra.flag == event.extra.FineTimestamp_flag)
      {
        precise_timestamp  = event.extra.precise_timestamp ;
      }
    }
    else timestamp = event.TRIGGER_TIME_TAG;

    if (detectorType == EAGLE && event.CH) detectorType = BGO;
    if (detectorType == NEDA) adc2 = event.EXTRAS;
  }

  Label      label              =  0;
  Label      glabel             =  0;
  ADC        adc                =  0;
  ADC        adc2               =  0;
  Energy     energy             =  0;
  int        detectorType       = -1;
  Timestamp  timestamp          =  0;
  Timestamp  extended_timestamp =  0;
  Timestamp  precise_timestamp  =  0;
  
  std::vector<CaenDataReader::Sample> samples;
  
  bool operator>(RawHit const & other) const {return precise_timestamp > other.precise_timestamp;}
};

/// @brief A RawHit buffer
template <size_t __size__>
class RawHitBuffer
{
public:
  RawHitBuffer(){}

  auto const & operator[](size_t const & i) const {return m_buffer[i];}
  auto         operator[](size_t const & i)       {return m_buffer[i];}

  auto         operator->()       {return &m_buffer;}
  auto const & operator->() const {return &m_buffer;}

  bool const & fill(CaenDataReader::BoardAggregate const & board) {
    if (full) {
      m_buffer.clear();
      full = false;
    }
    if (m_fill_size + board.nbEvents() > __size__) return (full = true);
    for (auto const & channel : board.channels) for (auto const & event : channel.events)
      m_buffer.emplace_back(board, channel, event);
    m_fill_size = 0;
    return full; // Returns false
  }
  
private:
  std::vector<RawHit> m_buffer;

  const size_t m_size = __size__;
  size_t m_fill_size = 0;
  bool full = false;
};


#endif //RAWHIT_HPP