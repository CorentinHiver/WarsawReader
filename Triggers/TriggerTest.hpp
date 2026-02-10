#pragma once

#pragma message "TriggerTest.hpp loaded" // This is a normal preprocessor message

#include "../CaenLib/CaenRootEventBuilder.hpp"

/**
 * @brief Allows a used-defined trigger. Copy this definition in a new .hpp file named as you whish, 
 * e.g. MyTrigger.hpp, and call it at compile time with -DTRIGGER=MyTrigger. But do not rename the class !
 */
class Trigger
{
public:
  Trigger(Caen1725EventBuilder * eventBuilder) : m_eventBuilder(eventBuilder) {}
  bool operator() (Caen1725EventBuilder::EventId const & event_id) 
  {
    // Example : DSSD labels are defined true in the m_lut (lookup table)
    // Here is an example for triggering on
    for (size_t i = 0; i<event_id.size(); ++i) 
    {
      auto const & boardID_i = (*m_eventBuilder)[event_id[i]].board_ID; // Simple aliasing
      if (m_RingLUT[boardID_i]) for (size_t j = 0; j<event_id.size(); ++j) 
      {
        auto const & boardID_j = (*m_eventBuilder)[event_id[j]].board_ID; // Simple aliasing
        if (m_SectorLUT[boardID_j]) return true;
      }
    }
    return false;
  }

private:
  Caen1725EventBuilder * m_eventBuilder = nullptr;
  auto static constexpr m_RingLUT = Colib::LUT<10000>([](UChar_t boardID)
  {
    return boardID == 6;
  });
  auto static constexpr m_SectorLUT = Colib::LUT<10000>([](UChar_t boardID)
  {
    return boardID == 6 || boardID == 7;
  });
};
