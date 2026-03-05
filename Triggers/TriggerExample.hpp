#pragma once

// Change this line to generate the correct preprocessor message
#pragma message "TriggerExample.hpp loaded" // This is a normal preprocessor message

#include "../CaenLib/CaenRootEventBuilder.hpp"

/**
 * @brief Allows a used-defined trigger. Copy this definition in a new .hpp file named as you whish, 
 * e.g. MyTrigger.hpp, and call it at compile time with -DTRIGGER=MyTrigger. But do not rename the class !
 */
class Trigger
{
public:
  Trigger(Caen1725::EventBuilder * eventBuilder) : m_eventBuilder(eventBuilder) {}
  bool operator() (Caen1725::EventBuilder::EventId const & event_id) 
  {
    switch (typeUsed)
    {
      // Here is an example for triggering on at least one ring and one sector, using user-defined lookup-tables (very efficient, but may requires c++17)
      case TriggerType::AtLeastOneRingOneSector:
      {
        for (size_t i = 0; i<event_id.size(); ++i) 
        {
          auto const & hit_i = (*m_eventBuilder)[event_id[i]];
          if (m_hitIsRingLUT[hit_i.board_ID]) for (size_t j = 0; j<event_id.size(); ++j) 
          {
            auto const & hit_j = (*m_eventBuilder)[event_id[j]];
            if (m_hitIsSectorLUT[hit_j.board_ID]) return true;
          }
          return false; // This event don't have sector
        }
        return false; // This event don't have rings
      }

      // Here is an example for triggering on exactly one ring and one sector, using == operator (less efficient, maybe more user-friendly ?)
      case TriggerType::ExactlyOneRingOneSector:
      {
        int nbRings   = 0;
        int nbSectors = 0;
        for (size_t hit_i = 0; hit_i<event_id.size(); ++hit_i) 
        {
          auto const & hit = (*m_eventBuilder)[event_id[hit_i]];
               if (hit.board_ID == 6) ++nbRings ;
          else if (hit.board_ID == 7 
                || hit.board_ID == 8) ++nbSectors;
        }
        return (nbRings == 1 && nbSectors == 1);
      }
      case TriggerType::None: default: return true;
    }
  }

private:

  enum TriggerType {AtLeastOneRingOneSector, ExactlyOneRingOneSector, None};
  static constexpr int typeUsed = AtLeastOneRingOneSector;
  Caen1725::EventBuilder * m_eventBuilder = nullptr;
  auto static constexpr m_hitIsRingLUT = Colib::LUT<10000>([](UChar_t boardID)
  {
    return boardID == 6;
  });
  auto static constexpr m_hitIsSectorLUT = Colib::LUT<10000>([](UChar_t boardID)
  {
    return boardID == 7 || boardID == 8;
  });
};
