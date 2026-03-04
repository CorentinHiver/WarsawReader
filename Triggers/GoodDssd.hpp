#pragma once

#pragma message (STRINGIFY(TRIGGER) ".hpp loaded") // This is a normal preprocessor message

#include "../CaenLib/CaenRootEventBuilder.hpp"

/**
 * @brief Allows a used-defined trigger. Copy this definition in a new .hpp file named as you whish, 
 * e.g. MyTrigger.hpp, and call it with -DTRIGGER=MyTrigger, or in the Makefile fill TRIGGER=MyTrigger. But do not rename the class !
 */
class Trigger
{
public:
  Trigger(Caen1725EventBuilder * eventBuilder) : m_eventBuilder(eventBuilder) {}
  bool operator() (Caen1725EventBuilder::EventId const & event_id) 
  {
    // Trigger on "good DSSD" events, i.e. exactly one ring and one sector in the event 
    m_nbRings = 0;
    m_nbSectors = 0;
    for (size_t i = 0; i<event_id.size(); ++i) 
    {
      auto const & hit = (*m_eventBuilder)[event_id[i]]; // Simple aliasing
      if (6 == hit.board_ID) ++m_nbRings;
      else if (7 == hit.board_ID || 8 == hit.board_ID) ++m_nbSectors;
    }
    return m_nbRings==1 && m_nbSectors==1;
  }

private:
  int m_nbRings = 0;
  int m_nbSectors = 0;
  Caen1725EventBuilder * m_eventBuilder = nullptr;
};
