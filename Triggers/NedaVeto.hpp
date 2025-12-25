#pragma once

#pragma message (STRINGIFY(TRIGGER) ".hpp loaded") // This is a normal preprocessor message

#include "../CaenLib/CaenRootEventBuilder.hpp"

/**
 * @brief Allows a used-defined trigger. Copy this definition in a new .hpp file named as you whish, 
 * e.g. MyTrigger.hpp, and call it with -DTRIGGER=MyTrigger. But do not rename the class !
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
      auto const & hit = (*m_eventBuilder)[event_id[i]]; // Simple aliasing
      if (hit.board_ID == 5 && hit.channel_ID < 5) return false;
    }
    return true;
  }

private:
  Caen1725EventBuilder * m_eventBuilder = nullptr;
};
