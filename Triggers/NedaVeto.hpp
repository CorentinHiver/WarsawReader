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
  Trigger(Caen1725::EventBuilder * eventBuilder) : m_eventBuilder(eventBuilder) {}
  bool operator() (Caen1725::EventBuilder::EventId const & event_id) 
  {
    // Here is an example to reject events with NEDA hits inside
    for (size_t i = 0; i<event_id.size(); ++i) 
    {
      auto const & hit = (*m_eventBuilder)[event_id[i]]; // Simple aliasing
      if (hit.board_ID == 5 && hit.channel_ID < 5) return false;
    }
    return true;
  }

private:
  Caen1725::EventBuilder * m_eventBuilder = nullptr;
};
