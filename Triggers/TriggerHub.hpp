#pragma once

#ifdef TRIGGER

#pragma message ("TRIGGER CHOSE TO BE " STRINGIFY(TRIGGER)) // This is a normal preprocessor message

#define STR_IMPL(x) #x
#define STR(x) STR_IMPL(x)

#define HEADER(x) STR(x.hpp)

// /**
//  * @brief Allows a used-defined trigger. Copy this definition in a new .hpp file named as you whish, 
//  * e.g. MyTrigger.hpp, and call it with -DTRIGGER=MyTrigger.
//  */
// class TriggerTemplate
// {
// public:
//   TriggerTemplate(CaenRootEventBuilder * eventBuilder) : m_eventBuilder(eventBuilder) {}
//   bool operator(CaenRootEventBuilder::EventId const & event_id) 
//   {
//     // Example : DSSD labels are defined true in the m_lut (lookup table)
//     // Here is an example for triggering on at least one ring and one sector
//     for (int i = 0; i<event_id.size(); ++i) 
//     {
//       auto const & label_i = (*m_eventBuilder)[event_id[i]].label;
//       if (m_RingLUT[label_i]) for (int j = 0; j<event_id.size(); ++j) 
//       {
//          auto const & label_j = (*m_eventBuilder)[event_id[j]];
//          if (m_SectorLUT[label_j]) return true;
//       }
//     }
//     return false;
//   }
// private:
//   CaenRootEventBuilder * m_eventBuilder = nullptr;
//   auto constexpr m_RingLUT = Colib::LUT<10000>([](UShort_t label)
//   {
//     return 6*16 < label <= 7*16;
//   });
//   auto constexpr m_SectorLUT = Colib::LUT<10000>([](UShort_t label)
//   {
//     return 7*16 < label <= 7*16;
//   });
// };

#include HEADER(TRIGGER) // This line triggers an error in VSCode but compiles just fine

#endif //TRIGGER