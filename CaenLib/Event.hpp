#pragma once
#include "Hit.hpp"

namespace Caen1725
{
  using EventID   = Long64_t; // Event ID
  using EventMult = Int_t   ; // Event multiplicity

  class Event 
  {
   public:
   
    EventID   eventID = 0;
    EventMult mult    = 0;

    constexpr static inline size_t maxSize = 1000;
    Label     label         [maxSize];
    BoardID   board_ID      [maxSize];
    ChannelID channel_ID    [maxSize];
    ChannelID subchannel_ID [maxSize];
    ADC       adc           [maxSize];
    ADC       qlong         [maxSize];
    Timestamp caen_time     [maxSize];
    Timestamp time          [maxSize];
    Time_rel  rel_time      [maxSize];
    Bool_t    wfa_success   [maxSize];
        
    Event(bool handle_traces = false)
    {
      if (handle_traces) Colib::throw_error("Caen1725::Event can't handle traces (TBD)");
    }

    virtual ~Event()
    {
    }

    void push_back(Hit const & hit) noexcept
    {
      if (maxSize <= static_cast<size_t>(mult)) return;

      auto const i = mult; // In principle, this helps the compiler optimize

      label         [i] = hit.label         ;
      board_ID      [i] = hit.board_ID      ;
      channel_ID    [i] = hit.channel_ID    ;
      subchannel_ID [i] = hit.subchannel_ID ;
      adc           [i] = hit.adc           ;
      qlong         [i] = hit.qlong         ;
      caen_time     [i] = hit.caen_time     ;
      time          [i] = hit.time          ;
      rel_time      [i] = ((i == 0) ? 0 : static_cast<Int_t>(hit.time - time[0]));
      wfa_success   [i] = hit.wfa_success;

      ++mult;
    }

    void clear() {mult = 0;}
    
    auto size () const {return static_cast<size_t>(mult);}
    auto const & multiplicity () const {return mult;}

    Hit operator[](size_t i) const
    {
      return Hit(
        label         [i],
        board_ID      [i],
        channel_ID    [i],
        subchannel_ID [i],
        adc           [i],
        qlong         [i],
        caen_time     [i],
        time          [i],
        rel_time      [i],
        wfa_success   [i]
      );
    }

    constexpr auto dT(int index1, int index2) const noexcept {return Long64_t(time[index1]-time[index2]);}
  };
}