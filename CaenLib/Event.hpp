#pragma once
#include "Hit.hpp"

namespace Caen1725
{
  class Event 
  {
   public:
        
    Event(bool handle_traces = false)
    {
      if (handle_traces) Colib::throw_error("Caen1725::Event can't handle traces (TBD)");
    }

    virtual ~Event()
    {
    }

    void push_back(Hit const & hit) noexcept
    {
      if (maxEvt <= mult) return;

      const auto i = mult; // In principle, this helps the compiler optimize
      label        [i] = hit.label         ;
      board_ID     [i] = hit.board_ID      ;
      channel_ID   [i] = hit.channel_ID    ;
      subchannel_ID[i] = hit.subchannel_ID ;
      adc          [i] = hit.adc           ;
      qlong        [i] = hit.qlong         ;
      caen_time    [i] = hit.caen_time     ;
      time         [i] = hit.time          ;
      rel_time     [i] = ((i == 0) ? 0 : static_cast<Int_t>(hit.time - time[0]));
      wfa_success  [i] = hit.wfa_success;

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

    Long64_t evtNb = 0;
    Int_t    mult  = 0;
    constexpr static inline Int_t maxEvt = 1000;
    UInt_t    label         [maxEvt];
    UShort_t  board_ID      [maxEvt];
    UChar_t   channel_ID    [maxEvt];
    UChar_t   subchannel_ID [maxEvt];
    Int_t     adc           [maxEvt];
    Int_t     qlong         [maxEvt];
    ULong64_t caen_time     [maxEvt];
    ULong64_t time          [maxEvt];
    Int_t     rel_time      [maxEvt];
    Bool_t    wfa_success   [maxEvt];
  };
}