#pragma once
#include "Hit.hpp"

namespace CaenDataReader1725
{
  class Event 
  {
   public:
        
    Event(bool handle_traces = false)
    {
      if (handle_traces) Colib::throw_error("CaenDataReader1725::Event can't handle traces (TBD)");
    }

    virtual ~Event()
    {
    }

    void push_back(Hit const & hit) noexcept
    {
      const auto i = mult; // In principle, this helps the compiler optimize
      label        [i] = hit.label         ;
      board_ID     [i] = hit.board_ID      ;
      channel_ID   [i] = hit.channel_ID    ;
      subchannel_ID[i] = hit.subchannel_ID ;
      adc          [i] = hit.adc           ;
      qlong        [i] = hit.qlong         ;
      timestamp    [i] = hit.timestamp     ;
      time         [i] = hit.time          ;
      rel_time     [i] = ((i == 0) ? 0 : static_cast<Int_t>(hit.time - time[0]));

      ++mult;
    }

    void clear()
    {
      mult = 0;
    }

    
    auto const & size () const {return mult;}

    Hit operator[](size_t i) const
    {
      return Hit(
        label         [i],
        board_ID      [i],
        channel_ID    [i],
        subchannel_ID [i],
        adc           [i],
        qlong         [i],
        timestamp     [i],
        time          [i],
        rel_time      [i]
      );
    }


    size_t evtNb = 0;
    int    mult  = 0;
    constexpr static inline size_t maxEvt = 1000;
    Int_t     label         [maxEvt];
    Int_t     board_ID      [maxEvt];
    Int_t     channel_ID    [maxEvt];
    Int_t     subchannel_ID [maxEvt];
    Int_t     adc           [maxEvt];
    Int_t     qlong         [maxEvt];
    ULong64_t timestamp     [maxEvt];
    ULong64_t time          [maxEvt];
    ULong64_t rel_time      [maxEvt];
  };
}