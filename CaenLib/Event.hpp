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

    virtual void push_back(Hit const & hit)
    {
      if (maxEvt < static_cast<size_t>(mult)) error("Event too big :", mult, ">", maxEvt);

      label        [mult] = hit.label         ;
      board_ID     [mult] = hit.board_ID      ;
      channel_ID   [mult] = hit.channel_ID    ;
      subchannel_ID[mult] = hit.subchannel_ID ;
      adc          [mult] = hit.adc           ;
      qlong        [mult] = hit.qlong         ;
      timestamp    [mult] = hit.timestamp     ;
      time         [mult] = hit.time          ;
      rel_time     [mult] = ((mult == 0) ? 0 : static_cast<Int_t>(hit.time - time[0]));

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