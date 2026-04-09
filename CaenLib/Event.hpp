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
      traces       [i] = hit.trace;

      ++mult;
    }

    void emplace_back(Hit && hit) noexcept
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
      traces       [i] = std::move(hit.trace);

      ++mult;

      hit.clear();
    }

    void clear() {mult = 0;}
    
    auto size () const {return static_cast<size_t>(mult);}
    auto const & multiplicity () const {return mult;}

    Hit getHit(size_t const i) const
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
        wfa_success   [i],
        traces        [i]
      );
    }

    Hit operator[](size_t const i) const
    {
      return getHit(i);
    }

    /// @brief More user-friendly interface : copies the data inside a vector of hits.
    /// @details This method can be very expensive in tight loops, but can be afforded in case of heavy data analysis.
    /// You really should avoid it in case of trace handling
    std::vector<Hit> getHits()
    {
      static thread_local std::vector<Hit> hits; // creating a static vector allows creating it only once per thread, avoiding systematic re-allocation
      hits.clear();
      hits.reserve(mult);
      for (int hit_i = 0; hit_i<mult; ++hit_i) hits.push_back(getHit(hit_i));
      return hits;
    }

    constexpr auto dT(int index1, int index2) const noexcept {return Long64_t(time[index1]-time[index2]);}

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
    Trace     traces        [maxEvt];
  };
}