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
    Trace     traces        [maxSize];
        
    Event(bool handle_traces = false) : m_handle_traces(handle_traces)
    {
      // if (handle_traces) Colib::throw_error("Caen1725::Event can't handle traces (TBD)");
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
      if (!hit.trace.empty()) traces [i] = hit.trace;
      ++mult;
    }

    void push_back(Hit && hit) noexcept
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
      if (!hit.trace.empty()) traces [i] = std::move(hit.trace);
      ++mult;
    }

    // void emplace_back(Hit && hit) noexcept
    // {
    //   if (maxSize <= static_cast<size_t>(mult)) return;

    //   const auto i = mult; // In principle, this helps the compiler optimize
    //   label        [i] = hit.label         ;
    //   board_ID     [i] = hit.board_ID      ;
    //   channel_ID   [i] = hit.channel_ID    ;
    //   subchannel_ID[i] = hit.subchannel_ID ;
    //   adc          [i] = hit.adc           ;
    //   qlong        [i] = hit.qlong         ;
    //   caen_time    [i] = hit.caen_time     ;
    //   time         [i] = hit.time          ;
    //   rel_time     [i] = ((i == 0) ? 0 : static_cast<Int_t>(hit.time - time[0]));
    //   wfa_success  [i] = hit.wfa_success;
    //   // traces       .emplace_back(std::move(hit.trace));

    //   ++mult;

    //   hit.clear();
    // }

    void clear() 
    {
      // for (int hit_i = 0; hit_i<mult; ++hit_i) delete (traces[hit_i]);
      mult = 0;
    }
    
    auto size () const {return static_cast<size_t>(mult);}
    auto const & multiplicity () const {return mult;}

    Hit getHit(size_t const i) const
    {
      // if (traces[i]) return Hit(
      //   label         [i],
      //   board_ID      [i],
      //   channel_ID    [i],
      //   subchannel_ID [i],
      //   adc           [i],
      //   qlong         [i],
      //   caen_time     [i],
      //   time          [i],
      //   rel_time      [i],
      //   wfa_success   [i],
        // traces[i]
      // );
      // else 
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

    Hit operator[](size_t const i) const
    {
      return getHit(i);
    }

    /// @brief More user-friendly interface : copies the data inside a vector of hits.
    /// @details This method can be very expensive in tight loops, but can be afforded in case of heavy data analysis.
    /// You really should avoid it in case of trace handling
    // std::vector<Hit> getHits()
    // {
    //   static thread_local std::vector<Hit> hits; // creating a static vector allows creating it only once per thread, avoiding systematic re-allocation
    //   hits.clear();
    //   hits.reserve(mult);
    //   for (int hit_i = 0; hit_i<mult; ++hit_i) hits.emplace_back(getHit(hit_i));
    //   return hits;
    // }

    constexpr auto dT(int index1, int index2) const noexcept {return Long64_t(time[index1]-time[index2]);}
    
  protected:
    bool m_handle_traces = false;
  };
}