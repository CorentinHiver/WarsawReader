#pragma once

#include "utils.hpp"
#include "BoardAggregate.hpp"
#include <RtypesCore.h>

/**
  UInt_t    label         : Global label (=board_ID*16 + channel_ID*2 + subchannel_ID)
  UShort_t  board_ID      : Board label [0;max board ID]
  UShort_t  channel_ID    : Channel label [0;8]
  UShort_t  subchannel_ID : Sub channel label [0,1]
  Int_t     adc           : PHA : ADC value.               PSD : qshort value.
  Int_t     qlong         : PHA : EXTRA[16:32] (see doc.). PSD : qlong value.
  ULong64_t caen_time     : Raw caen_time. Units : ps. Is equal to extended_timestamp if fine or extended caen_time is found in the data, or TRIGGER_TIME_TAG if no extended caen_time found
  ULong64_t time          : Absolute time. Units : ps. This field must be filled by the user (i.e., the program using this class).
*/

namespace Caen1725
{
  template<typename T>
  using Trace_t = std::vector<T>;
  using Trace = Trace_t<uint16_t>;
  using DP1_t = Trace_t<bool>;
  /**
   * @brief Interface between the caen binary data and root TTree.
   * @note  There are different time and timstamps, please read carefully the comments.
   * @todo Check weither when reading from root tree, one needs to delete the trace in the destructor.
   */
  class Hit
  {
  protected: // Options :
    bool handle_traces = true; // Option to skip trace reading - much faster data processing, but no trace analysis available.

  public:
    // Data fields :
    UInt_t    label         = {}; // Global label (=board_ID*16 + channel_ID*2 + subchannel_ID)
    UShort_t  board_ID      = {}; // Board label [0;max board ID]
    UChar_t   channel_ID    = {}; // Channel label [0;8]
    UChar_t   subchannel_ID = {}; // Sub channel label [0,1]
    Int_t     adc           = {}; // PHA : ADC value.               PSD : qshort value.
    Int_t     qlong         = {}; // PHA : EXTRA[16:32] (see doc.). PSD : qlong value.
    ULong64_t caen_time     = {}; // Raw caen_time (ts). Units : ps. Is equal to extended_timestamp in fine and extended timestamps modes, or TRIGGER_TIME_TAG if no extended caen_time found
    ULong64_t time          = {}; // Absolute time. Units : ps. This field must be filled by the user (i.e., the program using this class).
    Int_t     rel_time      = {}; // Relative time in the event. Units : ps. This field must be filled by the user (i.e., the program using this class).
    Bool_t    wfa_success   = {};
    
    // Internal variables :
    ULong64_t extended_ts   = {}; // Raw caen_time. Units : ps. Used only if fine caen_time and extended caen_time are found in the data.
    ULong64_t precise_ts    = {}; // Raw caen_time. Units : ps. Used only if fine caen_time mode is found in the data.
    
    // Trace-related fields :
    Trace trace; // Signal trace.
    DP1_t DP1 ; // Digital probe. Look at documentation. Use printTrace() to visualise it.
    uint16_t trigger_bin = {};     // Trigger position in the trace. Use printTrace() to visualise it.

    // Root io parameters :
    
    Hit(bool _handle_traces = true) noexcept:
      handle_traces (_handle_traces)
    {}

    void handleTraces(bool _handle_traces = true) noexcept {handle_traces = _handle_traces;}
    
    ~Hit() {}

    using Board1725 = Caen1725::BoardAggregate;
    using Channel1725 = Caen1725::ChannelAggregate;
    using Event1725 = Caen1725::CaenEvent;

    /**
     * @brief Loads the buffer in the hit. You can use skipTrace to choose hit-by-hit if trace analysis is needed (usefull to skip trace analysis only for some detectors)
     */
    void readCaenEvent(std::istream& data, Board1725 & board, Channel1725 & channel, Event1725 & caenEvent, bool loadTrace)
    {
      auto const pos_init = data.tellg();

      board_ID     = static_cast<UChar_t>(board  .BOARD_ID);
      channel_ID   = static_cast<UChar_t>(channel.ID);
      caenEvent.EX = static_cast<UChar_t>(channel.EX);

      // 1. Timestamp and subchannel_ID

      Caen1725::read_data(data, tmp_u32);
      debug("CH[31], TRIGGER_TIME_TAG[30:0] :", std::bitset<32>(tmp_u32));
      caenEvent.TRIGGER_TIME_TAG         = getBitField(tmp_u32, 30);
      subchannel_ID = static_cast<UChar_t>(getBit     (tmp_u32, 31));
      
      // 2. Trace : looping through the NUM_SAMPLES samples of the waveform of each event :
      if (handle_traces && loadTrace)
      {
        auto const & N = channel.NUM_SAMPLES;
        debug("Trace with :", N, "samples");
        trace.clear();
        DP1  .clear();
        trace.resize(N);
        DP1  .resize(N);
        for (size_t sample_i = 0; sample_i<N; ++sample_i)
        {
          auto & sample = trace[sample_i];      // Simple aliasing
          Caen1725::read_data(data, sample);  // Reading the buffer
          if(getBit(sample, 15)) trigger_bin = sample_i; // Gets the index of the sample where the trigger time tag have been measured
          DP1[sample_i] = getBit(sample, 14);   // Gets the digital probe DP1 value for this sample
          sample &= mask(13);                            // Removes the trigger_bin and DP1 bit values from the sample 
        }
      }
      else 
      { // If traces are not used, skip the data
        trace.clear();
        DP1  .clear();
        auto const & size_to_skip = channel.NUM_SAMPLES * sizeof(uint16_t);
        Caen1725::skip(data, size_to_skip);
      }

      // 3. EXTRAS2 field (EXTRAS for PSD)

      if (channel.E2) 
      {
        Caen1725::read_data(data, &caenEvent.EXTRAS2);
        debug("EXTRAS2", std::bitset<32>(caenEvent.EXTRAS2));

        caenEvent.extra = Caen1725::Extra2(caenEvent.EXTRAS2, caenEvent.TRIGGER_TIME_TAG, caenEvent.EX);
        
        switch (caenEvent.extra.flag)
        {
          case Caen1725::Extra2::ExtendedTimestamp_flag : 
            extended_ts = caenEvent.extra.extended_timestamp;
            caen_time   = caenEvent.extra.extended_timestamp; 
            break;

          case Caen1725::Extra2::FineTimestamp_flag : 
            extended_ts = caenEvent.extra.extended_timestamp;
            precise_ts  = caenEvent.extra.precise_timestamp ;
            caen_time   = caenEvent.extra.precise_timestamp;
            break;

          default: 
            caen_time = caenEvent.TRIGGER_TIME_TAG * ticks_to_ps;
            break;
        }
      }

      // 4. Energy, PU (analog probe) and EXTRAS (qlong for PSD)

      Caen1725::read_data(data, &tmp_u32);

      debug("EXTRAS, PU, ENERGY:", std::bitset<32>(tmp_u32));

      adc   = getBitField (tmp_u32, 14    );
      qlong = getBitField (tmp_u32, 31, 16);

      // 5. Label

      label = static_cast<UShort_t>(board_ID * 16 + channel_ID * 2 + subchannel_ID);

      // 6. Size of the hit in the data stream

      auto const read_size = data.tellg() - pos_init;
      board  .read_size += read_size;
      channel.read_size += read_size;
    }
    
    std::vector<int> getTraceBaselineRemoved(size_t nb_samples_baseline) const 
    {
      std::vector<int> traces;
      
      if (!handle_traces) {error("Trace not handled"); return traces;}
      if (trace.empty()) {error("No trace"); return traces;}
      if (trace.size() == 0) {error("trace has no samples"); return traces;}

      std::vector<int> _trace; _trace.reserve(trace.size());
      int baseline = 0;
      for (size_t i = 0; i<nb_samples_baseline; ++i) baseline += trace.at(i);
      baseline /= nb_samples_baseline;
      for (auto const & sample : trace) _trace . push_back(sample - baseline);
      return _trace;
    }

    bool hasTrace() const {return !trace.empty();}

    void clear()
    {
      label         = {};
      board_ID      = {};
      channel_ID    = {};
      subchannel_ID = {};
      adc           = {};
      qlong         = {};
      caen_time     = {};
      extended_ts   = {};
      precise_ts    = {};
      time          = {};
      rel_time      = {};
      wfa_success   = {};
      trace.clear();
      DP1.clear();
    }

    friend std::ostream& operator<<(std::ostream& out, Hit const & hit)
    {
     out << 
       " label "           <<  hit.label      << " " <<
       " board_ID "        <<  hit.board_ID   << " " <<
       " channel_ID "      <<  hit.channel_ID << " " <<
      //  " label "           <<  Colib::fill(hit.label     , 3, '0')todo: recoder Colib::fill, il a du se perdre dans les versionings...
      //  " board_ID "        <<  Colib::fill(hit.board_ID  , 3, '0') <<
      //  " channel_ID "      <<  Colib::fill(hit.channel_ID, 3, '0') <<
       " subchannel_ID "   <<              hit.subchannel_ID       <<
        std::scientific     ;
      if (hit.caen_time   != 0) out << " caen_time "   <<  std::setprecision(10) << double_cast(hit.caen_time)  ;
      if (hit.extended_ts != 0) out << " extended_ts " <<  std::setprecision(10) << double_cast(hit.extended_ts);
      if (hit.precise_ts  != 0) out << " precise_ts "  <<  std::setprecision(10) << double_cast(hit.precise_ts) ;
      if (hit.time        != 0) out << " time "        <<  std::setprecision(10) << double_cast(hit.time)       ;
      if (hit.rel_time    != 0) out << " rel_time "    <<  std::setprecision(10) << double_cast(hit.rel_time)   ;
      if (hit.adc         != 0) out << " adc "         <<                                       hit.adc         ;
      if (hit.qlong       != 0) out << " qlong "       <<                                       hit.qlong       ;
      if (!hit.trace.empty()  ) out << " trace "       <<  hit.trace.size() << " samples "                      ;
      if (!hit.wfa_success    ) out << " wave form analysis success "                                           ;
      out << std::setprecision(6);
      return out;
    }
 
    Hit(Hit      && other) noexcept = default; 
    Hit(Hit const & other) = delete; // Non copyable
    Hit& operator=(Hit&& other) noexcept = default;

    Hit (
      UInt_t    _label,
      UShort_t  _board_ID,
      UShort_t  _channel_ID,
      UShort_t  _subchannel_ID,
      Int_t     _adc,
      Int_t     _qlong,
      ULong64_t _timestamp,
      ULong64_t _time,
      Int_t     _rel_time,
      Bool_t    _wfa_success
    ) noexcept : 
      label         (_label        ),
      board_ID      (_board_ID     ),
      channel_ID    (_channel_ID   ),
      subchannel_ID (_subchannel_ID),
      adc           (_adc          ),
      qlong         (_qlong        ),
      caen_time     (_timestamp    ),
      time          (_time         ),
      rel_time      (_rel_time     ),
      wfa_success   (_wfa_success  )
    {
    }

    Hit (
      UInt_t    _label,
      UShort_t  _board_ID,
      UShort_t  _channel_ID,
      UShort_t  _subchannel_ID,
      Int_t     _adc,
      Int_t     _qlong,
      ULong64_t _timestamp,
      ULong64_t _time,
      Int_t     _rel_time,
      Trace     && _trace,
      Bool_t    _wfa_success
    ) noexcept : 
      label         (_label        ),
      board_ID      (_board_ID     ),
      channel_ID    (_channel_ID   ),
      subchannel_ID (_subchannel_ID),
      adc           (_adc          ),
      qlong         (_qlong        ),
      caen_time     (_timestamp    ),
      time          (_time         ),
      rel_time      (_rel_time     ),
      wfa_success   (_wfa_success  ),
      trace         (std::move(_trace))
    {
      if (!_trace.empty())
      {
        handleTraces(true);
      }
    }
  };
};