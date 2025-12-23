#pragma once

#include "utils.hpp"
#include "BoardAggregate.hpp"
#include <RtypesCore.h>

namespace CaenDataReader1725
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
    UShort_t  label         = 0; // Global label (=board_ID*16 + channel_ID*2 + subchannel_ID)
    UChar_t   board_ID      = 0; // Board label [0;max board ID]
    UChar_t   channel_ID    = 0; // Channel label [0;8]
    UChar_t   subchannel_ID = 0; // Sub channel label [0,1]
    UShort_t  adc           = 0; // PHA : ADC value.               PSD : qshort value.
    UShort_t  qlong         = 0; // PHA : EXTRA[16:32] (see doc.). PSD : qlong value.
    ULong64_t timestamp     = 0; // Raw timestamp (ts). Units : ps. Is equal to precise_ts if fine ts is found in the data, or extended_ts if only extended ts is found, or TRIGGER_TIME_TAG if no extended timestamp found
    ULong64_t extended_ts   = 0; // Raw timestamp. Units : ps. Used only if fine timestamp and extended timestamp are found in the data.
    ULong64_t precise_ts    = 0; // Raw timestamp. Units : ps. Used only if fine timestamp mode is found in the data.
    ULong64_t time          = 0; // Absolute time. Units : ps. This field must be filled by the user (i.e., the program using this class).
    Int_t     rel_time      = 0; // Relative time in the event. Units : ps. This field must be filled by the user (i.e., the program using this class).
    
    // Trace-related fields :
    std::unique_ptr<Trace> trace; // Signal trace.
    std::unique_ptr<DP1_t>  DP1 ; // Digital probe. Look at documentation. Use printTrace() to visualise it.
    uint16_t trigger_bin = 0;     // Trigger position in the trace. Use printTrace() to visualise it.

    // Root io parameters :
    
    Hit(bool _handle_traces = true) noexcept:
      handle_traces (_handle_traces)
    {}

    void handleTraces(bool _handle_traces = true) noexcept {handle_traces = _handle_traces;}
    
    ~Hit() {}

    bool hasTrace() const {return ( static_cast<bool>(trace) && (!trace -> empty()) );}

    using Board1725 = CaenDataReader1725::BoardAggregate;
    using Channel1725 = CaenDataReader1725::ChannelAggregate;
    using Event1725 = CaenDataReader1725::CaenEvent;

    void readCaenEvent(std::istream& data, Board1725 & board, Channel1725 & channel, Event1725 & caenEvent, std::vector<bool> const & boardReadTrace = {})
    {
      auto const pos_init = data.tellg();

      board_ID     = static_cast<UChar_t>(board  .BOARD_ID);
      channel_ID   = static_cast<UChar_t>(channel.ID);
      caenEvent.EX = static_cast<UChar_t>(channel.EX);

      // 1. Timestamp and subchannel_ID

      CaenDataReader1725::read_data(data, &tmp_u32);
      debug("CH[31], TRIGGER_TIME_TAG[30:0] :", std::bitset<32>(tmp_u32));
      caenEvent.TRIGGER_TIME_TAG         = getBitField(tmp_u32, 30);
      subchannel_ID = static_cast<UChar_t>(getBit     (tmp_u32, 31));
      
      // 2. Trace : looping through the NUM_SAMPLES samples of the waveform of each event :
      if (handle_traces)
      {
        if((board_ID < boardReadTrace.size()) ? true : boardReadTrace[board_ID])
        {
          debug("Trace with :", channel.NUM_SAMPLES, "samples");
          if (!trace) trace.reset(new Trace);
          if (!DP1  ) DP1  .reset(new DP1_t);
          trace -> resize(channel.NUM_SAMPLES);
          DP1   -> resize(channel.NUM_SAMPLES);
          for (size_t sample_i = 0; sample_i<trace->size(); ++sample_i)
          {
            auto & sample = (*trace.get())[sample_i];      // Simple aliasing
            CaenDataReader1725::read_data(data, &sample);  // Reading the buffer
            if(getBit(sample, 15)) trigger_bin = sample_i; // Gets the index of the sample where the trigger time tag have been measured
            (*DP1.get())[sample_i] = getBit(sample, 14);   // Gets the digital probe DP1 value for this sample
            sample &= mask(13);                            // Removes the trigger_bin and DP1 bit values from the sample 
          }
        }
      }
      else 
      { // If traces are not used, skip the data
        if (trace) trace -> clear();
        if (DP1  ) DP1   -> clear();
        auto const & size_to_skip = channel.NUM_SAMPLES * sizeof(uint16_t);
        CaenDataReader1725::skip(data, size_to_skip);
      }

      // 3. EXTRAS2 field (EXTRAS for PSD)

      if (channel.E2) CaenDataReader1725::read_data(data, &caenEvent.EXTRAS2);
      debug("EXTRAS2", std::bitset<32>(caenEvent.EXTRAS2));

      caenEvent.extra = CaenDataReader1725::Extra2(caenEvent.EXTRAS2, caenEvent.TRIGGER_TIME_TAG, caenEvent.EX);
      
      switch (caenEvent.extra.flag)
      {
        case CaenDataReader1725::Extra2::ExtendedTimestamp_flag : 
          extended_ts = caenEvent.extra.extended_timestamp;
          timestamp   = caenEvent.extra.extended_timestamp; 
          break;

        case CaenDataReader1725::Extra2::FineTimestamp_flag : 
          extended_ts = caenEvent.extra.extended_timestamp;
          precise_ts  = caenEvent.extra.precise_timestamp ;
          timestamp   = caenEvent.extra.extended_timestamp;
          break;

        default: 
          timestamp = caenEvent.TRIGGER_TIME_TAG;
          break;
      }

      // 4. Energy, PU (analog probe) and EXTRAS (qlong for PSD)

      CaenDataReader1725::read_data(data, &tmp_u32);

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
      if (!trace) {error("No trace"); return traces;}
      if (trace->size() == 0) {error("trace has no samples"); return traces;}

      std::vector<int> _trace; _trace.reserve(trace->size());
      int baseline = 0;
      for (size_t i = 0; i<nb_samples_baseline; ++i) baseline += trace->at(i);
      baseline /= nb_samples_baseline;
      for (auto const & sample : *trace) _trace . push_back(sample - baseline);
      return _trace;
    }

    /// @brief Get the trace.
    auto const & getTrace() const {return trace;}

    void clear()
    {
      label         = 0;
      board_ID      = 0;
      channel_ID    = 0;
      subchannel_ID = 0;
      adc           = 0;
      qlong         = 0;
      timestamp     = 0;
      extended_ts   = 0;
      precise_ts    = 0;
      time          = 0;
      rel_time      = 0;
      if (trace) trace->clear();
      if (DP1  ) trace->clear();
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
      if (hit.timestamp   != 0) out << " timestamp "   <<  std::setprecision(10) << double_cast(hit.timestamp)  ;
      if (hit.extended_ts != 0) out << " extended_ts " <<  std::setprecision(10) << double_cast(hit.extended_ts);
      if (hit.precise_ts  != 0) out << " precise_ts "  <<  std::setprecision(10) << double_cast(hit.precise_ts) ;
      if (hit.time        != 0) out << " time "        <<  std::setprecision(10) << double_cast(hit.time)       ;
      if (hit.rel_time    != 0) out << " rel_time "    <<  std::setprecision(10) << double_cast(hit.rel_time)   ;
      if (hit.adc         != 0) out << " adc "         <<                                       hit.adc         ;
      if (hit.qlong       != 0) out << " qlong "       <<                                       hit.qlong       ;
      if (hit.trace           ) out << " trace "       <<  hit.trace -> size() << " samples "                   ;
      out << std::setprecision(6);
      return out;
    }
 
    Hit(Hit      && other) noexcept = default; 
    Hit(Hit const & other) = delete; // Non copyable
    Hit& operator=(Hit&& other) noexcept = default;

    /// @brief Efficient copy of the trace. If trace is nullptr this function makes nothing.
    void traceEfficientCopy(Trace* _trace)
    {
      if (!_trace) return;
      trace->reserve(_trace->size());
      trace->clear();
      trace->resize(_trace->size());
      std::copy(_trace->begin(), _trace->end(), trace->begin());
    }

    // Hit (
    //   UInt_t    _label,
    //   UShort_t  _board_ID,
    //   UShort_t  _channel_ID,
    //   UShort_t  _subchannel_ID,
    //   Int_t     _adc,
    //   Int_t     _qlong,
    //   ULong64_t _timestamp,
    //   ULong64_t _time,
    //   Int_t     _rel_time,
    //   Trace     && _trace
    // ) noexcept : 
    //   label         (_label        ),
    //   board_ID      (_board_ID     ),
    //   channel_ID    (_channel_ID   ),
    //   subchannel_ID (_subchannel_ID),
    //   adc           (_adc          ),
    //   qlong         (_qlong        ),
    //   timestamp     (_timestamp    ),
    //   time          (_time         ),
    //   rel_time      (_rel_time     ),
    //   trace(std::make_unique<std::vector<uint16_t>>(std::move(_trace)))
    // {
    //   if (_trace.)
    //   {
    //     handleTraces(true);
    //   }
    // }

    Hit (
      UInt_t    _label,
      UShort_t  _board_ID,
      UShort_t  _channel_ID,
      UShort_t  _subchannel_ID,
      Int_t     _adc,
      Int_t     _qlong,
      ULong64_t _timestamp,
      ULong64_t _time,
      Int_t     _rel_time
    ) noexcept : 
      label         (_label        ),
      board_ID      (_board_ID     ),
      channel_ID    (_channel_ID   ),
      subchannel_ID (_subchannel_ID),
      adc           (_adc          ),
      qlong         (_qlong        ),
      timestamp     (_timestamp    ),
      time          (_time         ),
      rel_time      (_rel_time     )
    {
    }

    // Hit (
    //   UInt_t    _label,
    //   UShort_t  _board_ID,
    //   UShort_t  _channel_ID,
    //   UShort_t  _subchannel_ID,
    //   Int_t     _adc,
    //   Int_t     _qlong,
    //   ULong64_t _timestamp,
    //   ULong64_t _extended_ts,
    //   ULong64_t _precise_ts,
    //   ULong64_t _time,
    //   Int_t     _rel_time,
    //   Trace  && _trace
    // ) noexcept : 
    //   label         (_label        ),
    //   board_ID      (_board_ID     ),
    //   channel_ID    (_channel_ID   ),
    //   subchannel_ID (_subchannel_ID),
    //   adc           (_adc          ),
    //   qlong         (_qlong        ),
    //   timestamp     (_timestamp    ),
    //   extended_ts   (_extended_ts  ),
    //   precise_ts    (_precise_ts   ),
    //   time          (_time         ),
    //   rel_time      (_rel_time     ),
    //   trace(std::make_unique<std::vector<uint16_t>>(std::move(_trace)))
    // {
    //   // if (_trace)
    //   // {
    //     // handleTraces(true);
    //     // traceEfficientCopy(_trace);
    //   // }
    // }

    // Hit const & copy(Hit const & other, bool copyTrace = true)
    // {
    //   label         = other.label        ;
    //   board_ID      = other.board_ID     ;
    //   channel_ID    = other.channel_ID   ;
    //   subchannel_ID = other.subchannel_ID;
    //   adc           = other.adc          ;
    //   qlong         = other.qlong        ;
    //   timestamp     = other.timestamp    ;
    //   extended_ts   = other.extended_ts  ;
    //   precise_ts    = other.precise_ts   ;
    //   time          = other.time         ;
    //   rel_time      = other.rel_time     ;
    //   if (copyTrace && handle_traces && other.trace)
    //   {
    //     // if (!trace) trace.reset();
    //     // traceEfficientCopy(other.trace);
    //   }
    //   return other;
    // }
  };
};

using Caen1725Hit = CaenDataReader1725::Hit;