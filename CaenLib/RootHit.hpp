#ifndef ROOTHIT_HPP
#define ROOTHIT_HPP

#include "utils.hpp"
#include "RootUtils.hpp"
#include "BoardAggregate.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"

namespace CaenDataReader
{
  template<typename T>
  using Trace_t = std::vector<T>;
  using Trace = std::vector<uint16_t>;

  /**
   * @brief Interface between the caen binary data and root TTree.
   * @note  There are different time and timstamps, please read carefully the comments.
   * @attention When not using it in a root environnemnet, one must call RootHit::clean() to correctly delete the pointer to the trace
   * @todo Check weither when reading from root tree, one needs to delete the trace in the destructor.
   */
  class RootHit
  {
  private: // Options :
    bool handle_traces = true; // Option to skip trace reading - much faster data processing, but no trace analysis available.

  public:
    // Data fields :
    int      label         = 0; // Global label (=board_ID*16 + channel_ID*2 + subchannel_ID)
    u_short  board_ID      = 0; // Board label [0;max board ID]
    u_short  channel_ID    = 0; // Channel label [0;8]
    u_short  subchannel_ID = 0; // Sub channel label [0,1]
    int      adc           = 0; // PHA : ADC value. PSD : qshort value.
    int      qlong         = 0; // PHA : unused.    PSD : qlong value.
    uint64_t timestamp     = 0; // Raw timestamp (ts). Units : ps. Is equal to precise_ts if fine ts is found in the data, or extended_ts if only extended ts is found, or TRIGGER_TIME_TAG if no extended timestamp found
    uint64_t extended_ts   = 0; // Raw timestamp. Units : ps. Used only if fine timestamp and extended timestamp are found in the data.
    uint64_t precise_ts    = 0; // Raw timestamp. Units : ps. Used only if fine timestamp mode is found in the data.
    uint64_t time          = 0; // Absolute time. Units : ps. This field must be filled by the user (i.e., the program using this class).

    // Trace-related fields :
    Trace*         trace = nullptr; // Signal trace.
    Trace_t<bool>* DP1   = nullptr; // Digital probe. Look at documentation. Use printTrace() to visualise it.
    uint16_t trigger_bin = 0;       // Trigger position in the trace. Use printTrace() to visualise it.

    // Root io parameters :
    
    bool reading  = false; // Reading  to   a TTree
    bool writting = false; // Writting from a TTree

    void handleTraces(bool const & _handle_traces)
    {
      handle_traces = _handle_traces;
      if (_handle_traces)
      {
        if (!trace) trace = new Trace;
        if (!DP1  ) DP1   = new Trace_t<bool>;  
      }
      else
      {
        if (trace && trace -> size() > 0) delete trace;
        if (DP1   && DP1   -> size() > 0) delete DP1  ;
      }
    }
    
    RootHit(bool handle_traces = true) noexcept
    {
      if (handle_traces)
      {
        trace = new Trace;
        DP1   = new Trace_t<bool>;  
      }
    }
    
    ~RootHit() {delete trace; delete DP1;}

    TTree * writeTo(TTree * outTree)
    {
      writting = true;
      outTree->ResetBranchAddresses();
      outTree->Branch("label"        , &label        );
      outTree->Branch("board_ID"     , &board_ID     );
      outTree->Branch("channel_ID"   , &channel_ID   );
      outTree->Branch("subchannel_ID", &subchannel_ID);
      outTree->Branch("timestamp"    , &timestamp    );
      outTree->Branch("time"         , &time         );
      outTree->Branch("adc"          , &adc          );
      outTree->Branch("qlong"        , &qlong        );
      if (handle_traces) outTree->Branch("trace", &trace);
      // outTree->Branch("extended_ts"  , &extended_ts  );
      // outTree->Branch("precise_ts"   , &precise_ts   );
      
      return outTree;
    }

    TTree * readFrom(TTree * inTree)
    {
      reading = true;
      inTree->ResetBranchAddresses();
      inTree->SetBranchAddress("label"         , &label        );
      inTree->SetBranchAddress("board_ID"      , &board_ID     );
      inTree->SetBranchAddress("channel_ID"    , &channel_ID   );
      inTree->SetBranchAddress("subchannel_ID" , &subchannel_ID);
      inTree->SetBranchAddress("timestamp"     , &timestamp    );
      inTree->SetBranchAddress("time"          , &time         );
      inTree->SetBranchAddress("adc"           , &adc          );
      inTree->SetBranchAddress("qlong"         , &qlong        );
      // Special treatment for the trace : If not found in the data, then deactivate trace handling
      if (handle_traces) handleTraces(inTree->SetBranchAddress("trace", &trace) == 0);
      // inTree->SetBranchAddress("extended_ts"   , &extended_ts  );
      // inTree->SetBranchAddress("precise_ts"    , &precise_ts   );
      
      return inTree;
    }

    TTree * readFrom(TFile * file, std::string const & treename = "HIL")
    {
      auto tree = file->Get<TTree>(treename.c_str());
      if (!tree) error("Can't extract a valid", treename, "in file", file->GetName());
      else this->readFrom(tree);
      return tree;
    }
   
    void readCaenEvent(std::istream& data, CaenDataReader1725::BoardAggregate & board, CaenDataReader1725::ChannelAggregate & channel, CaenDataReader1725::CaenEvent & caenEvent)
    {
      auto const pos_init = data.tellg();

      board_ID     = board  .BOARD_ID;
      channel_ID   = channel.ID;
      caenEvent.EX = channel.EX;

      // 1. Timestamp and subchannel_ID

      CaenDataReader1725::read_buff(&tmp_u32, data);
      debug("TRIGGER_TIME_TAG, CH:", std::bitset<32>(tmp_u32));
      caenEvent.TRIGGER_TIME_TAG = getBitField(tmp_u32, 30);
      subchannel_ID     = int_cast(getBit     (tmp_u32, 31));
      
      // 2. Trace : looping through the NUM_SAMPLES samples of the waveform of each event :
      if (handle_traces)
      {
        debug("Trace with :", channel.NUM_SAMPLES, "samples");
        trace->resize(channel.NUM_SAMPLES);
        DP1  ->resize(channel.NUM_SAMPLES);
        for (size_t sample_i = 0; sample_i<trace->size(); ++sample_i)
        {
          auto & sample = (*trace)[sample_i];             // Simple aliasing
          CaenDataReader1725::read_buff(&sample, data);   // Reading the buffer
          if(getBit(sample, 15)) trigger_bin = sample_i;  // Gets the index of the sample where the trigger time tag have been measured
          (*DP1)[sample_i] = getBit(sample, 14);          // Gets the digital probe DP1 value for this sample
          sample &= mask(13);                             // Removes the trigger_bin and DP1 bit values from the sample 
        }
      }
      else 
      { // If traces are not used, skip the data
        trace->clear();
        DP1  ->clear();
        auto const & size_to_skip = channel.NUM_SAMPLES * sizeof(uint16_t);
        CaenDataReader1725::skip(data, size_to_skip);
      }

      // 3. EXTRAS2 field

      if (channel.E2) CaenDataReader1725::read_buff(&caenEvent.EXTRAS2, data);
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
          timestamp   = caenEvent.extra.precise_timestamp ;
          break;

        default: 
          timestamp = caenEvent.TRIGGER_TIME_TAG;
          break;
      }

      // 4. Energy, PU (analog probe) and EXTRAS

      CaenDataReader1725::read_buff(&tmp_u32, data);

      debug("EXTRAS, PU, ENERGY:", std::bitset<32>(tmp_u32));

      adc   = getBitField (tmp_u32, 14    );
      qlong = getBitField (tmp_u32, 25, 16);

      // 5. Label

      label = board_ID * 16 + channel_ID * 2 + subchannel_ID;

      // 6. Size of the hit in the data

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

    /// @brief Returns three graphs : ret = {trace, DP1, Trigger}
    std::vector<TGraph*> getTracesGraphs(size_t nb_samples_baseline = 0) const
    {
      std::vector<TGraph*> graphs;
      if (!handle_traces) {error("Trace not handled"); return graphs;}
      if (!trace) {error("No trace"); return graphs;}
      if (trace->size() == 0) {error("trace has no samples"); return graphs;}

      double baseline = 0;
      if (nb_samples_baseline != 0)
      {
        if (trace->size() < nb_samples_baseline) {error("trace has not enough samples for baseline(",nb_samples_baseline," samples required)"); return graphs;}
        for (size_t sample_i = 0; sample_i<nb_samples_baseline; ++sample_i) baseline += trace->at(sample_i);
        baseline /= nb_samples_baseline;
      }

      auto const & N = trace->size();
      auto const & maxS = Colib::maximum(*trace);
      std::vector<int> data_int(N, 0);
      std::vector<int> DP1_int (N, 0);
      std::vector<int> trig_int(N, 0);

      for (size_t i = 0; i<N; ++i) 
      {
        data_int[i] = int((*trace)[i]) - baseline;
        if (DP1_int[i]) DP1_int[i] = maxS;
        if (i == trigger_bin) trig_int[i] = maxS;
      }

      graphs.push_back(new TGraph(data_int.size(), Colib::linspace<int>(data_int.size(), 0, CaenDataReader1725::ticks_to_ns).data(), data_int.data()));
      graphs.push_back(new TGraph(DP1_int .size(), Colib::linspace<int>(DP1_int .size(), 0, CaenDataReader1725::ticks_to_ns).data(), DP1_int .data()));
      graphs.push_back(new TGraph(trig_int.size(), Colib::linspace<int>(trig_int.size(), 0, CaenDataReader1725::ticks_to_ns).data(), trig_int.data()));

      graphs[0] -> SetName("Trace"  ); graphs[0] -> SetTitle("Trace"  );
      graphs[1] -> SetName("DP1"    ); graphs[1] -> SetTitle("DP1"    );
      graphs[1] -> SetName("Trigger"); graphs[1] -> SetTitle("Trigger");

      return graphs;
    }

    /// @brief Returns the graph of the trace
    TGraph* getTraceGraph(size_t nb_samples_baseline = 0) const
    {
      if (!handle_traces) {error("Trace not handled"); return nullptr;}
      if (!trace) {error("No trace"); return nullptr;}
      if (trace->size() == 0) {error("trace has no samples"); return nullptr;}

      auto const & N = trace->size();
      double baseline = 0;
      if (0 < nb_samples_baseline)
      {
        if (trace->size() < nb_samples_baseline) {error("trace has not enough samples for baseline(",nb_samples_baseline," samples required)"); return nullptr;}
        for (size_t sample_i = 0; sample_i<nb_samples_baseline; ++sample_i) baseline += trace->at(sample_i) - baseline;
        baseline /= nb_samples_baseline;
      }
      std::vector<int> data_int; data_int.reserve(N);
      for (auto const sample : *trace) data_int.push_back(sample);
      auto graph = new TGraph(data_int.size(), Colib::linspace<int>(data_int.size(), 0, CaenDataReader1725::ticks_to_ns).data(), data_int.data());
      graph -> SetName ("Trace"); 
      graph -> SetTitle("Trace");
      return graph;
    }

    /// @brief Draws and returns three graphs : ret = {trace, DP1, Trigger}
    std::vector<TGraph*> drawTraces(std::string options = "") const
    {
      std::vector<TGraph*> graphs;
      if (!handle_traces) {error("Trace not handled"); return graphs;}
      if (!trace) {error("No trace"); return graphs;}
      if (trace->size() == 0) {error("trace has no samples"); return graphs;}

      graphs = getTracesGraphs();
      graphs[0] -> Draw(options.c_str());
      graphs[1] -> Draw((options+"same").c_str());
      graphs[2] -> Draw((options+"same").c_str());
      return graphs;
    }

    /// @brief Draws and returns the graph of the trace
    TGraph* drawTrace(std::string options = "") const
    {
      TGraph* graph = nullptr;
      if (!handle_traces) {error("Trace not handled"); return graph;}
      if (!trace) {error("No trace"); return graph;}
      if (trace->size() == 0) {error("trace has no samples"); return graph;}

      graph = getTraceGraph();
      graph -> Draw(options.c_str());
      return graph;
    }

    /// @brief Get the trace.
    auto const & getTrace() const {return trace;}

    friend std::ostream& operator<<(std::ostream& out, RootHit const & hit)
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
      if (hit.adc         != 0) out << " adc "         <<                                       hit.adc         ;
      if (hit.qlong       != 0) out << " qlong "       <<                                       hit.qlong       ;
      if (hit.trace           ) out << " trace "       <<  hit.trace -> size() << " samples "                   ;
      out << std::setprecision(6);
      return out;
    }
 
    bool operator>(RootHit const & other) const {return this->timestamp > other.timestamp;}

    RootHit (RootHit const & other) noexcept : 
      label         (other.label        ),
      board_ID      (other.board_ID     ),
      channel_ID    (other.channel_ID   ),
      subchannel_ID (other.subchannel_ID),
      adc           (other.adc          ),
      qlong         (other.qlong        ),
      timestamp     (other.timestamp    ),
      extended_ts   (other.extended_ts  ),
      precise_ts    (other.precise_ts   ),
      time          (other.time         )
    {
      if (other.trace)
      {
        trace         = new Trace;
        *trace        = *other.trace;
      }
    }

    RootHit const & copy(RootHit const & other, bool copyTrace = true)
    {
      label         = other.label        ;
      board_ID      = other.board_ID     ;
      channel_ID    = other.channel_ID   ;
      subchannel_ID = other.subchannel_ID;
      adc           = other.adc          ;
      qlong         = other.qlong        ;
      timestamp     = other.timestamp    ;
      extended_ts   = other.extended_ts  ;
      precise_ts    = other.precise_ts   ;
      time          = other.time         ;
      if (copyTrace && other.trace)
      {
        if (!trace) trace = new Trace;
        *trace = *other.trace;
      }
      return other;
    }

    RootHit const & operator=(RootHit const & other)
    {
      return this -> copy(other);      
    }    
  };

  class RootEvent 
  {
    // Helper functions :

    /// @brief Create a branch for a given array and name
    /// @param name_size: The name of the leaf that holds the size of the array
    template<class T>
    auto createBranchArray(TTree* tree, std::string const & name, T * array, std::string const & name_size, int buffsize = 64000)
    {
      auto const & type_root_format = name+"["+name_size+"]/"+typeRoot(**array);
      return (tree -> Branch(name.c_str(), array, type_root_format.c_str(), buffsize));
    }

  public:
    RootEvent(bool handle_traces = false)
    {
      if (handle_traces) Colib::throw_error("CaenDataReader1725::RootEvent can't handle traces (TBD)");
      // hits = new std::vector<RootHit>;
      // label          = new std::vector<int     >;
      // board_ID       = new std::vector<int     >;
      // channel_ID     = new std::vector<int     >;
      // subchannel_ID  = new std::vector<int     >;
      // adc            = new std::vector<int     >;
      // qlong          = new std::vector<int     >;
      // timestamp      = new std::vector<uint64_t>;
      // extended_ts    = new std::vector<uint64_t>;
      // precise_ts     = new std::vector<uint64_t>;
      // time           = new std::vector<uint64_t>;
    }

    ~RootEvent()
    {
      // delete label        ;
      // delete board_ID     ;
      // delete channel_ID   ;
      // delete subchannel_ID;
      // delete adc          ;
      // delete qlong        ;
      // delete timestamp    ;
      // delete extended_ts  ;
      // delete precise_ts   ;
      // delete time         ;
    }

    // void emplace_back(RootHit && hit)
    // {
    //   hits -> emplace_back(std::move(hit));
    //   ++mult;
    // }

    void push_back(RootHit const & hit)
    {
      // hits -> push_back(hit);

      // label         -> push_back(hit.label        );
      // board_ID      -> push_back(hit.board_ID     );
      // channel_ID    -> push_back(hit.channel_ID   );
      // subchannel_ID -> push_back(hit.subchannel_ID);
      // adc           -> push_back(hit.adc          );
      // qlong         -> push_back(hit.qlong        );
      // timestamp     -> push_back(hit.timestamp    );
      // extended_ts   -> push_back(hit.extended_ts  );
      // precise_ts    -> push_back(hit.precise_ts   );
      // time          -> push_back(hit.time         );

      label        [mult] = hit.label        ;
      board_ID     [mult] = hit.board_ID     ;
      channel_ID   [mult] = hit.channel_ID   ;
      subchannel_ID[mult] = hit.subchannel_ID;
      adc          [mult] = hit.adc          ;
      qlong        [mult] = hit.qlong        ;
      timestamp    [mult] = hit.timestamp    ;
      extended_ts  [mult] = hit.extended_ts  ;
      precise_ts   [mult] = hit.precise_ts   ;
      time         [mult] = hit.time         ;

      ++mult;
    }

    void clear()
    {
      // hits -> clear();

      // label         -> clear();
      // board_ID      -> clear();
      // channel_ID    -> clear();
      // subchannel_ID -> clear();
      // adc           -> clear();
      // qlong         -> clear();
      // timestamp     -> clear();
      // extended_ts   -> clear();
      // precise_ts    -> clear();
      // time          -> clear();

      mult = 0;
    }

    void writeTo(TTree* tree)
    {
      tree->ResetBranchAddresses();
      tree->Branch("evtNb"        , &evtNb        );
      tree->Branch("mult"         , &mult         );

      createBranchArray(tree, "label"        , &label        , "mult");
      createBranchArray(tree, "board_ID"     , &board_ID     , "mult");
      createBranchArray(tree, "channel_ID"   , &channel_ID   , "mult");
      createBranchArray(tree, "subchannel_ID", &subchannel_ID, "mult");
      createBranchArray(tree, "adc"          , &adc          , "mult");
      createBranchArray(tree, "qlong"        , &qlong        , "mult");
      createBranchArray(tree, "timestamp"    , &timestamp    , "mult");
      createBranchArray(tree, "time"         , &time         , "mult");

      // tree->Branch("hits"         , &hits         );

      // tree->Branch("label"        , &label        );
      // tree->Branch("board_ID"     , &board_ID     );
      // tree->Branch("channel_ID"   , &channel_ID   );
      // tree->Branch("subchannel_ID", &subchannel_ID);
      // tree->Branch("adc"          , &adc          );
      // tree->Branch("qlong"        , &qlong        );
      // tree->Branch("timestamp"    , &timestamp    );
      // tree->Branch("time"         , &time         );

      // tree->Branch("extended_ts", &extended_ts);
      // tree->Branch("precise_ts", &precise_ts);
    }

    TTree * readFrom(TTree* tree)
    {
      tree->ResetBranchAddresses();
      tree->SetBranchAddress("evtNb"        , &evtNb        );
      tree->SetBranchAddress("mult"         , &mult         );
      // tree->SetBranchAddress("hits"         , &hits         );
      tree->SetBranchAddress("label"        , &label        );
      tree->SetBranchAddress("board_ID"     , &board_ID     );
      tree->SetBranchAddress("channel_ID"   , &channel_ID   );
      tree->SetBranchAddress("subchannel_ID", &subchannel_ID);
      tree->SetBranchAddress("adc"          , &adc          );
      tree->SetBranchAddress("qlong"        , &qlong        );
      tree->SetBranchAddress("timestamp"    , &timestamp    );
      tree->SetBranchAddress("time"         , &time         );
      // tree->SetBranchAddress("extended_ts", &extended_ts);
      // tree->SetBranchAddress("precise_ts", &precise_ts);
      return tree;
    }

    TTree * readFrom(TFile * file, std::string const & treename = "HIL")
    {
      auto tree = file->Get<TTree>(treename.c_str());
      if (!tree) error("Can't extract a valid", treename, "in file", file->GetName());
      else this->readFrom(tree);
      return tree;
    }

    auto const & size () const {return mult;}

    RootHit operator[](size_t const & i) const
    {
      RootHit hit;

      // hit.label         = (*label         ) [i];
      // hit.board_ID      = (*board_ID      ) [i];
      // hit.channel_ID    = (*channel_ID    ) [i];
      // hit.subchannel_ID = (*subchannel_ID ) [i];
      // hit.adc           = (*adc           ) [i];
      // hit.qlong         = (*qlong         ) [i];
      // hit.timestamp     = (*timestamp     ) [i];
      // hit.time          = (*time          ) [i];

      hit.label         = label         [i];
      hit.board_ID      = board_ID      [i];
      hit.channel_ID    = channel_ID    [i];
      hit.subchannel_ID = subchannel_ID [i];
      hit.adc           = adc           [i];
      hit.qlong         = qlong         [i];
      hit.timestamp     = timestamp     [i];
      hit.time          = time          [i];

      // hit.extended_ts   = (*extended_ts)   [i];
      // hit.precise_ts    = (*precise_ts)    [i];

      // return (*hits)[i];
      return hit;
    }

    friend std::ostream& operator<<(std::ostream& out, RootEvent const & event)
    {
      out << "Event with " << event.mult << " hits :" << std::endl;
      for (int hit_i = 0; hit_i<event.mult; ++hit_i)
      {
        out << "\t";
        
        out << " label "         << event.label        [hit_i];
        out << " board_ID "      << event.board_ID     [hit_i];
        out << " channel_ID "    << event.channel_ID   [hit_i];
        out << " subchannel_ID " << event.subchannel_ID[hit_i];
        out << " adc "           << event.adc          [hit_i];
        out << " qlong "         << event.qlong        [hit_i];
        out << " timestamp "     << event.timestamp    [hit_i];
        out << " time "          << event.time         [hit_i];

        // out << " label "         << (*event.label        )[hit_i];
        // out << " board_ID "      << (*event.board_ID     )[hit_i];
        // out << " channel_ID "    << (*event.channel_ID   )[hit_i];
        // out << " subchannel_ID " << (*event.subchannel_ID)[hit_i];
        // out << " adc "           << (*event.adc          )[hit_i];
        // out << " qlong "         << (*event.qlong        )[hit_i];
        // out << " timestamp "     << (*event.timestamp    )[hit_i];
        // out << " time "          << (*event.time         )[hit_i];

        // out << " label "         << (*event.hits)[hit_i].label        ;
        // out << " board_ID "      << (*event.hits)[hit_i].board_ID     ;
        // out << " channel_ID "    << (*event.hits)[hit_i].channel_ID   ;
        // out << " subchannel_ID " << (*event.hits)[hit_i].subchannel_ID;
        // out << " adc "           << (*event.hits)[hit_i].adc          ;
        // out << " qlong "         << (*event.hits)[hit_i].qlong        ;
        // out << " timestamp "     << (*event.hits)[hit_i].timestamp    ;
        // out << " time "          << (*event.hits)[hit_i].time         ;
        // out << "extended_ts " << (*event.extended_ts)[hit_i];
        // out << "precise_ts " << (*event.precise_ts)[hit_i];
        out << std::endl;
      }
      return out;
    }

    size_t evtNb = 0;
    int    mult  = 0;
    // std::vector<RootHit> * hits = nullptr;
    constexpr static inline size_t maxEvt = 1000;
    int      label         [maxEvt];
    int      board_ID      [maxEvt];
    int      channel_ID    [maxEvt];
    int      subchannel_ID [maxEvt];
    int      adc           [maxEvt];
    int      qlong         [maxEvt];
    uint64_t timestamp     [maxEvt];
    uint64_t extended_ts   [maxEvt];
    uint64_t precise_ts    [maxEvt];
    uint64_t time          [maxEvt];
    // std::vector<int>      * label         = nullptr;
    // std::vector<int>      * board_ID      = nullptr;
    // std::vector<int>      * channel_ID    = nullptr;
    // std::vector<int>      * subchannel_ID = nullptr;
    // std::vector<int>      * adc           = nullptr;
    // std::vector<int>      * qlong         = nullptr;
    // std::vector<uint64_t> * timestamp     = nullptr;
    // std::vector<uint64_t> * extended_ts   = nullptr;
    // std::vector<uint64_t> * precise_ts    = nullptr;
    // std::vector<uint64_t> * time          = nullptr;
  };
};

using RootCaenHit   = CaenDataReader1725::RootHit  ;
using RootCaenEvent = CaenDataReader1725::RootEvent;

#endif //ROOTHIT_HPP
