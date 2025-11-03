#ifndef ROOTHIT_HPP
#define ROOTHIT_HPP

#include "utils.hpp"
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
   * @brief 
   * @attention When not using it in a root environnemnet, one must call RootHit::clean() to correctly delete the pointer to the trace
   * @todo Check wether when reading from root tree, one needs to delete the trace in the destructor.
   */
  struct RootHit
  {
    // Fields :
    int  label           = 0;
    int  board_ID        = 0;
    int  channel_ID      = 0;
    int  subchannel_ID   = 0;
    int  adc             = 0;
    int  qlong           = 0;
    uint64_t timestamp   = 0; // TRIGGER_TIME_TAG [32 bits]
    uint64_t extended_ts = 0; // TRIGGER_TIME_TAG + extended timestamp [47 bits]
    uint64_t precise_ts  = 0; // TRIGGER_TIME_TAG + extended timestamp + fine timestamp[47 bits]
    uint64_t cfd         = 0;
    Trace* trace = nullptr;
    
    // Options :
    uint16_t trigger_bin = 0;
    Trace_t<bool>* DP1 = nullptr;
    bool handle_traces   = true;
    
    RootHit() noexcept : trace (new Trace), DP1 (new Trace_t<bool>) {}
    
    ~RootHit() {delete trace; delete DP1;}

    TTree * writeTo(TTree * outTree)
    {
      outTree->Branch("label"         , &label        );
      outTree->Branch("board_ID"      , &board_ID     );
      outTree->Branch("channel_ID"    , &channel_ID   );
      outTree->Branch("subchannel_ID" , &subchannel_ID);
      outTree->Branch("timestamp"     , &timestamp    );
      outTree->Branch("extended_ts"   , &extended_ts  );
      outTree->Branch("precise_ts"    , &precise_ts   );
      outTree->Branch("cfd"           , &cfd          );
      outTree->Branch("adc"           , &adc          );
      outTree->Branch("qlong"         , &qlong        );
      outTree->Branch("trace"         , &trace        );
      
      return outTree;
    }

    TTree * readFrom(TTree * inTree)
    {
      inTree->SetBranchAddress("label"         , &label        );
      inTree->SetBranchAddress("board_ID"      , &board_ID     );
      inTree->SetBranchAddress("channel_ID"    , &channel_ID   );
      inTree->SetBranchAddress("subchannel_ID" , &subchannel_ID);
      inTree->SetBranchAddress("timestamp"     , &timestamp    );
      inTree->SetBranchAddress("extended_ts"   , &extended_ts  );
      inTree->SetBranchAddress("precise_ts"    , &precise_ts   );
      inTree->SetBranchAddress("cfd"           , &cfd          );
      inTree->SetBranchAddress("adc"           , &adc          );
      inTree->SetBranchAddress("qlong"         , &qlong        );
      inTree->SetBranchAddress("trace"         , &trace        );
      
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
        // if (m_NewTrace)
        // {
          // trace = new Trace;
        // }
        trace->resize(channel.NUM_SAMPLES);
        DP1  ->resize(channel.NUM_SAMPLES);
        for (size_t sample_i = 0; sample_i<trace->size(); ++sample_i)
        {
          auto & sample = (*trace)[sample_i];             // Simple aliasing
          CaenDataReader1725::read_buff(&sample, data);   // Reading the buffer
          if(getBit(sample, 15)) trigger_bin = sample_i;  // Gets the index of the sample where the trigger time tag have been measured
          (*DP1)[sample_i] = getBit(sample, 14);        // Gets the digital probe DP1 value for this sample
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

      if (channel.E2)  // TODO: do we need to skip EXTRAS2 if it is disabled, i.e. if channel.E2=false ?
        CaenDataReader1725::read_buff(&caenEvent.EXTRAS2, data);
      debug("EXTRAS2", std::bitset<32>(caenEvent.EXTRAS2));

      caenEvent.extra = CaenDataReader1725::Extra2(caenEvent.EXTRAS2, caenEvent.TRIGGER_TIME_TAG, caenEvent.EX);
      
      switch (caenEvent.extra.flag)
      {
        case CaenDataReader1725::Extra2::ExtendedTimestamp_flag : 
          extended_ts = caenEvent.extra.extended_timestamp;
          timestamp = caenEvent.extra.extended_timestamp; 
          break;

        case CaenDataReader1725::Extra2::FineTimestamp_flag     : 
          precise_ts = caenEvent.extra.precise_timestamp;
          extended_ts = caenEvent.extra.extended_timestamp;
          timestamp = caenEvent.extra.precise_timestamp ; 
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

    /// @brief Returns three graphs : ret = {trace, DP1, Trigger}
    std::vector<TGraph*> getTracesGraphs() const
    {
      std::vector<TGraph*> graphs;
      if (!trace) {return graphs;}
      if (trace->size() == 0) {error("trace has no samples"); return graphs;}

      auto const & N = trace->size();
      auto const & maxS = Colib::maximum(*trace);
      std::vector<int> data_int(N, 0);
      std::vector<int> DP1_int (N, 0);
      std::vector<int> trig_int(N, 0);

      for (size_t i = 0; i<N; ++i)
      {
        data_int[i] = int((*trace)[i]);
        if (DP1_int[i]) DP1_int[i] = maxS;
        if (i == trigger_bin) trig_int[i] = maxS;
      }

      graphs.push_back(new TGraph(data_int.size(), Colib::linspace<int>(data_int.size()).data(), data_int.data()));
      graphs.push_back(new TGraph(DP1_int .size(), Colib::linspace<int>(DP1_int .size()).data(), DP1_int .data()));
      graphs.push_back(new TGraph(trig_int.size(), Colib::linspace<int>(trig_int.size()).data(), trig_int.data()));

      graphs[0] -> SetName("Trace"  ); graphs[0] -> SetTitle("Trace"  );
      graphs[1] -> SetName("DP1"    ); graphs[1] -> SetTitle("DP1"    );
      graphs[1] -> SetName("Trigger"); graphs[1] -> SetTitle("Trigger");

      return graphs;
    }

    /// @brief Returns the graph of the trace
    TGraph* getTraceGraph() const
    {
      if (!trace) {return nullptr;}
      auto const & N = trace->size();
      if (N == 0) {error("trace has no samples"); return nullptr;}
      std::vector<int> data_int; data_int.reserve(N);
      for (auto const sample : data_int) data_int.push_back(sample);
      auto graph = new TGraph(data_int.size(), Colib::linspace<int>(data_int.size()).data(), data_int.data());
      graph -> SetName ("Trace"); 
      graph -> SetTitle("Trace");
      return graph;
    }

    /// @brief Draws and returns three graphs : ret = {trace, DP1, Trigger}
    std::vector<TGraph*> drawTraces(std::string options = "") const
    {
      auto graphs = getTracesGraphs();
      graphs[0] -> Draw(options.c_str());
      graphs[1] -> Draw((options+"same").c_str());
      graphs[2] -> Draw((options+"same").c_str());
      return graphs;
    }

    /// @brief Draws and returns the graph of the trace
    TGraph* drawTrace(std::string options = "") const
    {
      auto graph = getTraceGraph();
      graph -> Draw(options.c_str());
      return graph;
    }

    auto const & getTrace() const {return *trace;}

    friend std::ostream& operator<<(std::ostream& out, RootHit const & hit)
    {
     out << 
       " label "           <<  hit.label      << " " <<
       " board_ID "        <<  hit.board_ID   << " " <<
       " channel_ID "      <<  hit.channel_ID << " " <<
      //  " label "           <<  Colib::fill(hit.label     , 3, '0') << // TODO: recoder Colib::fill, il a du se perdre dans les versionings...
      //  " board_ID "        <<  Colib::fill(hit.board_ID  , 3, '0') <<
      //  " channel_ID "      <<  Colib::fill(hit.channel_ID, 3, '0') <<
       " subchannel_ID "   <<              hit.subchannel_ID       <<
        std::scientific     ;
      if (hit.timestamp   != 0) out << " timestamp "   <<  std::setprecision(10) << double_cast(hit.timestamp)  ;
      if (hit.extended_ts != 0) out << " extended_ts " <<  std::setprecision(10) << double_cast(hit.extended_ts);
      if (hit.precise_ts  != 0) out << " precise_ts "  <<  std::setprecision(10) << double_cast(hit.precise_ts) ;
      if (hit.cfd         != 0) out << " cfd "         <<  std::setprecision(10) << double_cast(hit.cfd)        ;
      if (hit.adc         != 0) out << " adc "         <<                                       hit.adc         ;
      if (hit.qlong       != 0) out << " qlong "       <<                                       hit.qlong       ;
      if (hit.trace       != 0) out << " trace "       <<  hit.trace -> size() << " samples "                   ;
      out << std::setprecision(6);
      return out;
    }
 
    bool operator>(RootHit const & other) const {return this->timestamp > other.timestamp;}

    RootHit (RootHit const & other) noexcept : 
      label         (other.label),
      board_ID      (other.board_ID),
      channel_ID    (other.channel_ID),
      subchannel_ID (other.subchannel_ID),
      adc           (other.adc),
      qlong         (other.qlong),
      timestamp     (other.timestamp),
      extended_ts   (other.extended_ts),
      precise_ts    (other.precise_ts),
      cfd           (other.cfd),
      trace         (new Trace)
    {
      *trace        = *other.trace;
    }

    RootHit const & copy(RootHit const & other)
    {
      label         = other.label;
      board_ID      = other.board_ID;
      channel_ID    = other.channel_ID;
      subchannel_ID = other.subchannel_ID;
      adc           = other.adc;
      qlong         = other.qlong;
      timestamp     = other.timestamp;
      extended_ts   = other.extended_ts;
      precise_ts    = other.precise_ts;
      cfd           = other.cfd;
      *trace        = *other.trace;
      
      return other;
    }

    RootHit const & operator=(RootHit const & other)
    {
      return this->copy(other);      
    }
    
  private:
    // bool m_NewTrace = false;
  };

  using RootEvent = std::vector<RootHit>;
};

using RootCaenHit = CaenDataReader1725::RootHit;
using RootCaenEvent = CaenDataReader1725::RootEvent;

#endif //ROOTHIT_HPP
