#ifndef ROOTHIT_HPP
#define ROOTHIT_HPP

#include "utils.hpp"
#include "BoardAggregate.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"

namespace CaenDataReader
{
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
    uint64_t timestamp   = 0;
    uint64_t precise_ts  = 0;
    uint64_t extended_ts = 0;
    uint64_t cfd         = 0;
    Trace* trace = nullptr;
    
    // Options :
    uint16_t trigger_bin = 0;
    bool handle_traces   = true;
    
      
    RootHit() noexcept : trace (new Trace) {}
    ~RootHit() {delete trace;}

    TTree * writeTo(TTree * outTree)
    {
      outTree->Branch("label"         , &label        );
      outTree->Branch("board_ID"      , &board_ID     );
      outTree->Branch("channel_ID"    , &channel_ID   );
      outTree->Branch("subchannel_ID" , &subchannel_ID);
      outTree->Branch("timestamp"     , &timestamp    );
      outTree->Branch("precise_ts"    , &precise_ts   );
      outTree->Branch("extended_ts"   , &extended_ts  );
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
      inTree->SetBranchAddress("precise_ts"    , &precise_ts   );
      inTree->SetBranchAddress("extended_ts"   , &extended_ts  );
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
   
    void readCaenEvent(std::istream& data, BoardAggregate & board, ChannelAggregate & channel, CaenEvent & caenEvent)
    {
      auto const pos_init = data.tellg();

      board_ID = board.BOARD_ID;
      channel_ID = channel.ID;
      
      caenEvent.EX = channel.EX;

      // 1. Timestamp and subchannel_ID

      CaenDataReader::read_buff(&tmp_u32, data);
      debug("TRIGGER_TIME_TAG, CH:", std::bitset<32>(tmp_u32));
      caenEvent.TRIGGER_TIME_TAG = getBitField(tmp_u32, 30);
      subchannel_ID     = int_cast(getBit     (tmp_u32, 31));
      
      // 2. Trace : looping through the NUM_SAMPLES samples of the waveform of each event :
      if (handle_traces) {
        debug("Trace with :", channel.NUM_SAMPLES, "samples");
        if (m_NewTrace)
        {
          trace = new Trace;
        }
        trace->resize(channel.NUM_SAMPLES);
        for (size_t sample_i = 0; sample_i<trace->size(); ++sample_i)
        {
          auto & sample = (*trace)[sample_i];
          CaenDataReader::read_buff(&sample, data);
          if(getBit(sample, 14)) trigger_bin = sample_i; // Gets the index of the sample where the trigger time tag have been measured
          sample &= mask(13);
        }
      }
      else 
      { // If traces are not used, skip the data
        trace->clear();
        auto const & size_to_skip = channel.NUM_SAMPLES * sizeof(uint16_t);
        CaenDataReader::skip(data, size_to_skip);
      }

      // 3. EXTRAS2 field

      if (channel.E2)  // TODO: do we need to skip EXTRAS2 if it is disabled, i.e. if channel.E2=false ?
        CaenDataReader::read_buff(&caenEvent.EXTRAS2, data);
      debug("EXTRAS2", std::bitset<32>(caenEvent.EXTRAS2));

      caenEvent.extra = Extra2(caenEvent.EXTRAS2, caenEvent.TRIGGER_TIME_TAG, caenEvent.EX);
      
      switch (caenEvent.extra.flag)
      {
        case Extra2::ExtendedTimestamp_flag : 
          extended_ts = caenEvent.extra.extended_timestamp;
          timestamp = caenEvent.extra.extended_timestamp; 
          break;

        case Extra2::FineTimestamp_flag     : 
          precise_ts = caenEvent.extra.precise_timestamp;
          extended_ts = caenEvent.extra.extended_timestamp;
          timestamp = caenEvent.extra.precise_timestamp ; 
          break;

        default: 
          timestamp = caenEvent.TRIGGER_TIME_TAG; 
          break;
      }

      // 4. Energy, PU (analog probe) and EXTRAS

      CaenDataReader::read_buff(&tmp_u32, data);

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

    TGraph* drawTrace() const
    {
      if (trace->size() == 0) {error("trace has no point"); return nullptr;}
      std::vector<int> data_int; data_int.reserve(trace->size());
      for (auto const & data : *trace) data_int.push_back(data);
      auto graph = new TGraph(data_int.size(), Colib::linspace<int>(data_int.size()).data(), data_int.data());
      graph -> SetName ("Hit");
      graph -> SetTitle("Hit");
      graph -> Draw();
      return graph;
    }

    auto const & getTrace() const {return *trace;}

    friend std::ostream& operator<<(std::ostream& out, RootHit const & hit)
    {
     out << 
       " label "           <<  Colib::fill(hit.label     , 3, '0') <<
       " board_ID "        <<  Colib::fill(hit.board_ID  , 3, '0') <<
       " channel_ID "      <<  Colib::fill(hit.channel_ID, 3, '0') <<
       " subchannel_ID "   <<              hit.subchannel_ID       <<
        std::scientific     ;
      if (hit.timestamp   != 0) out << " timestamp "   <<  std::setprecision(10) << double_cast(hit.timestamp)  ;
      if (hit.precise_ts  != 0) out << " precise_ts "  <<  std::setprecision(10) << double_cast(hit.precise_ts) ;
      if (hit.extended_ts != 0) out << " extended_ts " <<  std::setprecision(10) << double_cast(hit.extended_ts);
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
      precise_ts    (other.precise_ts),
      extended_ts   (other.extended_ts),
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
      precise_ts    = other.precise_ts;
      extended_ts   = other.extended_ts;
      cfd           = other.cfd;
      *trace        = *other.trace;
      
      return other;
    }

    RootHit const & operator=(RootHit const & other)
    {
      return this->copy(other);      
    }
    
    /// @brief This is for when the 
    void createNewTrace(bool const & b = true) {m_NewTrace = b;}

  private:
    bool m_NewTrace = false;
  };
};

using RootCaenHit = CaenDataReader::RootHit;

#endif //ROOTHIT_HPP
