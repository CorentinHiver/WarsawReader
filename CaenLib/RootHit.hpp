#ifndef ROOTHIT_HPP
#define ROOTHIT_HPP

#include "utils.hpp"
#include "BoardAggregate.hpp"

#include "TTree.h"
#include "TGraph.h"

namespace CaenDataReader
{
  using Trace = std::vector<uint16_t>;
  struct RootHit
  {
    bool handle_traces  = true;
    int  board_ID       = 0;
    int  channel_ID     = 0;
    int  subchannel_ID  = 0;
    int  adc            = 0;
    int  qlong          = 0;
    uint64_t timestamp = 0;
    uint64_t precise_ts = 0;
    uint64_t extended_ts = 0;
    uint64_t cfd       = 0;
    
    uint16_t trigger_bin = 0;

    Trace  outTrace;
    Trace* inTrace = nullptr;

    RootHit() noexcept = default;

    void writeTo(TTree * outTree)
    {
      outTree->Branch("board_ID"      , &board_ID     );
      outTree->Branch("channel_ID"    , &channel_ID   );
      outTree->Branch("subchannel_ID" , &subchannel_ID);
      outTree->Branch("trace"         , &outTrace     );
      outTree->Branch("timestamp"     , &timestamp    );
      outTree->Branch("adc"           , &adc          );
      outTree->Branch("qlong"         , &qlong        );
    }
    
    void readFrom(TTree * inTree)
    {
      inTree->SetBranchAddress("board_ID"      , &board_ID     );
      inTree->SetBranchAddress("channel_ID"    , &channel_ID   );
      inTree->SetBranchAddress("subchannel_ID" , &subchannel_ID);
      inTree->SetBranchAddress("trace"         , &inTrace      );
      inTree->SetBranchAddress("timestamp"     , &timestamp    );
      inTree->SetBranchAddress("adc"           , &adc          );
      inTree->SetBranchAddress("qlong"         , &qlong        );
    }

    bool operator>(RootHit const & other) const {return this->timestamp > other.timestamp;}

    void readCaenEvent(std::istream& data, BoardAggregate & board, ChannelAggregate & channel, CaenEvent & caenEvent)
    {
      auto const pos_init = data.tellg();

      board_ID = board.BOARD_ID;
      channel_ID = channel.ID;
      
      caenEvent.EX = channel.EX;

      read_buff(&tmp_u32, data);
      debug("TRIGGER_TIME_TAG, CH:", std::bitset<32>(tmp_u32));
      caenEvent.TRIGGER_TIME_TAG = getBitField(tmp_u32, 30);
      subchannel_ID    = int_cast(getBit     (tmp_u32, 31));
      
      // Looping through the NUM_SAMPLES samples of the waveform of each event :
      if (handle_traces) {
        debug("Trace with :", channel.NUM_SAMPLES, "samples");
        outTrace.resize(channel.NUM_SAMPLES);
        for (size_t sample_i = 0; sample_i<outTrace.size(); ++sample_i)
        {
          auto & sample = outTrace[sample_i];
          read_buff(&sample, data);
          if(getBit(sample, 14)) trigger_bin = sample_i; // Gets the index of the sample where the trigger time tag have been measured
          sample &= mask(13);
        }
      }
      else 
      {
        outTrace.clear();
        auto const & size_to_skip = channel.NUM_SAMPLES * sizeof(uint16_t);
        skip(data, size_to_skip);
      }

      if (channel.E2)  // TODO: do we need to skip EXTRAS2 if it is disabled, i.e. if channel.E2=false ?
        read_buff(&caenEvent.EXTRAS2, data);
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

      read_buff(&tmp_u32, data);

      debug("ENERGY, PU, EXTRAS:", std::bitset<32>(tmp_u32));

      adc   = getBitField (tmp_u32, 14    );
      qlong = getBitField (tmp_u32, 25, 16);

      auto const read_size = data.tellg() - pos_init;
      board  .read_size += read_size;
      channel.read_size += read_size;
    }

    TGraph* drawTrace()
    {
      std::vector<int> data_int; data_int.reserve(inTrace->size());
      for (auto const & data : *inTrace) data_int.push_back(data);
      auto graph = new TGraph(data_int.size(), linspace<int>(data_int.size()).data(), data_int.data());
      graph->Draw();
      return graph;
    }

    friend std::ostream& operator<<(std::ostream& out, RootHit const & hit)
    {
     out << 
      " board_ID "        <<  hit.board_ID          <<
      " channel_ID "      <<  hit.channel_ID        <<
      " subchannel_ID "   <<  hit.subchannel_ID     <<
      " timestamp "       <<  hit.timestamp         <<
      " adc "             <<  hit.adc               <<
      " qlong "           <<  hit.qlong             ;
      if (hit.inTrace)  out << " inTrace "  << hit.inTrace->size() << " samples ";
      else              out << " outTrace " << hit.outTrace.size() << " samples ";
      return out;
    }
  };
};

using RootCaenHit = CaenDataReader::RootHit;

#endif //ROOTHIT_HPP
