#pragma once
#include "Event.hpp"
#include "RootHit.hpp"

namespace Caen1725
{
  class RootEvent : public Event
  {
  public:
    /// @brief Inherit constructors from Hit
    template<class... ARGS>
    RootEvent(ARGS &&... args) : Event(std::forward<ARGS>(args)...)
    {}

    virtual ~RootEvent()
    {
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
      createBranchArray(tree, "caen_time"    , &caen_time    , "mult");
      createBranchArray(tree, "time"         , &time         , "mult");
      createBranchArray(tree, "rel_time"     , &rel_time     , "mult");
      createBranchArray(tree, "wfa_success"  , &wfa_success  , "mult");
    }

    TTree * readFrom(TTree* tree)
    {
      tree->ResetBranchAddresses();
      tree->SetBranchAddress("evtNb"        , &evtNb        );
      tree->SetBranchAddress("mult"         , &mult         );
      tree->SetBranchAddress("label"        , &label        );
      tree->SetBranchAddress("board_ID"     , &board_ID     );
      tree->SetBranchAddress("channel_ID"   , &channel_ID   );
      tree->SetBranchAddress("subchannel_ID", &subchannel_ID);
      tree->SetBranchAddress("adc"          , &adc          );
      tree->SetBranchAddress("qlong"        , &qlong        );
      tree->SetBranchAddress("caen_time"    , &caen_time    );
      tree->SetBranchAddress("time"         , &time         );
      tree->SetBranchAddress("rel_time"     , &rel_time     );
      tree->SetBranchAddress("wfa_success"  , &wfa_success  );
      return tree;
    }

    TTree * readFrom(TFile * file, std::string const & treename = "HIL")
    {
      auto tree = file->Get<TTree>(treename.c_str());
      if (!tree) error("Can't extract a valid", treename, "in file", file->GetName());
      else this->readFrom(tree);
      return tree;
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
        out << " caen_time "     << event.caen_time    [hit_i];
        out << " time "          << event.time         [hit_i];
        out << " rel_time "      << event.rel_time     [hit_i];
        out << " wfa_success "   << nicer_bool(event.wfa_success[hit_i]);

        out << std::endl;
      }
      return out;
    }
  };

};