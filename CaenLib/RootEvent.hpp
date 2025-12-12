#include "RootHit.hpp"

namespace CaenDataReader
{
  class RootEvent 
  {
    // Helper functions :

  public:
    RootEvent(bool handle_traces = false)
    {
      if (handle_traces) Colib::throw_error("CaenDataReader1725::RootEvent can't handle traces (TBD)");
    }

    ~RootEvent()
    {
    }

    void push_back(RootHit const & hit)
    {
      if (maxEvt < static_cast<size_t>(mult)) error("Event too big > 1000");

      label        [mult] = hit.label         ;
      board_ID     [mult] = hit.board_ID      ;
      channel_ID   [mult] = hit.channel_ID    ;
      subchannel_ID[mult] = hit.subchannel_ID ;
      adc          [mult] = hit.adc           ;
      qlong        [mult] = hit.qlong         ;
      timestamp    [mult] = hit.timestamp     ;
      time         [mult] = hit.time          ;
      rel_time     [mult] = static_cast<Int_t>(hit.time - time[0]);

      ++mult;
    }

    void clear()
    {
      for (int hit_i = 0; hit_i<mult; ++hit_i)
      { // might be useless for most fields !! Mandatory for at least qlong
        label        [hit_i] = 0;
        board_ID     [hit_i] = 0;
        channel_ID   [hit_i] = 0;
        subchannel_ID[hit_i] = 0;
        adc          [hit_i] = 0;
        qlong        [hit_i] = 0;
        timestamp    [hit_i] = 0;
        time         [hit_i] = 0;
        rel_time     [hit_i] = 0;
      }
      mult = 0;
    }

    void writeTo(TTree* tree)
    {
      tree->ResetBranchAddresses();
      tree->Branch     (      "evtNb"        , &evtNb                );
      tree->Branch     (      "mult"         , &mult                 );
      createBranchArray(tree, "label"        , &label        , "mult");
      createBranchArray(tree, "board_ID"     , &board_ID     , "mult");
      createBranchArray(tree, "channel_ID"   , &channel_ID   , "mult");
      createBranchArray(tree, "subchannel_ID", &subchannel_ID, "mult");
      createBranchArray(tree, "adc"          , &adc          , "mult");
      createBranchArray(tree, "qlong"        , &qlong        , "mult");
      createBranchArray(tree, "timestamp"    , &timestamp    , "mult");
      createBranchArray(tree, "time"         , &time         , "mult");
      createBranchArray(tree, "rel_time"     , &rel_time     , "mult");
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
      tree->SetBranchAddress("timestamp"    , &timestamp    );
      tree->SetBranchAddress("time"         , &time         );
      tree->SetBranchAddress("rel_time"     , &rel_time     );
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

      hit.label         = label         [i];
      hit.board_ID      = board_ID      [i];
      hit.channel_ID    = channel_ID    [i];
      hit.subchannel_ID = subchannel_ID [i];
      hit.adc           = adc           [i];
      hit.qlong         = qlong         [i];
      hit.timestamp     = timestamp     [i];
      hit.time          = time          [i];
      hit.rel_time      = rel_time      [i];

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
        out << " rel_time "      << event.rel_time     [hit_i];

        out << std::endl;
      }
      return out;
    }

    constexpr static inline size_t maxEvt = 1000;

    size_t evtNb = 0;
    int    mult  = 0;
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

  class RootEventVec 
  {
    // Helper functions :

  public:
    RootEventVec(bool handle_traces = false)
    {
      if (handle_traces) Colib::throw_error("CaenDataReader1725::RootEventVec can't handle traces (TBD)");
      label          = new std::vector<Int_t    >;
      board_ID       = new std::vector<Int_t    >;
      channel_ID     = new std::vector<Int_t    >;
      subchannel_ID  = new std::vector<Int_t    >;
      adc            = new std::vector<Int_t    >;
      qlong          = new std::vector<Int_t    >;
      timestamp      = new std::vector<ULong64_t>;
      time           = new std::vector<ULong64_t>;
      rel_time       = new std::vector<ULong64_t>;
    }

    ~RootEventVec()
    {
      delete label        ;
      delete board_ID     ;
      delete channel_ID   ;
      delete subchannel_ID;
      delete adc          ;
      delete qlong        ;
      delete timestamp    ;
      delete time         ;
      delete rel_time     ;
    }

    void push_back(RootHit const & hit)
    {
      if (maxEvt < static_cast<size_t>(mult)) error("Event too big > 1000");

      label         -> push_back(hit.label        );
      board_ID      -> push_back(hit.board_ID     );
      channel_ID    -> push_back(hit.channel_ID   );
      subchannel_ID -> push_back(hit.subchannel_ID);
      adc           -> push_back(hit.adc          );
      qlong         -> push_back(hit.qlong        );
      timestamp     -> push_back(hit.timestamp    );
      time          -> push_back(hit.time         );
      rel_time      -> push_back(hit.rel_time     );

      ++mult;
    }

    void clear()
    {
      label         -> clear();
      board_ID      -> clear();
      channel_ID    -> clear();
      subchannel_ID -> clear();
      adc           -> clear();
      qlong         -> clear();
      timestamp     -> clear();
      time          -> clear();
      rel_time      -> clear();
      mult = 0;
    }

    void writeTo(TTree* tree)
    {
      tree->ResetBranchAddresses();
      tree->Branch("evtNb"        , &evtNb        );
      tree->Branch("mult"         , &mult         );      
      tree->Branch("label"        , &label        );
      tree->Branch("board_ID"     , &board_ID     );
      tree->Branch("channel_ID"   , &channel_ID   );
      tree->Branch("subchannel_ID", &subchannel_ID);
      tree->Branch("adc"          , &adc          );
      tree->Branch("qlong"        , &qlong        );
      tree->Branch("timestamp"    , &timestamp    );
      tree->Branch("time"         , &time         );
      tree->Branch("rel_time"     , &rel_time     );
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
      tree->SetBranchAddress("timestamp"    , &timestamp    );
      tree->SetBranchAddress("time"         , &time         );
      tree->SetBranchAddress("rel_time"     , &rel_time     );
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

      hit.label         = (*label         ) [i];
      hit.board_ID      = (*board_ID      ) [i];
      hit.channel_ID    = (*channel_ID    ) [i];
      hit.subchannel_ID = (*subchannel_ID ) [i];
      hit.adc           = (*adc           ) [i];
      hit.qlong         = (*qlong         ) [i];
      hit.timestamp     = (*timestamp     ) [i];
      hit.time          = (*time          ) [i];
      hit.rel_time      = (*rel_time      ) [i];

      return hit;
    }

    friend std::ostream& operator<<(std::ostream& out, RootEventVec const & event)
    {
      out << "Event with " << event.mult << " hits :" << std::endl;
      for (int hit_i = 0; hit_i<event.mult; ++hit_i)
      {
        out << "\t";
        
        out << " label "         << (*event.label        )[hit_i];
        out << " board_ID "      << (*event.board_ID     )[hit_i];
        out << " channel_ID "    << (*event.channel_ID   )[hit_i];
        out << " subchannel_ID " << (*event.subchannel_ID)[hit_i];
        out << " adc "           << (*event.adc          )[hit_i];
        out << " qlong "         << (*event.qlong        )[hit_i];
        out << " timestamp "     << (*event.timestamp    )[hit_i];
        out << " time "          << (*event.time         )[hit_i];
        out << " rel_time "      << (*event.rel_time     )[hit_i];

        out << std::endl;
      }
      return out;
    }

    size_t evtNb = 0;
    int    mult  = 0;
    // std::vector<RootHit> * hits = nullptr;
    constexpr static inline size_t maxEvt = 1000;
    std::vector<Int_t>      * label         = nullptr;
    std::vector<Int_t>      * board_ID      = nullptr;
    std::vector<Int_t>      * channel_ID    = nullptr;
    std::vector<Int_t>      * subchannel_ID = nullptr;
    std::vector<Int_t>      * adc           = nullptr;
    std::vector<Int_t>      * qlong         = nullptr;
    std::vector<ULong64_t>  * timestamp     = nullptr;
    std::vector<ULong64_t>  * time          = nullptr;
    std::vector<ULong64_t>  * rel_time      = nullptr;
  };
};
using RootCaenEvent = CaenDataReader1725::RootEvent;