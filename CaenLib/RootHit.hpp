#pragma once

#include "Hit.hpp"
#include "RootUtils.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"

namespace Caen1725
{
  class RootHit : public Hit
  {
  public:
    
    /// @brief Inherit constructors from Hit
    template<class... ARGS>
    RootHit(ARGS &&... args) : Hit(std::forward<ARGS>(args)...)
    {}
    
    /// @brief Resets outTree and creates the TBranches with the Hit attributes. Sets RootHit in writting mode.
    inline TTree * writeTo(TTree * outTree) noexcept
    {
      ioMode = IOmode::Writting;
      outTree->ResetBranchAddresses();
      outTree->Branch("label"        , &label        );
      outTree->Branch("board_ID"     , &board_ID     );
      outTree->Branch("channel_ID"   , &channel_ID   );
      outTree->Branch("subchannel_ID", &subchannel_ID);
      outTree->Branch("timestamp"    , &timestamp    );
      outTree->Branch("time"         , &time         );
      outTree->Branch("rel_time"     , &rel_time     );
      outTree->Branch("adc"          , &adc          );
      outTree->Branch("qlong"        , &qlong        );
      if (Hit::handle_traces && trace) outTree->Branch("trace", &trace);
      
      return outTree;
    }

    /// @brief Resets inTree and set the TBranches addresses to the Hit attributes. Sets RootHit in reading mode.
    TTree * readFrom(TTree * inTree)
    {
      ioMode = IOmode::Reading;
      inTree->ResetBranchAddresses();
      inTree->SetBranchAddress("label"         , &label        );
      inTree->SetBranchAddress("board_ID"      , &board_ID     );
      inTree->SetBranchAddress("channel_ID"    , &channel_ID   );
      inTree->SetBranchAddress("subchannel_ID" , &subchannel_ID);
      inTree->SetBranchAddress("timestamp"     , &timestamp    );
      inTree->SetBranchAddress("time"          , &time         );
      inTree->SetBranchAddress("rel_time"      , &rel_time     );
      inTree->SetBranchAddress("adc"           , &adc          );
      inTree->SetBranchAddress("qlong"         , &qlong        );
      if (Hit::handle_traces) inTree->SetBranchAddress("trace", &trace);
      
      return inTree;
    }

    TTree * readFrom(TFile * file, std::string const & treename = "HIL")
    {
      auto tree = file->Get<TTree>(treename.c_str());
      if (!tree) error("Can't extract a valid", treename, "in file", file->GetName());
      else readFrom(tree);
      return tree;
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

      graphs.push_back(new TGraph(data_int.size(), Colib::linspace<int>(data_int.size(), 0, Caen1725::ticks_to_ns).data(), data_int.data()));
      graphs.push_back(new TGraph(DP1_int .size(), Colib::linspace<int>(DP1_int .size(), 0, Caen1725::ticks_to_ns).data(), DP1_int .data()));
      graphs.push_back(new TGraph(trig_int.size(), Colib::linspace<int>(trig_int.size(), 0, Caen1725::ticks_to_ns).data(), trig_int.data()));

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
      auto graph = new TGraph(data_int.size(), Colib::linspace<int>(data_int.size(), 0, Caen1725::ticks_to_ns).data(), data_int.data());
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

    enum class IOmode {None, Reading, Writting} ioMode = IOmode::None; // Read or Write a TTree
  };
};

using Caen1725RootHit = Caen1725::RootHit;