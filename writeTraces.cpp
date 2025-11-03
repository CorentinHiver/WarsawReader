// The point of this code is to write down the traces of the signal
// The additionnal point is to study the cfd and compare with the original crrc2 timestamp,
// but without any coincidences, simply comparison of both methods !

// #include "CaenLib/CaenRawReader.hpp"
#include "CaenLib/CaenRootReader.hpp"
#include "AnalysisLib/CFD.hpp"
#include "AnalysisLib/RawHit.hpp"
#include "LibCo/Timer.hpp"
#include "LibCo/libCo.hpp"

#include "TPad.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph.h"

using namespace std;

int writeTraces(string file, int nb_events_max = -1, int adcMin = 0, int adcMax = -1)
{
  Timer timer;
  bool max_events = (nb_events_max>0);
  bool finished = false;
  int nb_evts = 0;
  unordered_map<int, size_t> shifts = {
    {EAGLE     , 5},
    {DSSDSector, 2},
    {DSSDRing  , 2}
  };
  unordered_map<int, int> thresholds = {
    {EAGLE     , -50},
    {DSSDSector, -500},
    {DSSDRing  , -500}
  };

  // CaenRawReader725 reader(file);
  CaenRootReader1725 reader(file);
  reader.handleTraces(true);

  auto rootFile = TFile::Open("writeTraces.root", "recreate"); rootFile->cd();

  vector<string> folders_names;

  while(!finished && reader.readHit())
  // while(!finished && reader.readBoardAggregate())
  {
    // auto const & board = reader.getBoard();
    // for (auto const & channel : board.channels) for (auto const & event : channel.events){
      ++nb_evts;
      auto const & hit = reader.getHit();
      if (nb_evts % int(1e3) == 0) print(Colib::nicer_double(nb_evts, 1));
      if (max_events && nb_evts > nb_events_max) finished = true;

      // RawHit hit(board, channel, event);
      // print(DetectorName[hit.detectorType], int(hit.label), int(hit.glabel), hit.adc, hit.precise_timestamp);

      // if (event.samples.size() == 0) continue;
      if (!hit.trace || hit.trace->size() == 0) continue;

      if (hit.adc < adcMin) continue;

      if (0 < adcMax && adcMax < hit.adc) continue;

      // auto detectorName = DetectorName[getDetectorType(board, event)]+to_string(labels(board, channel, event));
      auto detectorName = getDetectorName(hit.board_ID, hit.channel_ID)+to_string(hit.label);

      TDirectory *dir = rootFile->GetDirectory(detectorName.c_str());

      if (!dir){
        rootFile->mkdir(detectorName.c_str());
        dir = rootFile->GetDirectory(detectorName.c_str());
      }
      dir->cd();
      string name = detectorName + "_" + to_string(nb_evts);
      // vector<uint16_t> raw_trace; raw_trace.reserve(event.samples.size());
      // vector<double> samplesT; samplesT.resize(event.samples.size(), 0);
      // vector<double> samplesDP1; samplesDP1.resize(event.samples.size(), 0);
      // for (size_t sample_i = 0; sample_i<event.samples.size(); ++sample_i) 
      // {
      //   auto const & sample = event.samples[sample_i];
      //   raw_trace.push_back(sample.sample);
      //   if (sample.T) samplesT[sample_i] = 1;
      //   if (sample.DP1) samplesDP1[sample_i] = 1;
      // }

      // auto const & maxSample = maximum(*hit.trace);
      // for (auto & Ti : samplesT) if (Ti > 0.5) Ti = maxSample;
      // for (auto & DP1 : samplesDP1) if (DP1 > 0.5) DP1 = maxSample;
      
      // CFD cfd(raw_trace);

      // auto graph  = new TGraph(cfd.trace.size(), Colib::linspace_for(cfd.trace, 0., 4.).data(), cfd.trace.data());
      // auto graph2 = new TGraph(samplesT .size(), Colib::linspace_for(samplesT , 0., 4.).data(), samplesT .data());
      // auto graph3 = new TGraph(samplesDP1 .size(), Colib::linspace_for(samplesDP1 , 0., 4.).data(), samplesDP1 .data());
      
      auto canvas = std::make_unique<TCanvas>(name.c_str(), name.c_str()); canvas->cd();
      // graph->Draw();
      // graph2->SetLineColor(kGreen);
      // graph->GetXaxis()->SetTitle("time [ns]");
      // graph->GetYaxis()->SetTitle("pulse height [ADC]");
      // graph2->Draw("same");
      // graph3->SetLineColor(kBlue);
      // graph3->Draw("same");

      auto graphs = hit.getTracesGraphs();

      graphs[0] -> SetLineColor(kBlack);
      graphs[0] -> GetXaxis() -> SetTitle("time [ns]");
      graphs[0] -> GetYaxis() -> SetTitle("pulse height [ADC]");
      graphs[1] -> SetLineColor(kBlue );
      graphs[2] -> SetLineColor(kGreen);
      graphs[0] -> Draw();
      graphs[1] -> Draw("same");
      graphs[2] -> Draw("same");

      CFD cfd(*hit.trace);

      auto const & detectorType = getDetectorType(hit.board_ID, hit.channel_ID);

      if (key_found(shifts, detectorType))
      {
        cfd.calculate(shifts[detectorType], 0.75);
        
        auto cfdGraph = new TGraph(cfd.cfd.size(), Colib::linspaceFor(cfd.cfd, 0., 4.).data(), cfd.cfd.data());
        cfdGraph->SetLineColor(kGray);
        cfdGraph->Draw("same");

        if (key_found(thresholds, detectorType))
        {
          auto const & zero = cfd.findZero(thresholds[detectorType]) * 4;
          if (zero != CFD::noSignal)
          {
            TMarker *zero_marker = new TMarker(zero, 0, 20);
            zero_marker->SetMarkerStyle(29);
            zero_marker->SetMarkerColor(kGreen);
            zero_marker->Draw("SAME");
          }
        }
      // }
      canvas->Write();
    }
  }
  rootFile->Close();
  print("writeTraces.root written");
  print(timer());
  return 1;
}

int main(int argc, char** argv)
{
  if (argc == 1) {
    error("No arguments !");
    print("writeTraces usage : ./writeTraces /path/to/file nbHits [[param value]]");
    print();
    print("Parameters : ");
    print("adcMin [int] : sets the minimum adc value to store the trace");
    print("adcMax [int] : sets the maximum adc value to store the trace");
  }
  istringstream iss(argv_to_string(argv));

  string file = "empty.caendat";
  iss >> file >> file;

  int nb_hits = 0;
  double tmp_d; iss >> tmp_d;
  nb_hits = tmp_d;

  print(file, nb_hits);
  
  int adcMin = 0;
  int adcMax = -1;
  string temp; 
  while(iss >> temp)
  {
    if (temp == "adcMin") 
    {
      iss >> adcMin;
    }
    if (temp == "adcMax") 
    {
      iss >> adcMax;
    }
  }

  print(file, nb_hits, adcMin, adcMax);
  writeTraces(file, nb_hits, adcMin, adcMax);
}

// g++ -o writeTraces writeTraces.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2 -std=c++17