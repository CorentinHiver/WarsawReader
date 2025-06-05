// The point of this code is to write down the traces of the signal
// The additionnal point is to study the cfd and compare with the original crrc2 timestamp,
// but without any coincidences, simply comparison of both methods !

#include "CaenLib/CaenRawReader.hpp"
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

int writeTraces(int nb_events_max = -1)
{
  Timer timer;
  bool max_events = (nb_events_max>0);
  bool finished = false;
  int nb_evts = 0;
  std::unordered_map<int, size_t> shifts = {
    {EAGLE     , 5},
    {DSSDSector, 2},
    {DSSDRing  , 2}
  };
  std::unordered_map<int, int> thresholds = {
    {EAGLE     , -50},
    {DSSDSector, -500},
    {DSSDRing  , -500}
  };

  // CaenRawReader reader("/home/corentin/coulexRU_i2097_3020_0000.caendat");
  // CaenRawReader reader("/home/corentin/60Co_data/eagleRU_i2608_0006_0000.caendat");
  // CaenRawReader reader("/home/corentin/60Co_data/eagleRU_i2514_0023_0000.caendat");
  CaenRawReader reader("../../data/coulexNov2024/coulexRU_i2097_3020_0000.caendat");
  reader.handleTraces(true);

  auto rootFile = TFile::Open("writeTraces.root", "recreate"); rootFile->cd();

  std::vector<std::string> folders_names;

  while(!finished && reader.readBoardAggregate())
  {
    auto const & board = reader.getBoard();
    for (auto const & channel : board.channels) for (auto const & event : channel.events){
      ++nb_evts;
      if (nb_evts % int(1e3) == 0) print(nicer_double(nb_evts, 1));
      if (max_events && nb_evts > nb_events_max) finished = true;

      RawHit hit(board, channel, event);
      // print(DetectorName[hit.detectorType], int(hit.label), int(hit.glabel), hit.adc, hit.precise_timestamp);

      if (event.samples.size() == 0) continue;

      auto detectorName = DetectorName[getDetectorType(board, event)]+std::to_string(labels(board, channel, event));
      TDirectory *dir = rootFile->GetDirectory(detectorName.c_str());

      if (!dir){
        rootFile->mkdir(detectorName.c_str());
        dir = rootFile->GetDirectory(detectorName.c_str());
      }
      dir->cd();
      std::string name = detectorName+"_"+std::to_string(nb_evts);
      std::vector<uint16_t> raw_trace; raw_trace.reserve(event.samples.size());
      std::vector<double> samplesT; samplesT.resize(event.samples.size(), 0);
      std::vector<double> samplesDP1; samplesDP1.resize(event.samples.size(), 0);
      for (size_t sample_i = 0; sample_i<event.samples.size(); ++sample_i) 
      {
        auto const & sample = event.samples[sample_i];
        raw_trace.push_back(sample.sample);
        if (sample.T) samplesT[sample_i] = 1;
        if (sample.DP1) samplesDP1[sample_i] = 1;
      }
      
      CFD cfd(raw_trace);

      auto const & maxSample = maximum(cfd.trace);
      for (auto & Ti : samplesT) if (Ti > 0.5) Ti = maxSample;
      for (auto & DP1 : samplesDP1) if (DP1 > 0.5) DP1 = maxSample;

      auto graph  = new TGraph(cfd.trace.size(), linspace_for(cfd.trace, 0., 4.).data(), cfd.trace.data());
      auto graph2 = new TGraph(samplesT .size(), linspace_for(samplesT , 0., 4.).data(), samplesT .data());
      auto graph3 = new TGraph(samplesDP1 .size(), linspace_for(samplesDP1 , 0., 4.).data(), samplesDP1 .data());
      
      auto canvas = new TCanvas(name.c_str(), name.c_str()); canvas->cd();
      graph->Draw();
      graph2->SetLineColor(kGreen);
      graph->GetXaxis()->SetTitle("time [ns]");
      graph->GetYaxis()->SetTitle("pulse height [ADC]");
      graph2->Draw("same");
      graph3->SetLineColor(kBlue);
      graph3->Draw("same");
      if (key_found(shifts, hit.detectorType))
      {
        cfd.calculate(shifts[hit.detectorType], 0.5);
        auto graph4 = new TGraph(cfd.cfd.size(), linspace_for(cfd.cfd, 0., 4.).data(), cfd.cfd.data());
        graph4->SetLineColor(kGray);
        graph4->Draw("same");

        if (key_found(thresholds, hit.detectorType))
        {
          auto const & zero = cfd.findZero(thresholds[hit.detectorType]) * 4;
          if (zero != CFD::noSignal)
          {
            TMarker *zero_marker = new TMarker(zero, 0, 20);
            zero_marker->SetMarkerStyle(29);
            zero_marker->SetMarkerColor(kGreen);
            zero_marker->Draw("SAME");
          }
        }
      }
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
  if (argc == 2) writeTraces(int_cast(std::floor(std::stod(argv[1]))));
  else return writeTraces();
}

// g++ -o exec writeTraces.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2