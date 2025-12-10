// The point of this code is to write down the traces of the signal
// The additionnal point is to study the cfd and compare with the original crrc2 timestamp,
// but without any coincidences, simply comparison of both methods !

// #include "CaenLib/CaenRawReader.hpp"
#include "CaenLib/CaenRootReader.hpp"
#include "AnalysisLib/CFD.hpp"
#include "AnalysisLib/RawHit.hpp"
#include "LibCo/Classes/Timer.hpp"
#include "LibCo/libCo.hpp"

#include "TPad.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph.h"

using namespace std;

// #define LOW_PASS
#ifdef LOW_PASS
  
  #include <fftw3.h>
  // #include "TFFTRealComplex.h"
  // #include "TFFTComplexReal.h"


  std::vector<int> lowpass_fft_root(const std::vector<int>& trace,
                                      double samplingRate,
                                      double cutoffFreq)
  {
      const int N = trace.size();
      if (N == 0) return {};

      // Convert trace to double
      std::vector<double> in(N);
      for (int i = 0; i < N; ++i)
          in[i] = static_cast<double>(trace[i]);

      // Forward FFT
      fftw_complex* freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
      fftw_plan fwd = fftw_plan_dft_r2c_1d(N, in.data(), freq, FFTW_ESTIMATE);
      fftw_execute(fwd);

      // Frequency bin width
      const double df = samplingRate / N;

      // Apply low-pass filter
      for (int i = 0; i < N / 2 + 1; ++i) {
          const double f = i * df;
          if (f > cutoffFreq) {
              freq[i][0] = 0.0;
              freq[i][1] = 0.0;
          }
      }

      // Inverse FFT
      fftw_plan inv = fftw_plan_dft_c2r_1d(N, freq, in.data(), FFTW_ESTIMATE);
      fftw_execute(inv);

      // Normalize inverse FFT
      for (int i = 0; i < N; ++i)
          in[i] /= N;

      fftw_destroy_plan(fwd);
      fftw_destroy_plan(inv);
      fftw_free(freq);

      std::vector<int> out; out.reserve(in.size());
      for (auto const & s : in) out.push_back(int_cast(s));

      return out; // keep as double for precision
  }
#endif //LOW_PASS

int writeTraces(string file, int nb_events_max = -1, int adcMin = 0, int adcMax = -1)
{
  auto rootFile = TFile::Open("writeTraces.root", "recreate"); rootFile->cd();

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


#ifdef LOW_PASS

      auto canvas0 = std::make_unique<TCanvas>((name+"filter").c_str(), (name+"filter").c_str()); canvas0->cd();
      
      auto dataGraph = hit.getTraceBaselineRemoved(10);

      auto graph1  = new TGraph(dataGraph  .size(), Colib::linspace<int>(dataGraph  .size(), 0, CaenDataReader1725::ticks_to_ns).data(), dataGraph  .data());
      auto dataGraphs2 = lowpass_fft_root(hit.getTraceBaselineRemoved(10), 1., 0.02);
      auto dataGraphs3 = lowpass_fft_root(hit.getTraceBaselineRemoved(10), 1., 0.03);
      auto dataGraphs4 = lowpass_fft_root(hit.getTraceBaselineRemoved(10), 1., 0.04);
      auto dataGraphs45 = lowpass_fft_root(hit.getTraceBaselineRemoved(10), 1., 0.045);
      auto dataGraphs5 = lowpass_fft_root(hit.getTraceBaselineRemoved(10), 1., 0.05);
      auto dataGraphs6 = lowpass_fft_root(hit.getTraceBaselineRemoved(10), 1., 0.06);

      auto graphs2 = new TGraph(dataGraphs2.size(), Colib::linspace<int>(dataGraphs2.size(), 0, CaenDataReader1725::ticks_to_ns).data(), dataGraphs2.data());
      auto graphs3 = new TGraph(dataGraphs3.size(), Colib::linspace<int>(dataGraphs3.size(), 0, CaenDataReader1725::ticks_to_ns).data(), dataGraphs3.data());
      auto graphs4 = new TGraph(dataGraphs4.size(), Colib::linspace<int>(dataGraphs4.size(), 0, CaenDataReader1725::ticks_to_ns).data(), dataGraphs4.data());
      auto graphs45 = new TGraph(dataGraphs45.size(), Colib::linspace<int>(dataGraphs45.size(), 0, CaenDataReader1725::ticks_to_ns).data(), dataGraphs45.data());
      auto graphs5 = new TGraph(dataGraphs5.size(), Colib::linspace<int>(dataGraphs5.size(), 0, CaenDataReader1725::ticks_to_ns).data(), dataGraphs5.data());
      auto graphs6 = new TGraph(dataGraphs6.size(), Colib::linspace<int>(dataGraphs6.size(), 0, CaenDataReader1725::ticks_to_ns).data(), dataGraphs6.data());

      graph1 -> SetTitle("Trace");
      graphs2 -> SetTitle("cutoff 0.02");
      graphs3 -> SetTitle("cutoff 0.03");
      graphs4 -> SetTitle("cutoff 0.04");
      graphs45 -> SetTitle("cutoff 0.045");
      graphs5 -> SetTitle("cutoff 0.05");
      graphs6 -> SetTitle("cutoff 0.06");
      
      graph1-> SetLineColor (kBlack);
      graphs2-> SetLineColor(kViolet);
      graphs3-> SetLineColor(kPink);
      graphs4-> SetLineColor(kBlue);
      graphs45-> SetLineColor(12);
      graphs5-> SetLineColor(kGreen);
      graphs6-> SetLineColor(8);

      graph1 -> Draw();
      graphs2 ->Draw("same");
      graphs3 ->Draw("same"); 
      graphs4 ->Draw("same"); 
      graphs45 ->Draw("same"); 
      graphs5 ->Draw("same"); 
      graphs6 ->Draw("same"); 
      
      canvas0->Write();

#endif //LOW_PASS

      auto canvas = std::make_unique<TCanvas>(name.c_str(), name.c_str()); canvas->cd();

      auto graphs = hit.getTracesGraphs(10);

      graphs[0] -> SetLineColor(kBlack);
      graphs[0] -> GetXaxis() -> SetTitle("time [ns]");
      graphs[0] -> GetYaxis() -> SetTitle("pulse height [ADC]");
      graphs[1] -> SetLineColor(kBlue );
      graphs[2] -> SetLineColor(kGreen);
      graphs[0] -> Draw();
      graphs[1] -> Draw("same");
      graphs[2] -> Draw("same");

      CFD cfd(*hit.trace, 2);

      auto const & detectorType = getDetectorType(hit.board_ID, hit.channel_ID);

      if (Colib::key_found(shifts, detectorType))
      {
        cfd.calculate(shifts[detectorType], 0.75);
        
        auto cfdGraph = new TGraph(cfd.cfd.size(), Colib::linspaceFor(cfd.cfd, 0., CaenDataReader1725::ticks_to_ns).data(), cfd.cfd.data());
        cfdGraph->SetLineColor(kGray);
        cfdGraph->Draw("same");

        if (Colib::key_found(thresholds, detectorType))
        {
          // auto const & zero = cfd.findZero() * CaenDataReader1725::ticks_to_ns;
          auto const & zero = cfd.findZero(thresholds[detectorType]) * CaenDataReader1725::ticks_to_ns;
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
    return 1;
  }
  istringstream iss(Colib::argv_to_string(argv));

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

// g++ -o writeTraces writeTraces.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2