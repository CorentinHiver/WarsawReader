#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "AnalysisLib/CFD.hpp"
#include "LibCo/Timer.hpp"

#include "CaenLib/CaenRootReader.hpp"
#include "CaenLib/CaenRootEventBuilder.hpp"

constexpr Long64_t time_window          = 2e6 ; // ps
constexpr size_t   reserved_buffer_size = 5000ul;
auto constexpr static ref_label = 81;

enum DetectorTypes{ EAGLE, BGO, NEDA, DSSDRing, DSSDSector, LaBr3, EMPTY};
static const std::array<std::string, 7> DetectorName = {"EAGLE", "BGO", "NEDA", "DSSDRing", "DSSDSector", "LaBr3", "EMPTY"};
static constexpr std::array<int, 10> Boards_map = {EAGLE, EAGLE, EMPTY, EMPTY, EMPTY, NEDA, DSSDRing, DSSDSector, DSSDSector, LaBr3};

int studyCFD(int nb_events_max = -1)
{
  Timer timer;
  bool max_events = (nb_events_max>0);
  std::unordered_map<int, int> cfd_shifts = {
    {0, 5},
    {1, 5},
    {6, 2},
    {7, 2},
    {8, 2}
  };

  std::unordered_map<int, double> cfd_thresholds = {
    {0, -50},
    {1, -50},
    {6, -500},
    {7, -500},
    {8, -500}
  };

  // std::vector<std::string> filenames = {"/home/corentin/60Co_data/eagleRU_i2514_0023_0000.caendat"};
  // std::vector<std::string> filenames = {"/home/corentin/60Co_data/eagleRU_i2607_0005_0000.caendat"};
  // std::vector<std::string> filenames = {
  //   "/home/corentin/60Co_data/eagleRU_i2606_0004_0000.caendat",
  //   "/home/corentin/60Co_data/eagleRU_i2606_0004_0001.caendat"
  // };

  auto constexpr static glabel = [](RootCaenHit const & hit){
    return hit.board_ID * 16 + hit.channel_ID * 2 + hit.subchannel_ID;
  };

  auto HitPattern   = new TH1F("HitPattern", "HitPattern", 200,0,200);
  auto HitPattern2D = new TH2F("HitPattern2D", "HitPattern2D", 200,0,200, 200,0,200);
  
  auto E_all = new TH2F("E_all", "E_all", 200,0,200, 10000,0,100000);
  auto neda_psd = new TH2F("neda_psd", "neda_psd", 16,0,16, 1000,-5,5);
  std::vector<TH2F*> neda_psds; neda_psds.reserve(16);
  for (int i = 0; i<16; ++i) 
  {
    TString name ("neda_psd"+std::to_string(i));
    neda_psds.push_back(new TH2F(name, name, 1000,-2*time_window,2*time_window, 1000,-5,5));
  }

  auto cfd_corrections = new TH2F("cfd_corrections", "cfd_corrections", 40,0,40, 400000,-2000000,2000000);

  auto dT_all = new TH1F("dT_all", "dT_all", 10000,-2*time_window,2*time_window);
  auto dT_Ge1_Ref = new TH1F("dT_Ge1_Ref", "dT_Ge1_Ref", 10000,-2*time_window,2*time_window);
  auto E_VS_dT_Ge1_Ref = new TH2F("E_VS_dT_Ge1_Ref", "E_VS_dT_Ge1_Ref", 1000,-2*time_window,2*time_window, 10000,0,100000);
  auto dT_Ref_VS_all = new TH2F("dT_Ref_VS_all", "dT_Ref_VS_all", 200,0,200, 10000,-2*time_window,2*time_window);
  
  auto dT_all_cfd = new TH1F("dT_all_cfd", "dT_all_cfd", 10000,-2*time_window,2*time_window);
  auto cfd_dT_Ge1_Ref = new TH1F("cfd_dT_Ge1_Ref", "cfd_dT_Ge1_Ref", 10000,-2*time_window,2*time_window);
  auto E_VS_dT_Ge1_Ref_cfd = new TH2F("E_VS_dT_Ge1_Ref_cfd", "E_VS_dT_Ge1_Ref_cfd", 10000,-2*time_window,2*time_window, 10000,0,100000);
  auto dT_Ref_VS_all_cfd = new TH2F("dT_Ref_VS_all_cfd", "dT_Ref_VS_all_cfd", 200,0,200, 10000,-2*time_window,2*time_window);

  std::vector<TH2F*> dT_all_vs_all; dT_all_vs_all.reserve(200);
  std::vector<TH2F*> dT_all_vs_all_cfd; dT_all_vs_all_cfd.reserve(200);
  for (size_t board_i = 0; board_i<Boards_map.size(); ++board_i) for (size_t channel_i = 0; channel_i<16; ++channel_i)
  {
    if (Boards_map[board_i] == EMPTY) 
    {
      dT_all_vs_all.push_back(nullptr);
      dT_all_vs_all_cfd.push_back(nullptr);
    }
    else
    {
      TString name ("dT_" + DetectorName[Boards_map[board_i]] + "_" + std::to_string(channel_i) + "_label_" + std::to_string(board_i*16+channel_i));
      dT_all_vs_all.push_back(new TH2F(name, name, 200,0,200, 2000,-time_window,time_window));
      TString name2 ("cfd_dT_" + DetectorName[Boards_map[board_i]] + "_" + std::to_string(channel_i) + "_label_" + std::to_string(board_i*16+channel_i));
      dT_all_vs_all_cfd.push_back(new TH2F(name2, name2, 200,0,200, 2000,-time_window,time_window));
    }
  }

  for (auto const & filename : filenames)
  {
    CaenRootReader reader(filename); //reader.handleTraces(false); //(not handling traces removes cfd handling)
    CaenRootEventBuilder event_builder(reserved_buffer_size);
    while(((max_events) ? (reader.nbHits() < nb_events_max) : (true)) && reader.readHit())
    {
      ////////////////////
      // Buffer Filling //
      ////////////////////

      auto & hit = reader.getHit();

      if (reader.nbHits() % int(1e5) == 0) print(nicer_double(reader.nbHits(), 1));

      // Correct timestamp with cfd :

      if (!hit.outTrace.empty() && key_found(hit.board_ID, cfd_shifts)) 
      {
        CFD cfd(hit.outTrace, cfd_shifts[hit.board_ID], 0.75);
        
        auto zero = cfd.findZero(cfd_thresholds[hit.board_ID]); 
        
        if (zero == CFD::noSignal || zero == CFD::noZero) continue;

        zero = zero * 4000.; // Convert from 4 ns ticks to ps

        cfd_corrections->Fill(glabel(hit), zero); 

        hit.cfd = hit.extended_ts + zero; 
      }
      else hit.cfd = hit.precise_ts;

      // Filter bad hits :

      // if (hit.adc < 10) continue;

      if (((max_events) ? (reader.nbHits() < nb_events_max) : (true)) && event_builder.fill_buffer(hit)) continue;

      ////////////////////
      // Event Building //
      ////////////////////

      event_builder.fast_event_building(time_window);

      ////////////////////
      // Event Analysis //
      ////////////////////

      for (auto const & event : event_builder)
      {
        auto const & first_hit = event_builder[event[0]];
        for (size_t hit_i = 0; hit_i<event.size(); ++hit_i) // Loop in the event
        {
          auto const & hit_index_i = event[hit_i];          // The index in the buffer
          auto const & hit_0 = event_builder[hit_index_i];  // The hit itself
          auto const & glabel_0 = glabel(hit_0);            // The global index (16*board_ID + 2*channel_ID + subchannel_ID)
          
          HitPattern->Fill(glabel_0);
          E_all->Fill(glabel_0, hit_0.adc);

          if (Boards_map[hit_0.board_ID] == NEDA && hit_0.qlong!=0)
          {
            auto const & channel_label = hit_0.channel_ID*2+hit_0.subchannel_ID;
            auto const & psd = hit_0.adc/float_cast(hit_0.qlong+hit_0.adc) - 1;
            neda_psd->Fill(channel_label, psd);
            neda_psds[channel_label]->Fill(hit_0.cfd - first_hit.cfd, psd);
          }

          for (size_t hit_j = hit_i+1; hit_j<event.size(); ++hit_j)
          {
            auto const & hit_index_j = event[hit_j];
            auto const & hit_1 = event_builder[hit_index_j];            
            auto const & glabel_1 = glabel(hit_1);

            HitPattern2D->Fill(glabel_0, glabel_1);
            HitPattern2D->Fill(glabel_1, glabel_0);

            dT_all_vs_all[glabel_0]->Fill(glabel_1, Long64_t(hit_1.timestamp - hit_0.timestamp));
            dT_all_vs_all[glabel_1]->Fill(glabel_0, Long64_t(hit_0.timestamp - hit_1.timestamp));

            dT_all_vs_all_cfd[glabel_0]->Fill(glabel_1, Long64_t(hit_1.cfd - hit_0.cfd));
            dT_all_vs_all_cfd[glabel_1]->Fill(glabel_0, Long64_t(hit_0.cfd - hit_1.cfd));

            dT_all->Fill(Long64_t(hit_1.timestamp - hit_0.timestamp));
            dT_all_cfd->Fill(Long64_t(hit_1.cfd - hit_0.cfd));
          }

          // Time reference histograms :
          if (glabel_0 == ref_label) for (size_t hit_j = 0; hit_j<event.size(); ++hit_j)
          {
            if (hit_j == hit_i) continue;

            auto const & hit_index_j = event[hit_j];
            auto const & hit_1 = event_builder[hit_index_j];            
            auto const & glabel_1 = glabel(hit_1);

            auto const & dT = Long64_t(hit_1.timestamp - hit_0.timestamp);
            auto const & dT_cfd = Long64_t(hit_1.cfd - hit_0.cfd);
            
            dT_Ref_VS_all     -> Fill(glabel_1, dT    );
            dT_Ref_VS_all_cfd -> Fill(glabel_1, dT_cfd);

            if (glabel_1 == 0) 
            {
              dT_Ge1_Ref -> Fill(dT);
              cfd_dT_Ge1_Ref -> Fill(dT_cfd);
              
              E_VS_dT_Ge1_Ref     -> Fill(dT    , hit_1.adc);
              E_VS_dT_Ge1_Ref_cfd -> Fill(dT_cfd, hit_1.adc);
            }            
          }
        }
      }
      event_builder.clear();
    }
  }

  auto rootfile = TFile::Open("studyCFD.root", "recreate"); rootfile->cd();

    HitPattern    ->  Write();
    HitPattern2D  ->  Write();

    dT_all           ->  Write();
    dT_Ref_VS_all    ->  Write();
    dT_Ge1_Ref       ->  Write();
    E_VS_dT_Ge1_Ref  ->  Write();

    dT_all_cfd          -> Write();
    cfd_dT_Ge1_Ref      -> Write();
    E_VS_dT_Ge1_Ref_cfd -> Write();
    dT_Ref_VS_all_cfd   -> Write();
    
    cfd_corrections ->  Write();
    
    E_all -> Write();
    neda_psd -> Write();
    
    auto filledHisto = []<class THist>(THist * histo){return histo && !histo->IsZombie() && histo->Integral()>0;};

    for (auto & histo : neda_psds         ) if (filledHisto(histo)) histo->Write();
    for (auto & histo : dT_all_vs_all     ) if (filledHisto(histo)) histo->Write();
    for (auto & histo : dT_all_vs_all_cfd ) if (filledHisto(histo)) histo->Write();

  rootfile->Close();
  print("studyCFD.root written");

  // This is not normal, normally the root file takes ownership of the histograms and delete them when Close is called
  delete HitPattern;
  delete HitPattern2D;
  delete E_all;
  delete neda_psd;
  delete dT_all;
  delete dT_Ref_VS_all;
  delete dT_Ge1_Ref;
  delete E_VS_dT_Ge1_Ref;
  delete dT_all_cfd;
  delete cfd_dT_Ge1_Ref;
  delete E_VS_dT_Ge1_Ref_cfd;
  delete dT_Ref_VS_all_cfd;

  for (auto & histo : neda_psds) delete histo;
  for (auto & histo : dT_all_vs_all) delete histo;
  for (auto & histo : dT_all_vs_all_cfd) delete histo;

  print(timer());
  return 0;
}

int main(int argc, char** argv)
{
  if (argc == 2) return studyCFD(int_cast(std::floor(std::stod(argv[1]))));
  else return studyCFD();
}


// g++ -o exec studyCFD.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2