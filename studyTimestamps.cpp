// The point of this code is to perform some event building on raw data, and play with the various information read in the data
// There is the possibility to play with traces as well, but studyCFD is a better code for this game
// Deprecated since the new RootCaenHit interface

#include "CaenReader.hpp"
#include "RawHit.hpp"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "../Nuball2/lib/Classes/Timer.hpp"

constexpr size_t nb_traces = 1000;

class TracesTH2F : public TH2F {
public:
  using TH2F::TH2F; // Inherit constructors

  void myFill(std::vector<CaenDataReader::Sample> const & samples) {
    if (m_y >= nb_traces) return;
    m_x = 0;
    for (auto const & sample : samples) TH2F::SetBinContent(m_x++, m_y, sample.sample);
    ++m_y;
  }

private:
  size_t m_x = 0; 
  size_t m_y = 0;
};

constexpr bool kBidims       = false;
constexpr bool kHandleTraces = true ;

constexpr Timestamp time_window          = 2e6 ; // ps
constexpr size_t    reserved_buffer_size = 5000ul;

int main(int argc, char** argv)
{
  Timer timer;
  int nb_events_max = -1;
  if (argc == 2) nb_events_max = int_cast(std::floor(std::stod(argv[1])));
  bool evts_limit = (nb_events_max>0);  

  /////////////////////////////
  // Prepare Root histograms //
  /////////////////////////////

  std::string rootFilename = "studyTimeshifts.root";
  auto rootFile = TFile::Open(rootFilename.c_str(), "recreate"); rootFile -> cd();

  auto fine_timestamp       = new TH1F("fine_timestamp"        , "fine_timestamp"     , 3000, 0, 3000);
  auto fine_timestamp_dssd  = new TH1F("fine_timestamp_dssd"   , "fine_timestamp_dssd", 4000, 0, 4000);
  auto fine_timestamp_used  = new TH1F("fine_timestamp_used"   , "fine_timestamp_used", 4000, 0, 4000);
  auto fine_timestamp_all   = new TH2F("fine_timestamp_all"    , "fine_timestamp_all" , 110,0,110, 4000,0,4000);

  auto GeHitPattern       = new TH1F("GeHitPattern"     , "GeHitPattern"      , 16, 0, 16);
  auto BGOHitPattern      = new TH1F("BGOHitPattern"    , "BGOHitPattern"     , 16, 0, 16);
  auto NedaHitPattern     = new TH1F("NedaHitPattern"   , "NedaHitPattern"    , 16, 0, 16);
  auto ringsHitPattern    = new TH1F("ringsHitPattern"  , "ringsHitPattern"   , 16, 0, 16);
  auto sectorsHitPattern  = new TH1F("sectorsHitPattern", "sectorsHitPattern" , 32, 0, 32);
  
  auto Neda_Qshort_Qlong = new TH2F("Neda_Qshort_Qlong", "Neda_Qshort_Qlong;short;long", 1000,0,100000, 1000,0,10000);
  
  auto HitPattern = new TH2F("HitPattern", "HitPattern", 110,0,110, 110,0,110);

  auto time_ref_extended  = new TH2F("time_ref_extended", "time_ref_extended;label;ps", 110,0,110, 4001, -2e6, 2e6);//1+int_cast(time_window*1e-3),-time_window,time_window);
  auto time_ref_precise   = new TH2F("time_ref_precise" , "time_ref_precise;label;ps" , 110,0,110, 4001, -2e6, 2e6);//1+int_cast(time_window*1e-3),-time_window,time_window);

  auto time_Ge_Ge_extended      = new TH1F("time_Ge_Ge_extended "     , "time_Ge_Ge_extended;ps"      , 1000,-10000,10000);
  auto time_Ge_Ge_precise       = new TH1F("time_Ge_Ge_precise"       , "time_Ge_Ge_precise;ps"       , 1000,-10000,10000);
  auto time_Neda_Neda_extended  = new TH1F("time_Neda_Neda_extended " , "time_Neda_Neda_extended ;ps" , 1000,-10000,10000);
  auto time_Neda_Neda_precise   = new TH1F("time_Neda_Neda_precise"   , "time_Neda_Neda_precise;ps"   , 1000,-10000,10000);

  std::array<std::vector<TH1F*>, 5> energies;
  std::array<std::vector<TracesTH2F*>, 5> traces;
  for (int label = 0; label<32; ++label)
  {
    auto const & label_str = std::to_string(label);
    energies[DSSDSector].push_back(new TH1F(("Sector_"+label_str+"_energy").c_str(), ("Sector_"+label_str+"_energy").c_str(), 10000,0,100000));
    if (kHandleTraces) traces[DSSDSector].push_back(new TracesTH2F(("Sector_"+label_str+"_traces").c_str(), ("Sector_"+label_str+"_traces").c_str(), nb_traces,0,nb_traces, 2000,0,2000));
    if (label > 15) continue;
    energies[EAGLE]   .push_back(new TH1F(("Ge_"  +label_str+"_energy").c_str(), ("Ge_"  +label_str+"_energy").c_str(), 10000,0,100000));
    energies[NEDA]    .push_back(new TH1F(("Neda_"+label_str+"_energy").c_str(), ("Neda_"+label_str+"_energy").c_str(), 10000,0,100000));
    energies[DSSDRing].push_back(new TH1F(("Ring_"+label_str+"_energy").c_str(), ("Ring_"+label_str+"_energy").c_str(), 10000,0,100000));
    energies[BGO]     .push_back(new TH1F(("BGO_" +label_str+"_energy").c_str(), ("BGO_" +label_str+"_energy").c_str(), 10000,0,100000));
    if (kHandleTraces){
      traces[EAGLE]   .push_back(new TracesTH2F(("Ge_"  +label_str+"_traces").c_str(), ("Ge_"  +label_str+"_traces").c_str(), nb_traces,0,nb_traces, 2000,0,2000));
      traces[NEDA]    .push_back(new TracesTH2F(("Neda_"+label_str+"_traces").c_str(), ("Neda_"+label_str+"_traces").c_str(), nb_traces,0,nb_traces, 2000,0,2000));
      traces[DSSDRing].push_back(new TracesTH2F(("Ring_"+label_str+"_traces").c_str(), ("Ring_"+label_str+"_traces").c_str(), nb_traces,0,nb_traces, 2000,0,2000));
      traces[BGO]     .push_back(new TracesTH2F(("BGO_" +label_str+"_traces").c_str(), ("BGO_" +label_str+"_traces").c_str(), nb_traces,0,nb_traces, 2000,0,2000));
    }
  }

  /////////////////////
  // Read Aggregates //
  /////////////////////

  size_t size_buffer = 0ul;
  std::vector<RawHit> hits_buffer; hits_buffer.reserve(reserved_buffer_size);
  std::vector<RawHit> events;

  int nb_evts = 0;
  bool finished = false;
  std::vector<std::string> filenames = 
  {
    // "/home/corentin/coulexRU_i2097_3020_0000.caendat",
    // "/home/corentin/coulexRU_i2097_3020_0001.caendat",
    // "/home/corentin/coulexRU_i2097_3020_0002.caendat"

    "/home/corentin/60Co_data/eagleRU_i2514_0023_0000.caendat"

  };

  for (auto const & filename : filenames)
  {
    if (finished) break;
    print("Reading", filename);
    RawCaenReader reader(filename);
    reader.handleTraces(kHandleTraces);
    while(reader.readBoardAggregate() && !finished)
    {
      auto const & board = reader.getBoard();
      
      //////////////////
      // Hits storing //
      //////////////////
      
      if (size_buffer + board.nbEvents() < reserved_buffer_size)
      {
        for (auto const & channel : board.channels) for (auto const & event : channel.events)
        {
          if (evts_limit && nb_evts > nb_events_max) finished = true;
          if (++nb_evts % int_cast(1.e6) == 0) print(nicer_double(nb_evts));
          auto extra = dynamic_cast<CaenDataReader::FineTimestamp*>(event.extra);
          if (extra) 
          {
            fine_timestamp -> Fill(extra->fine_timestamp / 4. * 1.024);
            fine_timestamp_used -> Fill(extra->fine_timestamp);
            if (board.BOARD_ID>=6) fine_timestamp_dssd -> Fill(extra->fine_timestamp);
            fine_timestamp_all -> Fill(gLabels(board, channel, event), extra->fine_timestamp);
          }
          if (kHandleTraces) {
            traces[getDetectorType(board, event)][labels(board, channel, event)] -> myFill(event.samples);
          }
          hits_buffer.emplace_back(board, channel, event);
        }
        size_buffer += board.nbEvents();
        if (kHandleTraces) for (auto & channel : reader.getBoard().channels) for (auto & event : channel.events) event.samples.clear();
        if (!finished) continue; // if the max number of iterations is reached, stop the buffer filling
      }

      // This part of the code is accessed only when the buffer is full, or the max number of iterations is reached
      size_buffer = 0; // Reinitialize for next hits_buffer filling
  
      ////////////////////
      // Event building //
      ////////////////////
  
      // Time sorting :
      auto ordered_index = bubble_sort(hits_buffer);
          
      // Some event building informations : 
      constexpr Label trigg_label = 48;
      bool trigger = false;
      size_t trigg_i = 0;
      uint64_t timestamp_extended     = 0;
      uint64_t timestamp_precise = 0;
      
      auto pre_info = [&trigg_label, &trigger, &timestamp_extended, &timestamp_precise, &trigg_i](RawHit const & hit, size_t const & i){
        if (hit.glabel == trigg_label){
          trigger = true;
          trigg_i = i;
          timestamp_extended = hit.extended_timestamp;
          timestamp_precise = hit.precise_timestamp;
        }
      };
      
      // Initialize iteration with first hit : 
      events.clear();
      events.emplace_back(hits_buffer[ordered_index[0]]);
      pre_info(events.front(), 0);
      for (size_t buff_loop_i = 1; buff_loop_i<ordered_index.size(); ++buff_loop_i)
      {
        auto const & buff_hit_i = ordered_index[buff_loop_i]; // Time ordered index
        auto const & buff_hit   = hits_buffer  [buff_hit_i ];  // Time ordered hit
  
        // Try to fill the event with the other time-coincident hits :
        if ((buff_hit.precise_timestamp - events.front().precise_timestamp) < time_window) {
          events.push_back(buff_hit);
          pre_info(buff_hit, buff_loop_i);
          continue; // If fill is successfull, iterate again with next hit
        }
  
        // This part of the code is accessed only if event filling failed, hence the current hit is out of event building time window
        // The hit will be filled at the end of the analysis section with the current hit, for next iteration to build next event
  
        ////////////////////
        // Event Analysis //
        ////////////////////
  
        for (size_t loop_i = 0; loop_i<events.size(); ++loop_i)
        {
          auto const & hit = events[loop_i];
          
          energies[hit.detectorType][hit.label] -> Fill(hit.adc);

          // if (buff_hit_i) print(int(hit.glabel));
          
          if (trigger && trigg_i != loop_i) 
          {
            auto const & dt_raw     = int_cast(timestamp_extended - hit.extended_timestamp);
            auto const & dt_precise = int_cast(timestamp_precise  - hit.precise_timestamp );
            time_ref_extended -> Fill(hit.glabel, dt_raw    );
            time_ref_precise  -> Fill(hit.glabel, dt_precise);
          }
          switch (hit.detectorType)
          {
          case EAGLE      : GeHitPattern      ->  Fill(hit.label); break;
          case BGO        : BGOHitPattern     ->  Fill(hit.label); break;
          case NEDA       : NedaHitPattern    ->  Fill(hit.label); break;
          case DSSDRing   : ringsHitPattern   ->  Fill(hit.label); break;
          case DSSDSector : sectorsHitPattern ->  Fill(hit.label); break;
          default: break;
          }

          if (hit.detectorType == NEDA){
            Neda_Qshort_Qlong->Fill(hit.adc, hit.adc2);
          }
          
          if (kBidims) for (size_t loop_j = loop_i + 1; loop_j<events.size(); ++loop_j)
          {
            auto const & event_j = events[loop_j];

            HitPattern->Fill(hit    .glabel, event_j.glabel);
            HitPattern->Fill(event_j.glabel, hit    .glabel);
  
            auto const & dT_exended = int_cast(hit.extended_timestamp - event_j.extended_timestamp)/1000.;
            auto const & dT_precise = int_cast(hit.precise_timestamp  - event_j.precise_timestamp )/1000.;
  
            switch (event_j.detectorType) {
            case EAGLE:
              if (hit.detectorType == EAGLE) 
              {
                time_Ge_Ge_extended ->  Fill(dT_exended);
                time_Ge_Ge_precise  ->  Fill(dT_precise);
              }
              break;
  
            case NEDA:
              if (hit.detectorType == NEDA) 
              {
                time_Neda_Neda_extended ->  Fill(dT_exended);
                time_Neda_Neda_precise  ->  Fill(dT_precise);
              }
              break;          
            
            default:
              break;
            }
          }
        }
        
        ////////////////////////
        // Prepare next event //
        ////////////////////////
        
        events.clear();
        events.push_back(buff_hit);
        trigger             = false;
        trigg_i             = 0;
        timestamp_extended  = 0;
        timestamp_precise   = 0;
      }
      // Filling the new hits_buffer with the current out-of-bound board aggregate :
      hits_buffer.clear();
      for (auto const & channel : board.channels) for (auto const & event : channel.events) 
        hits_buffer.emplace_back(board, channel, event);
    }  
  }

    // fine_timestamp  ->  Write();
    // fine_timestamp_used  ->  Write();

    // HitPattern        ->  Write();
    // ringsHitPattern   ->  Write();
    // sectorsHitPattern ->  Write();
    // GeHitPattern      ->  Write();
    // BGOHitPattern     ->  Write();
    // NedaHitPattern    ->  Write();

    // Neda_Qshort_Qlong    ->  Write();

    // time_ref_extended      ->  Write();
    // time_ref_precise  ->  Write();

    // time_Ge_Ge_extended           ->  Write();
    // time_Ge_Ge_precise      ->  Write();
    // time_Neda_Neda_extended       ->  Write();
    // time_Neda_Neda_precise  ->  Write();

    // for (auto & histos : energies) for (auto & histo : histos) histo->Write();
  rootFile->Write();
  rootFile->Close();
  print(rootFilename, "written");

  print(timer());
  return 0;
}

// g++ -o exec studyTimestamps.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs`




/*
for (auto const & board : boards)
    {
      auto const & boardType = Boards_map[board.BOARD_ID];

      for (auto const & channel : board.channels)
      {
        for (auto const & event : channel.events)
        {
          switch (boardType)
          {
          case EAGLE:{
            // ID is 0 to 16 : board.BOARD_ID [0;1], channel.ID [0,7]
            auto const & EagleID = board.BOARD_ID * 8 + channel.ID;
            if (event.CH) 
            { // BGOs :
              energies[4].at(EagleID)->Fill(event.ENERGY);
              if (trigger) time_ref[4].at(EagleID)->Fill(timestamp - event.TRIGGER_TIME_TAG);
              BGOHitPattern->Fill(EagleID);
            }
            else          
            { // Ge :
              energies[boardType].at(EagleID)->Fill(event.ENERGY);
              if (trigger) time_ref[boardType].at(EagleID)->Fill(timestamp - event.TRIGGER_TIME_TAG);
              GeHitPattern ->Fill(EagleID);
            }
            break;
          }

          case NEDA:{
            // ID is 0 to 16 : channel.ID [0,7], event.CH [0;1]
            auto const & NedaID = channel.ID * 2 + int_cast(event.CH);
            NedaHitPattern->Fill(NedaID);
            energies[boardType].at(NedaID)->Fill(event.ENERGY);
            if (trigger) time_ref[boardType].at(NedaID)->Fill(timestamp - event.TRIGGER_TIME_TAG);
            break;
          }

          case DSSDRing:{
            // ID is 0 to 16 : channel.ID [0,7], event.CH [0;1]
            auto const & RingID = channel.ID * 2 + int_cast(event.CH);
            ringsHitPattern->Fill(RingID);
            energies[boardType].at(RingID)->Fill(event.ENERGY);
            
            break;
          }

          case DSSDSector:{
            // ID is 0 to 32 : board.BOARD_ID [7;8], channel.ID [0,7], event.CH [0,1]
            auto const & SectorID = (board.BOARD_ID - 7) * 16 + channel.ID * 2 + int_cast(event.CH);
            sectorsHitPattern->Fill(SectorID);
            energies[boardType].at(SectorID)->Fill(event.ENERGY);
            if (trigger) time_ref[boardType].at(SectorID)->Fill(timestamp - event.TRIGGER_TIME_TAG);
            break;
          }

          default:
            error("Do not handle board", boardType);
            break;
          }
        }
      }
    }
    
*/