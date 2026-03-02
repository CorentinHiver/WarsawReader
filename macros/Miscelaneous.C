#include "../CaenLib/RootReader.hpp"
#include "../LibCo/Classes/Timeshifts.hpp"
#include "TH2F.h"

void GeSpectrums(TFile* file, std::string const & calibName = "")
{
  if (calibName != "") Colib::throw_error("Calibration option not coded yet!");
  Caen1725::RootReader reader(file);
  auto & event = reader.getEvent();
  auto h2 = new TH2F("GeSpectrumsH",";label;energy[ADC]",20,0,20, 3500,0,35000);
  while(reader.readNextEvent()) for (size_t hit_i = 0; hit_i < event.size(); ++hit_i)
  {// Loop through the hits
    if (event.board_ID[hit_i]<3 && event.subchannel_ID[hit_i]==0) h2->Fill(event.label[hit_i]/2, event.adc[hit_i]);
  }
  print("GeSpectrumsH created, GeSpectrumsH->Draw() to see it");
}

void dT(TFile* file, int ref_label=81)
{
  Caen1725::RootReader reader(file);
  auto & event = reader.getEvent();
  auto dT = new TH2F("hdT","dT;label;dT[ps]",200,0,200, 6000,-3000000,3000000);
  while(reader.readNextEvent())
  { // Looping through the events
    for (size_t trigger_i = 0; trigger_i < event.size(); ++trigger_i) if (event.label[trigger_i] == ref_label) 
    { // Looping through the hits of the event and proceed only if  a hit in the reference detector have been found
      for (size_t hit_i = 0; hit_i < event.size(); ++hit_i) if (hit_i !=trigger_i)
      {// Looping again through the hits of the event and calculate time diff between the reference and every other detector
        dT->Fill(event.label[hit_i], int(event.time[hit_i]-event.time[trigger_i]));
      }
      break;
    }
  }
  print("hdT created, hdT->Draw() to see it");
}

void dT(std::string const & filename, int ref_label = 81)
{
  dT(TFile::Open(filename.c_str(), "READ"), ref_label);
}

void coincMatrix(TFile* file)
{
  Caen1725::RootReader reader(file);
  auto & event = reader.getEvent();

  auto coinc = new TH2F("coincMatrixH",";label;label",200,0,200, 200,0,200);
  while(reader.readNextEvent())
  { // Looping through the events
    for (size_t hit_i = 0; hit_i < event.size(); ++hit_i)
    { // Looping through the hits of the event
      for (size_t hit_j = hit_i+1; hit_j < event.size(); ++hit_j)
      {// Looping through the other hits of the event
        coinc->Fill(event.label[hit_i], event.label[hit_j]);
        coinc->Fill(event.label[hit_j], event.label[hit_i]);
      }
      break;
    }
  }
  print("coincMatrixH created, coincMatrix->Draw() to see it");
  reader.resetCursor();
}

void forMarcin(std::string filename, std::string tsFile)
{
  auto file = TFile::Open(filename.c_str(), "READ");
  if (!file) Colib::throw_error(filename+" not found");
  Caen1725::RootReader reader(file);
  auto & event = reader.getEvent();
  auto dT = new TH2F("hdT","dT;HPGe_label;dT[ps]", 6000,-3000000,3000000, 20,0,20);
  Timeshifts ts(tsFile);
  while(reader.readNextEvent()) for (size_t dssd_i = 0; dssd_i < event.size(); ++dssd_i) if (7<event.board_ID[dssd_i])
  {// Loop through the hits
    auto const & label_dssd = event.label[dssd_i];
    event.time[dssd_i]+=ts.get(label_dssd);
    for (size_t Ge_i = 0; Ge_i < event.size(); ++Ge_i) if (event.board_ID[Ge_i] < 3 && event.subchannel_ID[Ge_i]==0)
      dT->Fill(event.time[Ge_i] - event.time[dssd_i], event.board_ID[Ge_i]*8 + event.channel_ID[Ge_i]);
  }
  print("hdT created, hdT->Draw() to see it");
}