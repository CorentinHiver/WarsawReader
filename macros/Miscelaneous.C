#include "../CaenLib/RootReader.hpp"
#include "TH2F.h"

void dT(TFile* file, int ref_label=81)
{
  RootReader reader(file);
  auto & event = reader.getEvent();
  // print("Sets the number of time bins between -2us and +2us. Default: 4000 (length = 1ns).")
  auto dT = new TH2F("hdT","dT;label;dT[ps]",200,0,200, 4000,-2000000,2000000);
  while(reader.readNextEvent())
  { // Looping through the events
    for (size_t trigger_i = 0; trigger_i < event.size(); ++trigger_i) if (event.label[trigger_i] == ref_label) 
    { // Looping through the hits of the event and proceed only if  a hit in the reference detector have been found
      for (size_t hit_i = 0; hit_i < event.size(); ++hit_i) if (hit_i !=trigger_i)
      {// Looping again through the hits of the events and calculate time diff between the reference and every other detector
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