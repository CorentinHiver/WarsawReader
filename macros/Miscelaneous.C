#include "../CaenLib/RootReader.hpp"
#include "TH2F.h"

void dT(TFile* file, int ref_label=81)
{
  RootReader reader(file);
  auto & event = reader.getEvent();
  auto dT = new TH2F("hdT","dT;label;dT[ps]",200,0,200, 4000,-2000000,2000000);
  while(reader.readNextEvent())
  {
    for (size_t trigger_i = 0; trigger_i < event.size(); ++trigger_i) if (event.label[trigger_i] == ref_label) 
    {
      for (size_t hit_i = 0; hit_i < event.size(); ++hit_i) if (hit_i !=trigger_i)
      {
        // print(event.label[hit_i], event.time[hit_i]-event.time[trigger_i]);
        dT->Fill(event.label[hit_i], int(event.time[hit_i]-event.time[trigger_i]));
      }
      break;
    }
  }
  print("hdT created");
}