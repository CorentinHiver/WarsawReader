#include "../CaenLib/RootReader.hpp"
#include "../Colib/lib/Classes/Timeshifts.hpp"

using namespace Colib;
using namespace Caen1725;
using namespace std;

void dTcalculate(std::string filename = "", int ref_label = 0)
{
  int maxLabel = 200;
  RootReader reader(filename);
  auto & event = reader.getEvent();
  auto dT = new TH2F("hdT","dT;label;dT[ps]",maxLabel,0,maxLabel, 6000,-3000000,3000000);
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
  map<int, double> means; 
  for (int label = 0; label<maxLabel; ++label)
  {
    string label_str = to_string(label);
    auto histo = dT->ProjectionY(label_str.c_str(), label+1, label+1);
    auto const nbHits = histo->GetEntries();
    if (nbHits < 1) continue;
    if (nbHits < int(1e3)) {error("For label", label, "there are only", nbHits); continue;}
    
    auto dTmaxY = histo->GetMaximum();
    auto dTmaxX = histo->GetMaximumBin();
    auto leftPeak = histo->FindFirstBinAbove(dTmaxY*0.7);
    auto rightPeak = histo->FindLastBinAbove(dTmaxY*0.7);

    histo->GetXaxis()->SetRange(leftPeak * 0.9, rightPeak * 1.1);
    
    // dTmaxX = histo->GetMaximumBin();
    leftPeak  = histo->FindFirstBinAbove(dTmaxY*0.5);
    rightPeak = histo->FindLastBinAbove (dTmaxY*0.5);
    auto proto_sigma = rightPeak-leftPeak;

    histo->GetXaxis()->SetRange(dTmaxX-proto_sigma, dTmaxX+proto_sigma);

    means.emplace(label, histo -> GetMean());
  }
  ofstream dTfile("unname.dT");
  for (auto const & [label, mean] : means) dTfile << label << " " << mean << "\n";
  dTfile.close();
}