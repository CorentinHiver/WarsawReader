// The point of this code is to develop CaenReaderLib

#include "CaenLib/CaenRootReader.hpp"
#include "TFile.h"

void readRaw()
{
  std::vector<std::string> filenames = {
    // "/home/corentin/data/60Co_Easter/eagleRU_i2608_0006_0000.caendat"
    "coulexRU_i3147_3005_0000.caendat"
  };
  for (auto const & filename : filenames)
  {
    TTree* outTree = new TTree("caenData", "caenData");
    Caen1725RootInterface reader(filename, outTree, true);
    while(reader.readHit() && reader.nbHits()<int(2e5)) 
    {
      if (reader.nbHits()%int(1e4) == 0) printsln(reader.nbHits()*100/int(2e5), "%");
      outTree->Fill();
    }

    auto outFile = TFile::Open("readRaw.root", "recreate");
      outTree->Write();
    outFile->Close();
    print("readRaw.root written");
  }
}

int main()
{
  readRaw();
}

// Regular :
// g++ -o readRaw readRaw.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2
// g++ -o readRaw readRaw.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -g -DDEBUG