// The point of this code is to develop CaenReaderLib

#include "CaenLib/CaenRootReader.hpp"
#include "TFile.h"

void readRaw()
{
  std::vector<std::string> filenames = {
    "/home/corentin/60Co_data/eagleRU_i2514_0023_0000.caendat"
  };
  for (auto const & filename : filenames)
  {
    TTree* outTree = new TTree("caenData", "caenData");
    CaenRootReader reader(filename, outTree);
    int i = 0;
    while(reader.readHit() && ++i<int(1e5)) outTree->Fill();

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
// g++ -o exec readRaw.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2
// g++ -o exec readRaw.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -g -DDEBUG

// Traces : 
// g++ -o exec readRaw.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2 -DTRACES
// g++ -o exec readRaw.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -g -DDEBUG -DTRACES