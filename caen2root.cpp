#include "AnalysisLib/CFD.hpp"
#include "CaenLib/CaenRootReader.hpp"
#include "CaenLib/CaenRootEventBuilder.hpp"
#include "LibCo/FilesManager.hpp"
#include "LibCo/Timer.hpp"

constexpr int reader_version = 1;

//////////////////////////////////////
// Version 1.                       //
// Is now included :                //
//    - Modularized CFD parameters  //
// Todo :                           //
//    - A trigger logic             //
//                                  //
//////////////////////////////////////


/*
How to use the trigger : 
- First, create your own trigger in the folder Triggers by following the instructions
- Then, include your trigger in the compilation command found at the end of this file : 
  g++ (the rest of the line) -DTRIGGER="\"Triggers/TheNameOfYourTriggerFile.hpp\""
  
  Attention : do not forget the ' \" ', if you know how to get rid of this requirement contact me
*/

constexpr Long64_t time_window          = 2e6 ; // ps
constexpr size_t   reserved_buffer_size = 500000ul;
constexpr bool applyCFD = true;

///////////////////////////////////////////////////////////////
// So far, the cfd parameters are set board by board.        //
// An improvement would be to apply it detector by detector, //
// and most importantly, to read a parameter file rather     //
// that these hard-coded values                              //
///////////////////////////////////////////////////////////////

CFD::Shifts cfd_shifts = {
  {0, 5},
  {1, 5},
  {6, 2},
  {7, 2},
  {8, 2}
};

CFD::Thresholds cfd_thresholds = {
  {0, -50},
  {1, -50},
  {6, -500},
  {7, -500},
  {8, -500}
};

CFD::Fractions cfd_fractions = {
  {0, 0.5},
  {1, 0.5},
  {6, 0.5},
  {7, 0.5},
  {8, 0.5}
};

int main(int argc, char** argv)
{
  Timer timer;

  // Parameters :
  std::vector<std::string> filenames;
  std::string outpath = "./";
  bool handleTraces = true;
  Long64_t nbHitsMax = -1;
  auto printHelp = [](){ 
    print("-f [file name]");
    print("-F [folder name]");
    print("-h");
    print("-n [number of hits]");
    print("-o [output path]");
    print("--no-traces");
  };
  if (argc < 3) {printHelp(); return 1;}
  else
  {
    std::istringstream iss(argv_to_string(argv));
    std::string temp;
    while(iss >> temp)
    {
           if (temp == "-f")
      {
        iss >> temp;
        filenames.push_back(temp);
      }
      else if (temp == "-F")
      {
        iss >> temp;
        FilesManager files; files.addFolder(temp, -1, {"caendat"});
        for (auto const & file : files) filenames.push_back(file);
      }
      else if (temp == "-L")
      {
        throw_error("list mode -L not implemented yet");
      }
      else if (temp == "-o")
      {
        iss >> outpath;
      }
      else if (temp == "--no-traces")
      {
        handleTraces = false;
      }
      else if (temp == "-n")
      {
        double tmp_d;
        iss >> tmp_d;
        nbHitsMax = tmp_d;
      }
      else if (temp == "-h")
      {
        printHelp();
      }
    }

    for (auto const & filename : filenames)
    {
      File file(filename);
      if (!file) {error("can't find ", file); continue;}
      
      CaenRootReader reader(filename);
      reader.handleTraces(handleTraces);
      CaenRootEventBuilder eventBuilder(reserved_buffer_size);

      auto rootFilename = outpath + file.shortName()+".root";
      auto rootFile = TFile::Open(rootFilename.c_str(), "recreate");
      auto tree = new TTree("HIL", ("WarsawReader_v"+std::to_string(reader_version)).c_str());

      auto & inHit = reader.getHit();

      int evtNb = 0;
      int evtMult = 0;
      tree -> Branch("evtNb", &evtNb);
      tree -> Branch("evtMult", &evtMult);

      RootCaenHit outHit;
      outHit.writeTo(tree);

      // Pre-declaration of this piece of code
      auto fillTree = [&]()
      {
        // 3. Perform the event building 

        eventBuilder.fast_event_building(time_window);

        // 4. Write the hits to the tree

        for (auto const & event : eventBuilder)
        {
          evtMult = event.size();

        #ifdef TRIGGER
          
        #endif //TRIGGER

          for (auto const & hit_i : event)
          {
            outHit.copy(eventBuilder[hit_i]);
            tree -> Fill();
          }
          ++evtNb;
        }

        // 5. Reset event builder

        eventBuilder.clear();
      };

      while(reader.readHit())
      {
        if (nbHitsMax>0 && reader.nbHits()>nbHitsMax) break;
        if (reader.nbHits() > 0 && reader.nbHits() % int(1e5) == 0) print(nicer_double(reader.nbHits(), 1));

        // 1. Apply the cfd

        if (applyCFD && !inHit.getTrace().empty() && key_found(CFD::sShifts, inHit.board_ID)) 
        {
          CFD cfd(inHit.getTrace(), CFD::sShifts[inHit.board_ID], CFD::sFractions[inHit.board_ID]);
          
          auto zero = cfd.findZero(CFD::sThresholds[inHit.board_ID]); 
          if (zero == CFD::noSignal) inHit.cfd = inHit.precise_ts;
          else
          {
            zero = zero * 4000.; // Convert from 4 ns ticks to ps, !! might change depending on the daq setup !!  
            inHit.cfd = inHit.extended_ts + zero;
          }
        }
        else inHit.cfd = inHit.precise_ts;

        // 2. Fill the event builder buffer

        if (eventBuilder.fill_buffer(inHit)) continue; // Continue the loop as long as the buffer is not filled
        
        // (this piece of code is reached only when the buffer is full)

        fillTree();
      }

      fillTree(); // Fill the tree with the last event

      rootFile->cd();

        print(tree->GetEntries());
        tree->Write();

      rootFile->Close();
      print(rootFile->GetName(), "written");
    }
  }
  print(timer());
  return 0;
}

// g++ -o caen2root caen2root.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2