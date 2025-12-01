#include "CaenLib/utils.hpp"
#include "AnalysisLib/CFD.hpp"
#include "CaenLib/CaenRootReader.hpp"
#include "CaenLib/CaenRootEventBuilder.hpp"
#include "LibCo/FilesManager.hpp"
#include "LibCo/Timer.hpp"
#include "LibCo/libCo.hpp"

constexpr int reader_version = 1;

using namespace Colib;

//////////////////////////////////////
// Version 1.02                     //
// Is now included :                //
//    - Label-base trigger logic    //
//    - Wildcard file gathering     //
//                                  //
//////////////////////////////////////

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
  size_t nbHitsMax = -1;
  bool hitsMaxSet = false;

  std::vector<int> trigger_labels;
  std::vector<int> trigger_boards;

  auto printHelp = [](){ 
    print("-f --files            [caen file name (include wildcards * and ?)] ");
    print("-F --files_and_number [caen file name (include wildcards * and ?)] [nb_files (-1 = all, 1.e3 (=1000) format accepted)]");
    print("-h --help");
    print("-n                    [number of hits (-1 = all, 1.e3 (=1000) format accepted)]");
    print("-o --output           [output path]");
    print("-t --trigger          [--label [global_label(16 x boardID + channelID)]] [--board [boardID]] [--file [filename]]");
    print("   --no-traces");
  };
  if (argc < 3) {printHelp(); return 1;}
  else
  {
    std::istringstream iss(argv_to_string(argv));
    std::string temp;
    while(iss >> temp)
    {
           if (temp == "-f" || temp ==  "--files")
      {
        iss >> temp;
        for (auto const & file : Colib::findFilesWildcard(temp)) filenames.push_back(file);
      }
      else if (temp == "-F" || temp == "--files_and_number")
      {
        iss >> temp;
        double nb = -1; iss >> nb;
        size_t nb_i = nb; // cast to size_t ()
        if (nb_i == 0) continue;
        auto const & files = Colib::findFilesWildcard(temp);
        for (size_t file_i = 0; file_i<files.size() && file_i<nb_i; ++file_i) filenames.push_back(files[file_i]);
      }
      else if (temp == "-h" || temp == "--help")
      {
        printHelp();
      }
      else if (temp == "-n")
      {
        double tmp_d; iss >> tmp_d;
        nbHitsMax = tmp_d;
        hitsMaxSet = true;
      }
      else if (temp == "-o")
      {
        iss >> outpath;
      }
      else if (temp == "--no-traces")
      {
        handleTraces = false;
      }
      else if (temp == "-t" || temp == "--trigger")
      {
        iss >> temp;
        if (temp == "--label")
        {
          int label;
          iss >> label;
          trigger_labels.push_back(label);
        }
        else if (temp == "--board")
        {
          int boardID;
          iss >> boardID;
          for (int label = boardID*16; label < (boardID+1)*16; ++label) trigger_labels.push_back(label);
        }
        else if (temp == "--file")
        {
          iss >> temp;
          std::ifstream trigger_file(temp);
          if (!trigger_file.is_open()) throw_error("Can't open trigger file "+temp+" !!");
          std::string line;
          while(std::getline(trigger_file, line))
          {
            std::istringstream iss(line);
            int label;
            while(iss >> label) trigger_labels.push_back(label);
          }
        }
      }
    }

    if (filenames.empty()) throw_error("No files !!");

    for (auto const & filename : filenames)
    {
      File file(filename);
      if (!file) {error("can't find ", file); continue;}
      
      CaenRootReader1725 reader(filename);
      reader.handleTraces(handleTraces);
      CaenRootEventBuilder1725 eventBuilder(reserved_buffer_size);

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

        // 4. Write the events to the ROOT tree

        for (auto const & event : eventBuilder)
        {
          evtMult = event.size();
          ++evtNb;

        #ifdef TRIGGER
          // TODO:
        #endif //TRIGGER

          static thread_local auto trigger_label = !trigger_labels.empty();
          bool trigger = true;
          if (trigger_label) 
          {
            trigger = false;
            for (auto const & hit_i : event) 
            {
              if (Colib::found(trigger_labels, eventBuilder[hit_i].label)) trigger = true;
            }
          }

          if (trigger) for (auto const & hit_i : event)
          {
            outHit.copy(eventBuilder[hit_i]);
            tree -> Fill();
          }
        }

        // 5. Reset event builder

        eventBuilder.clear();
      };

      while(reader.readHit())
      {
        if (hitsMaxSet && nbHitsMax < reader.nbHits()) break;
        if (reader.nbHits() > 0 && reader.nbHits() % int(1e5) == 0) print(nicer_double(reader.nbHits(), 1));

        // 1. Apply the cfd

        if (applyCFD && !inHit.getTrace().empty() && key_found(CFD::sShifts, inHit.board_ID)) 
        {
          CFD cfd(inHit.getTrace(), CFD::sShifts[inHit.board_ID], CFD::sFractions[inHit.board_ID]);
          
          auto zero = cfd.findZero(CFD::sThresholds[inHit.board_ID]); 
          if (zero == CFD::noSignal) inHit.cfd = inHit.precise_ts;
          else
          {
            zero = zero * CaenDataReader1725::ticks_to_ns; // Convert from 4 ns ticks to ps, !! might change depending on the daq setup !!  
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