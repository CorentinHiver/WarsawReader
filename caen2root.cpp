// g++ -o caen2root caen2root.cpp -Wall -Wextra $(root-config --cflags) $(root-config --glibs) -O2 -std=c++17

#include "AnalysisLib/CFD.hpp"
#include "CaenLib/utils.hpp"
#include "CaenLib/RootEvent.hpp"
#include "CaenLib/CaenRootReader.hpp"
#include "CaenLib/CaenRootEventBuilder.hpp"
#include "CaenLib/RootHit.hpp"
#include "LibCo/Classes/Timer.hpp"
#include "LibCo/Classes/Timeshifts.hpp"
#include "LibCo/libCo.hpp"

constexpr int reader_version = 110;

using namespace Colib;

/**
 * Structure of the output tree :
 * 
    UInt_t    label         : Global label (=board_ID*16 + channel_ID*2 + subchannel_ID)
    UShort_t  board_ID      : Board label [0;max board ID]
    UShort_t  channel_ID    : Channel label [0;8]
    UShort_t  subchannel_ID : Sub channel label [0,1]
    Int_t     adc           : PHA : ADC value. PSD : qshort value.
    Int_t     qlong         : PHA : unused.    PSD : qlong value.
    ULong64_t timestamp     : Raw timestamp (ts). Units : ps. Is equal to precise_ts if fine ts is found in the data, or extended_ts if only extended ts is found, or TRIGGER_TIME_TAG if no extended timestamp found
    ULong64_t extended_ts   : Raw timestamp. Units : ps. Used only if fine timestamp and extended timestamp are found in the data.
    ULong64_t precise_ts    : Raw timestamp. Units : ps. Used only if fine timestamp mode is found in the data.
    ULong64_t time          : Absolute time. Units : ps. This field must be filled by the user (i.e., the program using this class).
    Int_t     rel_time      : Relative time in the event. Units : ps. This field must be filled by the user (i.e., the program using this class).
 */

//////////////////////////////////////
// Version 1.03                     //
// Is now included :                //
//    - Writting events in .root    //
//  Todo :                          //
//    - Include traces in events    //
//////////////////////////////////////

//////////////////////////////////////
// Version 1.02                     //
// Is now included :                //
//    - Label-base trigger logic    //
//    - Wildcard file gathering     //
//////////////////////////////////////

/* Do not read the following, this trigger logic is stil TBD :
How to use the trigger : 
- First, create your own trigger in the folder Triggers by following the instructions
- Then, include your trigger in the compilation command found at the end of this file : 
  g++ (the rest of the line) -DTRIGGER="\"Triggers/TheNameOfYourTriggerFile.hpp\""
  
  Attention : do not forget the ' \" ', if you know how to get rid of this requirement contact me
*/

constexpr size_t   reserved_buffer_size = 500000ul;
constexpr bool applyCFD = true;

///////////////////////////////////////////////////////////////
// So far, the cfd parameters are set board by board.        //
// An improvement would be to apply it detector by detector, //
// and most importantly, to read a parameter file rather     //
// that these hard-coded values                              //
///////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  CFD::sShifts = { // BOARD_ID, SAMPLES
    {0, 7},
    {1, 7},
    {6, 2},
    {7, 2},
    {8, 2}
  };

  CFD::sFractions = { // BOARD_ID, fraction
    {0, 0.75},
    {1, 0.75},
    {6, 0.75},
    {7, 0.75},
    {8, 0.75}
  };

  auto useCFD = LUT<9*16>([&](int boardID)
  {
    return key_found(CFD::sShifts, boardID);
  });

  Timer timer;

  // Parameters :
  std::string outpath = "./";
  bool analyseTraces = true;
  bool storeTraces = false;
  size_t nbHitsMax = -1;
  bool hitsMaxSet = false;
  bool ts_evt_build = false;
  bool group = true;
  
  std::vector<std::string> filenames;
  Timeshifts timeshifts;

  uint64_t time_window = 2e6 ; // ps

  std::vector<u_int> trigger_labels;
  std::vector<u_int> trigger_boards;

  auto printHelp = [](){ 
    print("caen2root usage");
    print("-e --ts-evt-build      [0 or 1] (default 0). Perform event building based on the raw timestamp instead of the absolute time (usually corrected by cfd).");
    print("-f --files             [caen file name] : File to convert. Include wildcards * and ?, but ONLY IF the name is guarded by quotes (i.e. -f \"/path/to/file/names*.caendat\") ");
    print("-F --files-nb          [caen file name] [nb_files] : Files to convert. Include wildcards * and ?, but ONLY IF the name is guarded by quotes (i.e. -f \"/path/to/file/names*.caendat\"). For nb_files, -1 = all, scientific format accepted (e.g. 1.e3 (=1000))]");
    print("-g --group             [0 or 1] (default 1) : Sets the output format. 0 : plain tree with additionnal event number and multiplicity fields. 1 : each leaf is a vector.");
    print("-h --help              print this help");
    print("-n                     [number of hits (-1 = all, 1.e3 (=1000) format accepted)]");
    print("-o --output            [output path]");
    print("   --trace-analysis    [0 or 1] (default 1). Perform trace analysis (so far, only cfd is implemented).");
    print("   --trace-storing     [0 or 1] (default 0). Store trace in the root tree.");
    print("-T --timeshifts        [filename] : List of timestamp shifts. Format : in each line : global_label timeshift[ns].");
    print("-t --trigger :");
    print("            -l --label [global_label(16 x boardID + channelID)]]");
    print("            -b --board [boardID]] ");
    print("            -f --file  [filename (containing a list of global_labels)]]");
    print("-tw --time-window      [time_window (float, in ns)] (default 2000 ns) ");
    print();
    print("example of a command line including all the above options :");
    print();
    print("./caen2root \\");
    print("\t -f /data/experiment1/run0/*.caendat \\");
    print("\t -F /data/experiment1/run1/*.caendat 3 \\");
    print("\t -n 1e9 \\");
    print("\t -o /root_data/experiment1/ \\");
    print("\t -t --label 2 \\");
    print("\t -t --board 5 -t --board 6 \\");
    print("\t trace-analysis 1 \\");
    print("\t trace-storing 1 \\");
    print();
    print("This made the folllowing :");
    print("Gather all the .caendat files found in /data/experiment1/run0/ and the three first found in the run0/ folder.");
    print(" The code will stop after processing 1e9 hits. The output root files will be written in the /root_data/experiment1/ folder.");
    print(" Only event containing at least one detector with label 5, or any detector of board_ID==5 or 6 are kept.");
    print(" The traces will be kept for analysis and written in the root file");
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
        for (auto const & file : Colib::findFilesWildcard(temp)){ print(file); filenames.push_back(file);}
      }
      else if (temp == "-F" || temp == "--files-nb")
      {
        iss >> temp;
        double nb = -1; iss >> nb;
        size_t nb_i = static_cast<size_t>(nb); // cast to size_t (). if nb<0 -> overflow -> nb_i is very very big
        if (nb_i == 0) continue;
        auto const & files = Colib::findFilesWildcard(temp);
        for (size_t file_i = 0; file_i < (files.size() && file_i<nb_i); ++file_i) filenames.push_back(files[file_i]);
      }
      else if (temp == "-h" || temp == "--help")
      {
        printHelp();
      }
      else if (temp == "-n")
      {
        double tmp_d = 0; iss >> tmp_d;
        nbHitsMax = tmp_d;
        hitsMaxSet = true;
      }
      else if (temp == "-o")
      {
        iss >> outpath;
      }
      else if (temp == "-g" || temp == "--group")
      {
        iss >> group;
      }
      else if (temp == "-T" || temp == "--timeshifts")
      {
        iss >> temp;
        timeshifts.load(temp);
        // timeshifts*=1000;
      }
      else if (temp == "--trace-analysis")
      {
        iss >> analyseTraces;
      }
      else if (temp == "--trace-storing")
      {
        iss >> storeTraces;
      }
      else if (temp == "--ts-evt-build")
      {
        iss >> ts_evt_build;
      }
      else if (temp == "-t" || temp == "--trigger")
      {
        iss >> temp;
        if (temp == "-l" || temp == "--label")
        {
          int label = 0;
          iss >> label;
          trigger_labels.push_back(label);
        }
        else if (temp == "-b" || temp == "--board")
        {
          int boardID = 0; iss >> boardID;
          for (int label = boardID*16; label < (boardID+1)*16; ++label) trigger_labels.push_back(label);
        }
        else if (temp == "-f" || temp == "--file")
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
        else if (temp == "-tw" || temp == "time-window")
        {
          double e; iss >> e;
          time_window = e*1000;
        }
      }
    }

    if (filenames.empty()) throw_error("No files !!");

    for (auto const & filename : filenames)
    {
      File file(filename);
      if (!file) {error("can't find ", file); continue;}
      
      CaenRootReader1725 reader(filename, analyseTraces);
      CaenRootEventBuilder1725 eventBuilder(reserved_buffer_size);
      eventBuilder.buildOnTimestamp(ts_evt_build);

      auto rootFilename = outpath + file.shortName()+".root";
      auto rootFile = TFile::Open(rootFilename.c_str(), "recreate");
      TString treeName = "HIL";
      if (!group) treeName = "HILplain";
      auto tree = new TTree(treeName, ("WarsawReader_v"+std::to_string(reader_version)).c_str());

      auto & inHit = reader.getHit(); // Aliasing the internal hit of the reader

      RootCaenEvent outEvent(storeTraces);
      RootCaenHit   outHit  (storeTraces);
      size_t  evtNb   = 0;
      int     evtMult = 0;

      if (group) outEvent.writeTo(tree);
      else
      {
        tree -> Branch("evtNb", &evtNb  );
        tree -> Branch("mult" , &evtMult);
        outHit.writeTo(tree);
      }

      double timeRead = 0;
      double timeCFD = 0;
      double timeTShift = 0;
      double timeEvtBuild = 0;
      double timeCopy = 0;
      double timeFill = 0;
      double timeWrite = 0;

      // To investigate : might not be optimized
      // Pre-declaration of this piece of code
      auto fillTree = [&]() -> void
      {
        // 3. Perform the event building 
          Timer timerEvtBuild;
        eventBuilder.fast_event_building(time_window);
          timeEvtBuild += timerEvtBuild.Time();

        // 4. Write the events to the ROOT tree
        for (auto const & event : eventBuilder)// Loop over all the event in buffer :
        {
          evtMult = event.size();

        #ifdef TRIGGER
          // TODO:
        #endif //TRIGGER

          // 4.1 Apply the trigger
          static thread_local bool trigger; trigger = true;
          static thread_local bool trigger_label = !trigger_labels.empty();
          if (trigger_label) 
          {
            trigger = false;
            for (auto const & hit_i : event) if (Colib::found(trigger_labels, eventBuilder[hit_i].label)) trigger = true;
          }
          
          if (trigger) 
          {
            // 4.2 Write the event
            for (auto const & hit_i : event)
            {
              // Calculate the relative time
              if (group)
              {
                  Timer timerCopy;
                outEvent.push_back(eventBuilder[hit_i]);
                  timeCopy += timerCopy.Time();
              }
              else
              {
                eventBuilder[hit_i].rel_time = eventBuilder[hit_i].time - eventBuilder[0].time;
                  Timer timerCopy;
                outHit.copy(eventBuilder[hit_i], false);
                  timeCopy += timerCopy.Time();
                  Timer timerFill;
                tree -> Fill();
                  timeFill += timerFill.Time();
                ++evtNb;
              }
            }
            if (group)
            {
              // Colib::printPause(outEvent);
                Timer timerFill;
              tree -> Fill();
                timeFill += timerFill.Time();
              outEvent.clear();
              ++outEvent.evtNb;
            }
          }
        }

        // 5. Reset event builder

        eventBuilder.clear();
      };

      while(true)
      // while(reader.readHit())
      {
          Timer timerRead;
        if (!reader.readHit()) break;
          timeRead += timerRead.Time();
        if (hitsMaxSet && nbHitsMax < reader.nbHits()) break;
        if (reader.nbHits() > 0 && reader.nbHits() % int(1e5) == 0) printsln(nicer_double(reader.nbHits(), 1), "    ");

        // 1. Apply the cfd
          Timer timerCFD;
        if (applyCFD && inHit.getTrace() && !inHit.getTrace() -> empty() && useCFD[inHit.board_ID])
        {
          CFD cfd(*inHit.getTrace(), CFD::sShifts[inHit.board_ID], CFD::sFractions[inHit.board_ID]);
          auto zero = cfd.findZero(CFD::sThresholds[inHit.board_ID]);
               if (zero==CFD::noSignal) {inHit.time = inHit.precise_ts; debug("noSignal");}
          else if (zero==CFD::noZero  ) {inHit.time = inHit.precise_ts; debug("noZero");}
          else                          inHit.time = inHit.extended_ts + zero * CaenDataReader1725::ticks_to_ps;
        }
        else inHit.time = inHit.precise_ts;
          timeCFD += timerCFD.Time();

        //1.1 Apply the time shifts if registered
          Timer timerTShift;
        if (static_cast<size_t>(inHit.label) < timeshifts.size()) 
        {
          inHit.timestamp += timeshifts[inHit.label];
          inHit.time      += timeshifts[inHit.label];
        }
          timeTShift += timerTShift.Time();

        // 2. Fill the event builder buffer
        
          Timer timerCopy;
        if (eventBuilder.fill_buffer(inHit)) 
        {
            timeCopy += timerCopy.Time();
          continue; // Continue the loop as long as the buffer is not filled
        }
        // (this piece of code is reached only when the buffer is full)

        fillTree();
      }

      fillTree(); // Fill the tree with the last event

      rootFile->cd();

      print();
      if (group) print(tree->GetEntries(), "events in the tree");
      else       print(tree->GetEntries(), "hits in the tree");
        Timer timerFill;
      tree->Write();
        timeFill += timerFill.Time();
      rootFile->Close();
      
      if (false)
      {
        print("timeRead", timeRead/1000, "ms");
        print("timeCFD", timeCFD/1000, "ms");
        print("timeEvtBuild", timeEvtBuild/1000, "ms");
        print("timeCopy", timeCopy/1000, "ms");
        print("timeFill", timeFill/1000, "ms");
        print("timeWrite", timeWrite/1000, "ms");
      }

      print(rootFile->GetName(), "written");
    }
    
  }
  print(timer());
  return 0;
}