// g++ -o caen2root caen2root.cpp -Wall -Wextra $(root-config --cflags) $(root-config --glibs) -O2 -std=c++17

#include "AnalysisLib/CFD.hpp"
#include "CaenLib/utils.hpp"
#include "CaenLib/RootEvent.hpp"
#include "CaenLib/CaenRootInterface.hpp"
#include "CaenLib/CaenRootEventBuilder.hpp"
#include "CaenLib/RootHit.hpp"
#include "LibCo/Classes/Timer.hpp"
#include "LibCo/Classes/Timeshifts.hpp"
#include "LibCo/libCo.hpp"


constexpr int reader_version = 110;
constexpr size_t LUT_size = 10000;

#ifdef TRIGGER
  #include "Triggers/TriggerHub.hpp"
#endif //TRIGGER

using namespace Colib;

/**
 * Structure of the output tree :
 * 
    UInt_t    label         : Global label (=board_ID*16 + channel_ID*2 + subchannel_ID)
    UShort_t  board_ID      : Board label [0;max board ID]
    UShort_t  channel_ID    : Channel label [0;8]
    UShort_t  subchannel_ID : Sub channel label [0,1]
    UShort_t  adc           : PHA : ADC value.               PSD : qshort value.
    UShort_t  qlong         : PHA : EXTRA[16:32] (see doc.). PSD : qlong value.
    ULong64_t extended_ts   : Raw timestamp. Units : ps. Used only if fine timestamp and extended timestamp are found in the data. Internal use only.
    ULong64_t timestamp     : Raw timestamp (ts). Units : ps. Is equal to precise_ts if fine ts is found in the data, or extended_ts if only extended ts is found, or TRIGGER_TIME_TAG if no extended timestamp found
    ULong64_t precise_ts    : Raw timestamp. Units : ps. Used only if fine timestamp mode is found in the data. Internal use only.
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
    {0, 5},
    {1, 5},
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

  auto useCFD = LUT<LUT_size>([&](int boardID)
  {
    return key_found(CFD::sShifts, boardID);
  });

  Timer timer;

  // Parameters :
  std::string outpath = "./";
  std::vector<bool> boardReadTrace(100, true);
  size_t nbHitsMax = -1;
  size_t nbHitsMaxTot = -1;
  bool readTraces = true;
  bool writeTraces = false;
  bool ts_evt_build = false;
  bool group = true;
  bool inMemory = true;
  
  std::vector<std::string> filenames;
  Timeshifts timeshifts;

  uint64_t time_window = 2e6 ; // ps

  std::vector<UShort_t> trigger_labels;

  auto printHelp = [](){
    print("caen2root usage");
    print("Note: if accepted, it means that e.g. 1e3 is a valid shorthand for 1000)");
    print("-e --ts-evt-build      [0 or 1] (default 0). Perform event building based on the raw timestamp instead of the absolute time (usually corrected by cfd).");
    print("-f --files             [caen_filename] : File to convert. Include wildcards * and ?, but ONLY IF the name is guarded by quotes (i.e. -f \"/path/to/file/names*.caendat\") ");
    print("-F --files-nb          [caen_filename] [nb_files] : Same as -f. Additionally, can select the number of files (-1 = all, scientific format accepted]");
    print("-g --group             [0 or 1] (default 1) : Sets the output format. 0 : plain tree with additionnal event number and multiplicity fields. 1 : each leaf is a vector.");
    print("-h --help              : print this help");
    print("-i --in-memory         [0 or 1] (default 1) : Choose weither the tree is built in memory (faster but RAM-consuming) or in file (may be much slower)");
    print("-n --hits-nb           [nb_hits (-1 = all, scientific format accepted)] : maximum number of hits to be read by the programm IN EACH FILE");
    print("-N --hits-tot-nb       [nb_hits (-1 = all, scientific format accepted)] : maximum number of hits to be read by the programm (if multithread : not thread safe)");
    print("-o --output            [output_path]");
    print("   --read-traces       [0 or 1] (default 1) : Read the traces for all boards.");
    print("   --board-skip-trace  [boardID] : Skip the trace for this boad ID. Used only if --read-traces is true.");
    print("   --write-traces      [0 or 1] (default 0) : Write the trace in the output root file.");
    print("-T --timeshifts        [filename] : List of timestamp shifts. Format : in each line : global_label timeshift[ns].");
    print("-t --trigger [option] : trigger on the given global label or board ID, or a user-defined file with a list of labels. Example : -t -l 0 to trigger on label 0");
    print("            -l --label  [global_label(16 x boardID + channelID)]]");
    print("            -L --labels [nbLabels] [[global_label(16 x boardID + channelID)]]");
    print("            -b --board  [boardID]] ");
    print("            -B --boards [nbLabels] [[boardID]] ");
    print("            -f --file   [filename (containing a list of labels)]]");
    print("-tw --time-window      [time_window] (default 2000 ns): Coincidence time window for event building (scientific format accepted) ");
    print();
  };
  if (argc < 3) {printHelp(); return 1;}
  std::istringstream iss(argv_to_string(argv));
  std::string temp;
  iss >> temp; // Skipping the first parameter because it is the name of the executable
  while(iss >> temp)
  {
          if (temp == "-f" || temp ==  "--files")
    {
      iss >> temp;
      for (auto const & file : Colib::findFilesWildcard(temp)) {filenames.push_back(file);}
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
    else if (temp == "-i" || temp == "--in-memory")
    {
      iss >> inMemory;
      print("tree built", (inMemory) ? "in memory" : "in file");
    }
    else if (temp == "-n")
    {
      double tmp_d = 0; iss >> tmp_d;
      nbHitsMax = static_cast<size_t>(tmp_d);
      print("Maximum hits per file =", nicer_double(nbHitsMax));
    }
    else if (temp == "-N")
    {
      double tmp_d = 0; iss >> tmp_d;
      nbHitsMaxTot = static_cast<size_t>(tmp_d);
      print("Maximum hits total", nicer_double(nbHitsMaxTot));
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
      print(timeshifts.size(), "timeshifts loaded from", temp);
    }
    else if (temp == "--read-traces")
    {
      iss >> readTraces;
    }
    else if (temp == "--board-skip-trace")
    {
      int boardID; iss >> boardID;
      boardReadTrace[boardID] = false;
    }
    else if (temp == "--write-traces")
    {
      iss >> writeTraces;
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
      if (temp == "-L" || temp == "--labels")
      {
        int nbLabels = 0; iss >> nbLabels;
        int label = 0;
        for (int label_i = 0; label_i < nbLabels; ++label_i)
        {
          iss >> label;
          trigger_labels.push_back(label);
        }
      }
      else if (temp == "-b" || temp == "--board")
      {
        int boardID = 0; iss >> boardID;
        for (int label = boardID*16; label < (boardID+1)*16; ++label) trigger_labels.push_back(label);
      }
      else if (temp == "-B" || temp == "--boards")
      {
        int nbLabels = 0; iss >> nbLabels;
        int boardID = 0;
        for (int label_i = 0; label_i < nbLabels; ++label_i)
        {
          iss >> boardID;
          for (int label = boardID*16; label < (boardID+1)*16; ++label) trigger_labels.push_back(label);
        }
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
    }
    else if (temp == "-tw" || temp == "time-window")
    {
      double e; iss >> e;
      time_window = e*1000.;
    }
    else
    {
      throw_error("Unkown " + temp + " command");
    }
  }

  // -------------------- //
  // Check the parameters //
  // -------------------- //

  const bool hitsMaxSet    = 0 < nbHitsMax   ;
  const bool hitsMaxTotSet = 0 < nbHitsMaxTot;
  const bool useTimeShifts = !timeshifts    .empty();
  const bool trigger_label = !trigger_labels.empty();

  // if (hitsMaxTotSet && multithreadSet) throw-error("Can't have -N and -M at the same time !");
  if (trigger_label) print(trigger_labels);
#ifdef TRIGGER
  if (trigger_label) print(STRINGIFY(TRIGGER), "chosen along with label-based trigger. The behavior is AND (contact dev if you want to switch).");
#endif //TRIGGER

  // Look-up tables (LUT) :
  auto triggerLUT = Colib::LUT<LUT_size>([&trigger_labels](UShort_t label){
    return Colib::found(trigger_labels, label);
  });

  if (filenames.empty()) throw_error("No files !!");
  
  print(filenames.size(), " file"+std::string((filenames.size() > 1) ? "s" : ""));

  size_t nbHitsTot = 0;

  for (auto const & filename : filenames)
  {
    File file(filename);

    if (!file) {error("can't find ", file); continue;}
    
    Caen1725RootInterface reader(filename, readTraces);
    reader.setBoardReadTrace(boardReadTrace);
    Caen1725EventBuilder eventBuilder(reserved_buffer_size);
    eventBuilder.buildOnTimestamp(ts_evt_build);
  #ifdef TRIGGER
    Trigger trigger(&eventBuilder);
  #endif //TRIGGER

    auto rootFilename = outpath + file.shortName()+".root";
    auto rootFile = TFile::Open(rootFilename.c_str(), "recreate");
    if (inMemory) gROOT->cd();
    TString treeName = "HIL";
    if (!group) treeName = "HILplain";
    auto tree = new TTree(treeName, ("WarsawReader_v"+std::to_string(reader_version)).c_str());

    auto & inHit = reader.getHit(); // Aliasing the internal hit of the reader

    Caen1725RootEvent outEvent(writeTraces);
    Caen1725RootHit   outHit  (writeTraces);
    size_t  evtNb   = 0;
    int     evtMult = 0;

    if (group) outEvent.writeTo(tree);
    else
    {
      tree -> Branch("evtNb", &evtNb  );
      tree -> Branch("mult" , &evtMult);
      outHit.writeTo(tree);
    }

    Timer timerRead;
    Timer timerCFD;
    Timer timerTShift;
    Timer timerCopy;
    Timer timerTrigger;
    Timer timerEvtBuild;
    Timer timerFill;

    // To investigate : might not be optimized
    // Pre-declaration of this piece of code
    auto fillTreeAndClear = [&]() -> void
    {
      // 3. Perform the event building 
        timerEvtBuild.Start(); // For profiling
      eventBuilder.fast_event_building(time_window);
        timerEvtBuild.Stop(); // For profiling

      bool triggerBool = true;
      // 4. Write the events to the ROOT tree
      for (auto const & event : eventBuilder)// Loop over all the event in buffer :
      {
        if (!group) evtMult = static_cast<int>(event.size());

          timerTrigger.Start(); // For profiling
        // 4.1 Apply the trigger
        if (trigger_label)
        {
          triggerBool = false;
          for (auto const & hit_i : event) if (triggerLUT[eventBuilder[hit_i].label]) triggerBool = true;
        }
        
      #ifdef TRIGGER
        triggerBool = triggerBool && trigger(event);
      #endif //TRIGGER

          timerTrigger.Stop(); // For profiling

        if (triggerBool) 
        {
          if (group)
          {
            // 4.2 Write the event
              timerCopy.Start(); // For profiling
            for (auto const & hit_i : event)
            {
              outEvent.push_back(eventBuilder[hit_i]);
            }
              timerCopy.Stop(); // For profiling

              timerFill.Start(); // For profiling
            tree -> Fill();
              timerFill.Stop(); // For profiling

            outEvent.clear();
            ++outEvent.evtNb;
          }
          else for (auto const & hit_i : event)
          {
            eventBuilder[hit_i].rel_time = eventBuilder[hit_i].time - eventBuilder[0].time;

              timerCopy.Start(); // For profiling
            outHit = std::move(eventBuilder[hit_i]);
              timerCopy.Stop(); // For profiling
              
              timerFill.Start(); // For profiling
            tree -> Fill();
              timerFill.Stop(); // For profiling
            
            ++evtNb;
          }
        }
      }

      // 5. Reset event builder

      eventBuilder.clear();
    };


    std::vector<int> nbNoZero  (1000,0);
    std::vector<int> nbNoSignal(1000,0);
    std::vector<int> nb        (1000,0);

    // while(reader.readHit())
    while(true)
    {
      // 0. Read the data 
        timerRead.Start(); // For profiling
      auto breaking = !reader.readHit();
        timerRead.Stop(); // For profiling
      
      // 0.1 Manage the end of the reading (end of file or max number of hits reached)
      if (breaking) break;
      if ((hitsMaxSet && nbHitsMax < reader.nbHits()) || (hitsMaxTotSet && nbHitsMaxTot < nbHitsTot)) break;

      ++nbHitsTot;

      // 0.2 Print the hits numbers
      if (reader.nbHits() > 0 && reader.nbHits() % int(1e5) == 0) 
      {
        if (group) printsln(nicer_double(reader.nbHits(), 1), "hits in .caendat", nicer_double(tree->GetEntries(), 1), "evts in .root       ");
        else       printsln(nicer_double(reader.nbHits(), 1), "hits in .caendat", nicer_double(tree->GetEntries(), 1), "hits in .root       ");
      }

      // 1. Apply the cfd
        timerCFD.Start(); // For profiling
      if (applyCFD && readTraces && inHit.hasTrace() && useCFD[inHit.board_ID])
      {
        ++nb[inHit.label];
        CFD cfd(*inHit.getTrace(), CFD::sShifts[inHit.board_ID], CFD::sFractions[inHit.board_ID], 10);
        auto zero = cfd.findZero();
        // auto zero = cfd.findZero(CFD::sThresholds[inHit.board_ID]); // Not using thresholds anymore
             if (zero==CFD::noSignal) {inHit.time = inHit.precise_ts; ++nbNoSignal[inHit.label];}
        else if (zero==CFD::noZero  ) {inHit.time = inHit.precise_ts; ++nbNoZero  [inHit.label];}
        else                           inHit.time = inHit.extended_ts + zero * Caen1725::ticks_to_ps;
      }
      else inHit.time = inHit.precise_ts;
        timerCFD.Stop(); // For profiling

      //1.1 Apply the time shifts if registered
        timerTShift.Start(); // For profiling
      if (useTimeShifts && static_cast<size_t>(inHit.label) < timeshifts.size()) 
      {
        inHit.timestamp += timeshifts[inHit.label];
        inHit.time      += timeshifts[inHit.label];
      }
        timerTShift.Stop(); // For profiling

      // 2. Fill the event builder buffer
      
      if (!eventBuilder.fill_buffer(std::move(inHit))) fillTreeAndClear(); // If buffer is full, fill the tree and clear the buffer
    }
    
    print();
    print("nbNoSignal");
    for (size_t label_i = 0; label_i<nbNoSignal.size(); ++label_i) if (nbNoSignal[label_i]>0) print(label_i, nbNoSignal[label_i], 100.*nbNoSignal[label_i]/double(nb[label_i]), "%");
    print("nbNoZero");
    for (size_t label_i = 0; label_i<nbNoZero  .size(); ++label_i) if (nbNoZero  [label_i]>0) print(label_i, nbNoZero  [label_i], 100.*nbNoZero  [label_i]/double(nb[label_i]), "%");

    fillTreeAndClear(); // Fill the tree with the last event

    rootFile->cd();

    print();
    if (group) print(nicer_double(tree->GetEntries()), "evts in .root       ");
    else       print(nicer_double(tree->GetEntries()), "hits in .root       ");
      Timer timerWrite; // For profiling
    tree->Write();
      timerWrite.Stop(); // For profiling
    rootFile->Close();
    
    if (true)
    {
      print("timeRead", timerRead.timeElapsed());
      print("timeCFD", timerCFD.timeElapsed());
      print("timeTShift", timerTShift.timeElapsed());
      print("timeEvtBuild", timerEvtBuild.timeElapsed());
      print("timeCopy", timerCopy.timeElapsed());
      print("timerTrigger", timerTrigger.timeElapsed());
      print("timeFill", timerFill.timeElapsed());
      print("timeWrite", timerWrite.timeElapsed());
    }

    print(rootFile->GetName(), "written");
  }

  print(timer());
  return 0;
}