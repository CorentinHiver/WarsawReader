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
    Int_t     adc           : PHA : ADC value.               PSD : qshort value.
    Int_t     qlong         : PHA : EXTRA[16:32] (see doc.). PSD : qlong value.
    ULong64_t caen_time     : Raw timestamp. Units : ps. Is equal to extended_timestamp if fine or extended timestamp is found in the data, or TRIGGER_TIME_TAG if no extended timestamp found
    ULong64_t time          : Absolute time. Units : ps. This field must be filled by the user (i.e., the program using this class).
    // Int_t     rel_time      : Relative time in the event. Units : ps. This field must be filled by the user (i.e., the program using this class).
 */

//////////////////////////////////////
// Version 1.04                     //
// Validated on real data           //
//////////////////////////////////////

//////////////////////////////////////
// Version 1.03                     //
// Is now included :                //
//    - Writting events in .root    //
//  Todo :                          //
//    - Include traces in events    //
//////////////////////////////////////

/* Do not read the following, this trigger logic is stil TBD :
How to use the trigger : 
- First, create your own trigger in the folder Triggers by following the instructions
- Then, include your trigger in the compilation command found at the end of this file : 
  g++ (the rest of the line) -DTRIGGER="\"Triggers/TheNameOfYourTriggerFile.hpp\""
  
  Attention : do not forget the ' \" ', if you know how to get rid of this requirement contact me


  TODO :
  Ajouter un booléen wfa_success // wave form analysis successful
  renommer timestamp board_time
  renommer time best_time
*/

class Profiler : public Timer
{
public:
  Profiler() noexcept : Timer() {}
  auto StartProfiling()
  {
  #ifdef PROFILE
    return Start();
  #else
    return;
  #endif //PROFILE
  }
  auto StopProfiling()
  {
  #ifdef PROFILE
    return Stop();
  #else
    return;
  #endif //PROFILE
  }
  void printElapsedTime(std::string const & prepend = "", std::string const & append = "")
  {
    print(prepend, timeElapsed(), append);
  }
};


constexpr size_t   reserved_buffer_size = 500000ul;

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
  bool applyCFD = true;
  
  std::vector<std::string> filenames;
  Timeshifts timeshifts;

  uint64_t time_window = 2e6 ; // ps

  std::vector<UShort_t> trigger_labels;

  auto printHelp = [](){
    print("caen2root usage");
    print("Note: if \"scientific format accepted\", it means that e.g. 1e3 is a valid shorthand for 1000)");
    print("   --cfd               [0 or 1] (default 1). Use CFD timestamp correction (hard-coded parameters(shift, fraction, nb samples for baseline...). Parameter file incoming).");
    print("-e --ts-evt-build      [0 or 1] (default 0). Perform event building based on : [0] the absolute time (usually corrected by cfd) [1] the raw timestamp.");
    print("-f --files             [caen_filename] : File to convert. Include wildcards * and ?, but ONLY IF the name is guarded by quotes (i.e. -f \"/path/to/file/names*.caendat\") ");
    print("-F --files-nb          [caen_filename] [nb_files] : Same as -f. Additionally, can select the number of files (-1 = all, scientific format accepted]");
    print("-g --group             [0 or 1] (default 1) : Sets the output format. 0 : plain tree with additionnal event number and multiplicity fields. 1 : each leaf is a vector representing the event.");
    print("-h --help              : print this help");
    print("-i --in-memory         [0 or 1] (default 1) : Choose weither the tree is built in memory (faster but RAM-consuming) or in file (may be much slower)");
    print("-n --hits-nb           [nb_hits (-1 = all, scientific format accepted)] : maximum number of hits to be read by the programm IN EACH FILE");
    print("-N --hits-tot-nb       [nb_hits (-1 = all, scientific format accepted)] : maximum number of hits to be read by the programm");
    print("-o --output            [output_path]");
    print("   --read-traces       [0 or 1] (default 1) : Read the traces for all boards.");
    print("   --board-skip-trace  [boardID] : Do not read the trace for this board ID, hence no trace analysis (e.g. CFD) is done. Used only if --read-traces is true.");
    print("   --write-traces      [0 or 1] (default 0) : Write the trace in the output root file.");
    print("-T --timeshifts        [filename] : List of timeshifts. Format : in each line : global_label timeshift[ns].");
    print("-t --trigger [option] : trigger on the given global label or board ID, or a user-defined file with a list of labels. Example : -t -l 0 to trigger on label 0");
    print("            -l --label  [global_label(16 x boardID + channelID)]]");
    print("            -L --labels [nb labels] [[global_label(16 x boardID + channelID)]]");
    print("            -b --board  [boardID]] ");
    print("            -B --boards [nb labels] [[boardID]] ");
    print("            -f --file   [filename (containing a list of labels)]]");
    print("-tw --time-window      [time_window ns] (default 2000) : Coincidence time window for event building (scientific format accepted) ");
    print();
  };
  if (argc < 3) {printHelp(); return 1;}
  std::istringstream iss(argv_to_string(argv));
  std::string temp;
  iss >> temp; // Skipping the first parameter because it is the name of the executable
  while(iss >> temp)
  {
         if (temp ==  "--cfd")
    {
      iss >> applyCFD;
    }
    else if (temp == "-f" || temp ==  "--files")
    {
      iss >> temp;
      for (auto const & file : Colib::findFilesWildcard(temp)) {filenames.push_back(file);}
    }
    else if (temp == "-F" || temp == "--files-nb")
    {
      iss >> temp;
      double temp_nb = -1; iss >> temp_nb;
      size_t nb_i = static_cast<size_t>(temp_nb); // cast to size_t (). if nb<0 -> overflow -> nb_i is very very big
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
      int temp_nb = 0; 
      for (size_t i = 0; i<timeshifts.size(); ++i) if (0 < timeshifts[i]) ++temp_nb;
      print(temp_nb, "timeshifts loaded from", temp);
    }
    else if (temp == "--read-traces")
    {
      iss >> readTraces;
    }
    else if (temp == "--board-skip-trace")
    {
      int boardID = 0; iss >> boardID;
      boardReadTrace[boardID] = false;
    }
    else if (temp == "--write-traces")
    {
      iss >> writeTraces;
    }
    else if (temp == "-e" || temp == "--ts-evt-build")
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
  if (trigger_label) print(STRINGIFY(TRIGGER), "chosen along with label-based trigger. The behavior is AND (possibility to develop it if you really need a OR logic).");
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

    Caen1725::RootEvent outEvent(writeTraces);
    Caen1725::RootHit   outHit  (writeTraces);
    Long64_t evtNb   = 0;
    Int_t    evtMult = 0;

    if (group) outEvent.writeTo(tree);
    else
    {
      tree -> Branch("evtNb", &evtNb  );
      tree -> Branch("mult" , &evtMult);
      outHit.writeTo(tree);
    }

    Profiler timerRead;
    Profiler timerCFD;
    Profiler timerTShift;
    Profiler timerCopy;
    Profiler timerTrigger;
    Profiler timerEvtBuild;
    Profiler timerFill;

    CFD cfd;

    // To investigate : might not be optimized
    // Pre-declaration of this piece of code
    auto dumpBufferInTree = [&]() -> void
    {
      // 3. Perform the event building 
        timerEvtBuild.StartProfiling();
      eventBuilder.fast_event_building(time_window);
        timerEvtBuild.StopProfiling();

      // 4. Write the events to the ROOT tree
      for (auto const & event : eventBuilder)// Loop over all the event in buffer :
      {
        bool triggerBool = true;
        if (!group) evtMult = static_cast<int>(event.size());

          timerTrigger.StartProfiling();
        // 4.1 Apply the trigger
        if (trigger_label)
        {
          triggerBool = false;
          for (auto const & hit_i : event) if (triggerLUT[eventBuilder[hit_i].label]) triggerBool = true;
        }
        
      #ifdef TRIGGER
        triggerBool = trigger(event) && triggerBool;
      #endif //TRIGGER

          timerTrigger.StopProfiling();

        if (triggerBool) 
        {
          if (group)
          {
            // 4.2 Write the event
              timerCopy.StartProfiling();
            for (auto const & hit_i : event)
            {
              outEvent.push_back(eventBuilder[hit_i]);
            }
              timerCopy.StopProfiling();

              if (outEvent.mult == Caen1725::Event::maxEvt-1)
              {
                error("Event with multiplicity > 1000 for event n",evtNb,", not written");
                continue;
              }

              timerFill.StartProfiling();
            tree -> Fill();
              timerFill.StopProfiling();

            outEvent.clear();
            ++outEvent.evtNb;
          }
          else for (auto const & hit_i : event)
          {
            eventBuilder[hit_i].rel_time = eventBuilder[hit_i].time - eventBuilder[0].time;

              timerCopy.StartProfiling();
            outHit = std::move(eventBuilder[hit_i]);
              timerCopy.StopProfiling();
              
              timerFill.StartProfiling();
            tree -> Fill();
              timerFill.StopProfiling();
            
            ++evtNb;
          }
        }
      }

      // 5. Reset event builder

      eventBuilder.clear();
    };


    constexpr int labelMax = 1000;
    int nbNoZero  [labelMax];
    int nbNoSignal[labelMax];
    int nbHit     [labelMax];

    // while(reader.readHit())
    while(true)
    {
      // 0. Read the data 
        timerRead.StartProfiling();
      auto breaking = !reader.readHit();
        timerRead.StopProfiling();
      
      // 0.1 Manage the end of the reading (end of file or max number of hits reached)
      if (breaking) break;
      if ((hitsMaxSet && nbHitsMax < reader.nbHits()) || (hitsMaxTotSet && nbHitsMaxTot < nbHitsTot)) break;

      ++nbHitsTot;

      // 0.2 Print the hits number
      if (reader.nbHits() > 0 && reader.nbHits() % int(1e5) == 0) 
      {
        if (group) printsln(nicer_double(reader.nbHits(), 1), "hits in .caendat", nicer_double(tree->GetEntries(), 1), "evts in .root       ");
        else       printsln(nicer_double(reader.nbHits(), 1), "hits in .caendat", nicer_double(tree->GetEntries(), 1), "hits in .root       ");
      }

      // 1. Apply the cfd
        timerCFD.StartProfiling();
      if (applyCFD && inHit.hasTrace() && useCFD[inHit.board_ID])
      {
        ++nbHit[inHit.label];
        cfd.generate(inHit.trace, CFD::sShifts[inHit.board_ID], CFD::sFractions[inHit.board_ID], 10);
        // CFD cfd(std::move(inHit.getTrace()), CFD::sShifts[inHit.board_ID], CFD::sFractions[inHit.board_ID], 10);
        auto zero = cfd.findZero();
        // auto zero = cfd.findZero(CFD::sThresholds[inHit.board_ID]); // Not using thresholds anymore
             if (zero==CFD::noSignal) {inHit.time = inHit.precise_ts; ++nbNoSignal[inHit.label];}
        else if (zero==CFD::noZero  ) {inHit.time = inHit.precise_ts; ++nbNoZero  [inHit.label];}
        else                          {inHit.time = inHit.extended_ts + zero*Caen1725::ticks_to_ps; inHit.wfa_success = true;}
      }
      else inHit.time = inHit.precise_ts;
        timerCFD.StopProfiling();

      //1.1 Apply the time shifts if registered
        timerTShift.StartProfiling();
      if (useTimeShifts && static_cast<size_t>(inHit.label) < timeshifts.size()) 
      {
        inHit.caen_time += timeshifts[inHit.label];
        inHit.time      += timeshifts[inHit.label];
      }
        timerTShift.StopProfiling();

      // 2. Fill the event builder buffer
      timerCopy.StartProfiling();
      if (!writeTraces) inHit.trace.clear();
      bool const bufferFull = !eventBuilder.fill_buffer(std::move(inHit));
      timerCopy.StopProfiling();
      if (bufferFull) dumpBufferInTree(); // If buffer is full, fill the tree and clear the buffer
    }
    
    print();
    print("nbNoSignal");
    for (size_t label_i = 0; label_i<labelMax; ++label_i) if (nbNoSignal[label_i]>0) print(label_i, nbNoSignal[label_i], 100.*nbNoSignal[label_i]/double(nbHit[label_i]), "%");
    print("nbNoZero");
    for (size_t label_i = 0; label_i<labelMax; ++label_i) if (nbNoZero  [label_i]>0) print(label_i, nbNoZero  [label_i], 100.*nbNoZero  [label_i]/double(nbHit[label_i]), "%");

    dumpBufferInTree(); // Fill the tree with the last event

    rootFile->cd();

    print();
    if (group) print(nicer_double(tree->GetEntries()), "evts in .root       ");
    else       print(nicer_double(tree->GetEntries()), "hits in .root       ");
      Profiler timerWrite;
    tree->Write();
      timerWrite.StopProfiling();
    rootFile->Close();
    
    #ifdef PROFILE
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
    #endif //PROFILE

    print(rootFile->GetName(), "written");
  }

  print(timer());
  return 0;
}