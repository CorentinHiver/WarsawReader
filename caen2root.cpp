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
    {0, 7},
    {1, 7}/*,
    {6, 2},
    {7, 2},
    {8, 2}*/
  };

  CFD::sFractions = { // BOARD_ID, fraction
    {0, 0.75},
    {1, 0.75}/*,
    {6, 0.75},
    {7, 0.75},
    {8, 0.75}*/
  };

  auto useCFD = LUT<9*16>([&](int boardID)
  {
    return key_found(CFD::sShifts, boardID);
  });

  Timer timer;

  // Parameters :
  std::string outpath = "./";
  bool readTraces = true;
  std::vector<bool> boardReadTrace(100, true);
  bool writeTraces = false;
  size_t nbHitsMax = -1;
  size_t nbHitsMaxTot = -1;
  bool hitsMaxSet = false;
  bool hitsMaxTotSet = false;
  bool ts_evt_build = false;
  bool group = true;
  
  std::vector<std::string> filenames;
  Timeshifts timeshifts;

  uint64_t time_window = 2e6 ; // ps

  std::vector<UShort_t> trigger_labels;

  auto printHelp = [](){ 
    print("caen2root usage");
    print("-e --ts-evt-build      [0 or 1] (default 0). Perform event building based on the raw timestamp instead of the absolute time (usually corrected by cfd).");
    print("-f --files             [caen_filename] : File to convert. Include wildcards * and ?, but ONLY IF the name is guarded by quotes (i.e. -f \"/path/to/file/names*.caendat\") ");
    print("-F --files-nb          [caen_filename] [nb_files] : Files to convert. Include wildcards * and ?, but ONLY IF the name is guarded by quotes (i.e. -f \"/path/to/file/names*.caendat\"). For nb_files, -1 = all, scientific format accepted (e.g. 1.e3 (=1000))]");
    print("-g --group             [0 or 1] (default 1) : Sets the output format. 0 : plain tree with additionnal event number and multiplicity fields. 1 : each leaf is a vector.");
    print("-h --help              : print this help");
    print("-n --hits-nb           [nb_hits (-1 = all, 1.e3 (=1000) format accepted)] : maximum number of hits to be read by the programm IN EACH FILE");
    print("-N --hits-tot-nb       [nb_hits (-1 = all, 1.e3 (=1000) format accepted)] : maximum number of hits to be read by the programm (if multithread : not thread safe)");
    print("-o --output            [output_path]");
    print("   --read-traces       [0 or 1] (default 1) : Read the traces for all boards.");
    print("   --board-skip-trace  [boardID] : Skip the trace for this boad ID. Used only if --read-traces is true.");
    print("   --write-traces      [0 or 1] (default 0) : Write the trace in the output root file.");
    print("-T --timeshifts        [filename] : List of timestamp shifts. Format : in each line : global_label timeshift[ns].");
    print("-t --trigger [option] : trigger on the given global label or board ID, or a user-defined file with a list of labels. Example : -t -l 0 to trigger on label 0");
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
    print("\t read-traces 1 \\");
    print("\t write-traces 1 \\");
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
      else if (temp == "-n")
      {
        double tmp_d = 0; iss >> tmp_d;
        nbHitsMax = static_cast<size_t>(tmp_d);
        hitsMaxSet = true;
        print("Maximum hits per file =", nbHitsMax);
      }
      else if (temp == "-N")
      {
        double tmp_d = 0; iss >> tmp_d;
        nbHitsMaxTot = static_cast<size_t>(tmp_d);
        hitsMaxTotSet = true;
        print("Maximum hits total", nbHitsMaxTot);
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
      }
      else if (temp == "-tw" || temp == "time-window")
      {
        double e; iss >> e;
        time_window = e*1000;
      }
      else
      {
        throw_error("Unkown " + temp + " command");
      }
    }

    // -------------------- //
    // Check the parameters //
    // -------------------- //

    // if (hitsMaxTotSet && multithreadSet) throw-error("Can't have -N and -M at the same time !");
    print(trigger_labels);

    if (filenames.empty()) throw_error("No files !!");

    size_t nbHitsTot = 0;

    for (auto const & filename : filenames)
    {
      File file(filename);

      if (!file) {error("can't find ", file); continue;}
      
      Caen1725RootInterface reader(filename, readTraces);
      reader.setBoardReadTrace(boardReadTrace);
      Caen1725EventBuilder eventBuilder(reserved_buffer_size);
      eventBuilder.buildOnTimestamp(ts_evt_build);

      auto rootFilename = outpath + file.shortName()+".root";
      auto rootFile = TFile::Open(rootFilename.c_str(), "recreate");
      gROOT->cd();
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
      Timer timerBufferFill;
      Timer timerTrigger;
      Timer timerEvtBuild;
      Timer timerFill;
      Timer timerClear;

      // To investigate : might not be optimized
      // Pre-declaration of this piece of code
      auto fillTreeAndClear = [&]() -> void
      {
        // 3. Perform the event building 
          timerEvtBuild.Start();
        eventBuilder.fast_event_building(time_window);
          timerEvtBuild.Stop();

        // 4. Write the events to the ROOT tree
        for (auto const & event : eventBuilder)// Loop over all the event in buffer :
        {
          evtMult = event.size();

        #ifdef TRIGGER
          // TODO:
        #endif //TRIGGER

          // 4.1 Apply the trigger
            timerTrigger.Start();
          static thread_local bool trigger; trigger = true;
          static thread_local bool trigger_label = !trigger_labels.empty();
          if (trigger_label)
          {
            trigger = false;
            for (auto const & hit_i : event) if (Colib::found(trigger_labels, eventBuilder[hit_i].label)) trigger = true;
          }
            timerTrigger.Stop();

          if (trigger) 
          {
            if (group)
            {
              // 4.2 Write the event
              for (auto const & hit_i : event)
              {
                  timerCopy.Start();
                outEvent.push_back(eventBuilder[hit_i]);
                  timerCopy.Stop();
              }

                timerFill.Start();
              tree -> Fill();
                timerFill.Stop();

                timerClear.Start();
              outEvent.clear();
                timerClear.Stop();
              ++outEvent.evtNb;
            }
            else for (auto const & hit_i : event)
            {
              eventBuilder[hit_i].rel_time = eventBuilder[hit_i].time - eventBuilder[0].time;
                timerCopy.Start();
              outHit = std::move(eventBuilder[hit_i]);
                timerCopy.Stop();
                
                timerFill.Start();
              tree -> Fill();
                timerFill.Stop();
              ++evtNb;
            }
          }
        }

        // 5. Reset event builder

        eventBuilder.clear();
      };

      int nbNoSignal = 0;
      int nbNoZero = 0;

      // while(reader.readHit())
      while(true)
      {
        // 0. Read the data 
          timerRead.Start();
        auto breaking = !reader.readHit();
          timerRead.Stop();
        
        // 0.1 Manage the end of the reading (end of file or max number of hits reached)
        if (breaking) break;
        ++nbHitsTot;
        if (hitsMaxSet && nbHitsMax < reader.nbHits()) break;
        if (hitsMaxTotSet && nbHitsMaxTot < nbHitsTot) break;
        if (reader.nbHits() > 0 && reader.nbHits() % int(1e5) == 0) printsln(nicer_double(reader.nbHits(), 1), "    ");

        // 1. Apply the cfd
          timerCFD.Start();
        if (applyCFD && readTraces && inHit.hasTrace() && useCFD[inHit.board_ID])
        {
          CFD cfd(*inHit.getTrace(), CFD::sShifts[inHit.board_ID], CFD::sFractions[inHit.board_ID]);
          auto zero = cfd.findZero();
          // auto zero = cfd.findZero(CFD::sThresholds[inHit.board_ID]); // Not using thresholds anymore
               if (zero==CFD::noSignal) {inHit.time = inHit.precise_ts; ++nbNoSignal;}
          else if (zero==CFD::noZero  ) {inHit.time = inHit.precise_ts; ++nbNoZero  ;}
          else                           inHit.time = inHit.extended_ts + zero * CaenDataReader1725::ticks_to_ps;
        }
        else inHit.time = inHit.precise_ts;
          timerCFD.Stop();

        //1.1 Apply the time shifts if registered
          timerTShift.Start();
        if (static_cast<size_t>(inHit.label) < timeshifts.size()) 
        {
          inHit.timestamp += timeshifts[inHit.label];
          inHit.time      += timeshifts[inHit.label];
        }
          timerTShift.Stop();

        // 2. Fill the event builder buffer
        
        if (!eventBuilder.fill_buffer(std::move(inHit))) fillTreeAndClear(); // If buffer is full, fill the tree and clear the buffer
      }
      
      print();
      print("noSignal", nbNoSignal);
      print("noZero"  , nbNoZero  );

      fillTreeAndClear(); // Fill the tree with the last event

      rootFile->cd();

      print();
      if (group) print(tree->GetEntries(), "events in the tree");
      else       print(tree->GetEntries(), "hits in the tree"  );
        Timer timerWrite;
      tree->Write();
        timerWrite.Stop();
      rootFile->Close();
      
      if (true)
      {
        print("timeRead", timerRead.timeElapsed());
        print("timeCFD", timerCFD.timeElapsed());
        print("timeTShift", timerTShift.timeElapsed());
        print("timeEvtBuild", timerEvtBuild.timeElapsed());
        print("timeCopy", timerCopy.timeElapsed());
        print("timerTrigger", timerTrigger.timeElapsed());
        print("timeBufferFill", timerBufferFill.timeElapsed());
        print("timeFill", timerFill.timeElapsed());
        print("timeClear", timerClear.timeElapsed());
        print("timeWrite", timerWrite.timeElapsed());
      }

      print(rootFile->GetName(), "written");
    }
    
  }
  print(timer());
  return 0;
}