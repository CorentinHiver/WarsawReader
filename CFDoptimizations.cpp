#include "Colib/lib/CoMT.hpp"
#include "Colib/lib/Classes/Arguments.hpp"
#include "AnalysisLib/CFDOptimizer.hpp"
#include "CaenLib/RootReader.hpp"

using namespace std;
using namespace Colib;
using namespace Caen1725;

constexpr Label refLabel = 81;
constexpr Label refBoard = 5;
int main(int argc, char** argv)
{
  Arguments args(argc, argv);
  auto nb_events_max = max<size_t>();
  vector<string> filenames;
  string parameterFile = "cfd.opt.param";
  bool onlyMin{}; 
  std::string rootfilename;
  // Calibration calib;
  while(args.next())
  {
    if (args == "-n") nb_events_max = size_cast(args.load<double>());
    else if (args == "-f") filenames = findFilesWildcard(args.load<string>());
    // else if (args == "-c") calib.load(args.load<string>());
    else if (args == "-p") parameterFile = args.load<string>();
    else if (args == "min") 
    {
      onlyMin = true;
      rootfilename = args.load<std::string>();
    }
    else throw_error("Unkown argument "+args.getArg());
  }
  if (onlyMin)
  {
    CFDOptimizer optimizer;
    optimizer.findMinima(rootfilename);
    optimizer.write_dT("cfd.params");
    return 0;
  }
  bool const max_events = nb_events_max < max<size_t>();
  if (filenames.empty()) throw_error("No files !! Use -f options to feed me.");

  std::vector<Label> labels_to_study;
  for (Label board = 0; board<2; ++board) for (Label channel = 0; channel<8; ++channel) labels_to_study.push_back(board*16+channel*2);
  auto const isStudied = LUT<200>([&labels_to_study](Label const & label){return found(labels_to_study, label);});
  CFDOptimizer optimizer(labels_to_study);
  CFDMinimisationParameters parameters;
  parameters.set(parameterFile);
  print(parameters);
  optimizer.setParameters(parameters);

  auto distributed_filenames = MT::distribute(filenames);
  MT::parallelise_function([&](){
    for (auto const & filename : filenames)
    {
      RootReader reader(filename);
      printsln(filename);
      while(reader.readNextEvent())
      {
        if (reader.getCursor()%1000 == 0) printsln(getShortname(filename), 
          nicer_double((100.*reader.getCursor())/reader.getTree()->GetEntries(), 1), " %");
        if (max_events && nb_events_max < reader.getCursor()) break;
        auto event = reader.getEvent();
        for (int hit_i = 0; hit_i<event.mult; ++hit_i) if (event.board_ID[hit_i] == refBoard) 
        {
          auto const & refHit = event.getHit(hit_i);
          for (int hit_j = 0; hit_j<event.mult; ++hit_j)
          {
            auto const & label = event.label[hit_j];
            auto const & trace = event.traces[hit_j];
            if (isStudied[label] && event.board_ID[hit_j] != refBoard)
            {
              if (trace.empty()) throw_error("No trace for detector ", label);
              optimizer.calculate_dT(refHit, label, event.time[hit_j], trace);
            }
          }
        }
      }
    }
  });

  printsln("dT calculated, finding optimal cfd parameters");

  optimizer.calculateResolutions();
  optimizer.writeRoot("cfdOpti.root");
  optimizer.findMinima("cfdOpti.root");
  optimizer.write_dT("cfd.params");

  return 0;
}