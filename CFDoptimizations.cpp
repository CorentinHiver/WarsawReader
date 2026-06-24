#include "Colib/lib/CoMT.hpp"
#include "Colib/lib/Classes/Arguments.hpp"

#include "AnalysisLib/CFDOptimizer.hpp"

#include "CaenLib/RootReader.hpp"

using namespace std;
using namespace Colib;
using namespace Caen1725;

int main(int argc, char** argv)
{
  int refLabelI{-1};
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
    else if (args == "-t") refLabelI = args.load<Label>();
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
  if (refLabelI<0) throw_error("No reference label !! Use option -t to give it to me.");
  if (filenames.empty()) throw_error("No files !! Use -f options to feed me.");
  auto refLabel = size_cast(refLabelI);

  CFDMinimisationParameters parameters(parameterFile);
  print(parameters);

  // CFDOptimizer optimizer;
  // optimizer.setParameters(parameters);
  CFDOptimizer optimizer(parameters);

  auto distributed_filenames = MT::distribute(filenames);
  MT::parallelise_function([&](){
    for (auto const & filename : filenames)
    {
      RootReader reader(filename);
      reader.printEvery(1000);
      reader.setMaxHits(nb_events_max);
      printsln(filename);
      while(reader.readNextEvent())
      {
        auto const & event = reader.getEvent();
        for (int hit_i = 0; hit_i<event.mult; ++hit_i) if (event.label[hit_i] == refLabel) 
        {
          for (int hit_j = 0; hit_j<event.mult; ++hit_j)
          {
            auto const & label = event.label[hit_j];
            auto const & trace = event.traces[hit_j];
            if (label != refLabel)
            {
              optimizer.calculate_dT(event.time[hit_i], event.time[hit_j], label, trace, ticks_to_ps);
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