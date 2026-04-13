#include "AnalysisLib/CFDOptimizer.hpp"
#include "CaenLib/RootReader.hpp"
// #include "Arguments.hpp"
#include "Colib/lib/Classes/Arguments.hpp"

using namespace std;
using namespace Colib;
using namespace Caen1725;

constexpr Label refLabel = 81;
int main(int argc, char** argv)
{
  Arguments args(argc, argv);
  auto nb_events_max = max<size_t>();
  vector<string> filenames;
  string parameterFile = "cfd.opt.param";
  while(args.next())
  {
    if (args == "-n") nb_events_max = size_cast(args.load<double>());
    else if (args == "-f") filenames = findFilesWildcard(args.load<string>());
    else if (args == "-p") parameterFile = args.load<string>();
    else throw_error("Unkown argument "+args.getArg());
  }
  bool const max_events = nb_events_max < max<size_t>();
  print(nb_events_max);
  if (filenames.empty()) throw_error("No files !! Use -f options to feed me.");

  CFDOptimizer optimizer({0, 2, 4});
  CFDMinimisationParameters parameters;
  parameters.set(parameterFile);
  print(parameters);
  optimizer.setParameters(parameters);

  for (auto const & filename : filenames)
  {
    RootReader reader(filename);
    printsln(filename);
    while(reader.readNextEvent())
    {
      if (max_events && nb_events_max < reader.getCursor()) break;
      auto event = reader.getEvent();
      Hit refHit;
      bool hasRef = false;
      for (int hit_i = 0; hit_i<event.mult; ++hit_i) if (event.label[hit_i] == refLabel) 
      {
        hasRef = true;
        refHit = event.getHit(hit_i);
        break;
      }
      if (!hasRef) continue;
      for (int hit_i = 0; hit_i<event.mult; ++hit_i)
      {
        if (event.label[hit_i] == refLabel) continue;
        auto const & trace = event.traces[hit_i];
        optimizer.calculate_dT(refHit, event.label[hit_i], event.time[hit_i], trace);
      }
    }
  }

  printsln("dT calculated, finding optimal");

  optimizer.calculateResolutions();
  optimizer.write("cfdOpti.root");

  return 0;
}