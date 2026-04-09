#include "AnalysisLib/CFDOptimizer.hpp"
#include "Arguments.hpp"
#include "CaenLib/RootReader.hpp"

using namespace std;
using namespace Colib;
using namespace Caen1725;
int main(int argc, char** argv)
{
  Arguments args(argc, argv);
  auto nb_events_max = max<size_t>();
  vector<string> filenames;
  while(args.next())
  {
    if (args == "-n") nb_events_max = args.load<size_t>();
    else if (args == "-f") filenames = findFilesWildcard(args.load<string>());
    else throw_error("Unkown argument "+args.getArg());
  }
  bool const max_events = nb_events_max < max<size_t>();
  if (filenames.empty()) throw_error("No files !! Use -f options to feed me.");

  CFDOptimizer optimizer;
  for (auto const & filename : filenames)
  {
    RootReader reader(filename);
    while(reader.readNextEvent())
    {
      if (max_events && nb_events_max < reader.getCursor()) break;
      auto event = reader.getEvent();
      for (int hit_i = 0; hit_i<event.mult; ++hit_i)
      {
        auto const & trace = event.traces[hit_i];
        
      }
      print();
    }
  }

  return 0;
}