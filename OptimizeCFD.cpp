#include "AnalysisLib/CFDOptimizer.hpp"
#include "Arguments.hpp"

int main(int argc, char** argv)
{
  Arguments args(argc, argv);
  auto nb_events_max = Colib::max<size_t>();
  while(args.next())
  {
    if (args == "-n") nb_events_max = args.load<size_t>();
  }
  bool const max_events = nb_events_max < Colib::max<size_t>();

  OptimizeCFD optimizer;
  for (auto const & filename : filenames)
  {
    Caen1725::RootReader reader(filename);
    while(reader.readNextEvent())
    {
      if (max_events && nb_events_max < reader.getCursor()) break;
      auto event = reader.getEvent();
      for (auto const & hit : event)
    }
  }

  return 0;
}