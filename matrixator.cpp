#include "CaenLib/CaenRootReader.hpp"
#include "CaenLib/RootReader.hpp"

double matrixator(std::vector<std::string> filenames, int nbHitsMax = -1)
{
  for (auto const & filename : filenames)
  {
    if (extension(filename) == "root")
    {
      RootReader reader(filename);
      while(reader.fillEvent())
      {
        if (reader.eventReady())
        {
          
        }
      }
    }
  }
}