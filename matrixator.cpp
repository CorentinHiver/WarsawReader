#include "CaenLib/CaenRootReader.hpp"
#include "CaenLib/RootReader.hpp"

double matrixator(std::vector<std::string> filenames, int nbHitsMax = -1)
{
  for (auto const & filename : filenames)
  {
    if (Colib::extension(filename) == "root")
    {
      RootReader reader(filename, nbHitsMax);
      while(reader.readNextEvent())
      {
        
      }
    }
  }
}