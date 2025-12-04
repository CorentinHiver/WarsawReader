// g++ -o rootReaderExample rootReaderExample.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -g

#include "CaenLib/RootReader.hpp" // You can simply source CaenLib folder for a global access , no installation required

void rootReaderExample(std::vector<std::string> filenames, int nbHitsMax = -1)
{
  for (auto const & filename : filenames)
  {
    if (Colib::extension(filename) != "root")
    {
      print("Skipping file", filename, "that is without .root extension");
      continue;
    }

    // --------------------------------------------- //
    // A. Attaching the data reader to the root file //
    // --------------------------------------------- //
    RootReader reader(filename, nbHitsMax); 

    // ---------------------------------------------- //
    // B. Looping through the events in the root file //
    // ---------------------------------------------- //
    while(reader.readNextEvent()) 
    {
      // ---------------------------------------- //
      // C. Looping through the hits in the event //
      // ---------------------------------------- //
    
      // Two ways to loop through the hits of the event :
      // 1. range-based loop (prettier)
      // for (auto const & hit : reader.getEvent())
      // 2. Classic loop, get access to the hit index in the event hit_i :
      for (size_t hit_i = 0; hit_i<reader.getEvent().size(); ++hit_i)
      {
        auto const & hit = reader.getEvent()[hit_i]; // Alias to the hit to read
        print(hit);                                  // Prints the hit to console
        auto const & Ecal   = hit.adc  * 1 + 0;       // Calibration example
        auto const & time_s = hit.time * 1e12 ;       // Converting the hit absolute time in ps
      }
    }
  }
}

int main(int argc, char** argv)
{
  std::istringstream iss(Colib::argv_to_string(argv));
  if (argc == 1)
  {
    print("rootReaderExample usage : ./rootReaderExample [[parameters]]");
    print("-f [filename.root or filenames*.root]");
    print("-n [number of hits]");
  }
  std::string command;
  std::vector<std::string> filenames;
  int nb_hits = -1;
  while(iss >> command)
  {
    if (command == "-f") 
    {
      iss >> command;
      for (auto const & file : Colib::findFilesWildcard(command)) filenames.push_back(file);
    }
    else if (command == "-n")
    {
      double e; iss >> e;
      nb_hits = static_cast<int>(e);
    }
  }
  rootReaderExample(filenames, nb_hits);
}
