// g++ -o rootReaderExample rootReaderExample.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -g

#include "CaenLib/RootReader.hpp" // You can simply source CaenLib folder for a global access , no installation required
#include "LibCo/libCo.hpp"

// using namespace Colib // This is recommended if you want to keep a code easier to read, although use it only in .cpp or .C files

void rootReaderExample(std::vector<std::string> filenames, std::string program_str = "print", int nbHitsMax = -1)
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
      print(reader.getEvent()); // Printing the whole event
      // ---------------------------------------- //
      // C. Looping through the hits in the event //
      // ---------------------------------------- //
    
      // Two ways to loop through the hits of the event :
      // 1. range-based loop (prettier)
      // for (auto const & hit : reader.getEvent())
      // 2. Classic loop, get access to the hit index in the event hit_i :
      for (size_t hit_i = 0; hit_i<reader.getEvent().size(); ++hit_i) 
      { // Printing hit by hit
        auto const & hit = reader.getEvent()[hit_i]; // Alias to the hit to read (not efficient because creating a hit)
        print(hit);                                  // Prints the hit to console
        // auto const & Ecal   = hit.adc  * 1 + 0;      // Calibration example
        // auto const & time_s = hit.time * 1e12 ;      // Converting the absolute time in ps to seconds
      }
      Colib::pause();
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
    print("-p [name]: What program to run. Possibilities : print dt. Default : print");
  }
  std::string command;
  std::string program;
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
  if (filenames.empty()) Colib::throw_error("No file !!");
  rootReaderExample(filenames, program, nb_hits);
}
