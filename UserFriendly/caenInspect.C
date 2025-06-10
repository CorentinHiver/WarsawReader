#include "../CaenLib/CaenRootReader.hpp"
#include "../CaenLib/RootReader.hpp"
#include "TSystem.h"

#include <iostream>
#include <readline/readline.h>
#include <readline/history.h>


void caenInspect(std::string filename = "")
{
  if (filename == "")
  {
    char* input;
    input = readline("enter the name of the .caendat file to inspect: ");
    add_history(input); // TODO : add a / instead of a space when tab completion on a folder....
    filename = input;
    filename = removeBlankSpace(filename);
  }

  std::string user_input;
  File file(filename);
       if (file.extension() == "caendat")
  {
    CaenRootReader reader(filename);
    while(reader.readHit())
    {
      auto const & nbLines = Colib::getTerminalRows();
      if (reader.nbHits() % nbLines == 0) {
        std::getline(std::cin, user_input);
        if (user_input == "q") break;
      }
      print(reader.getHit());
    }
  }
  else if (file.extension() == "root")
  {
    RootReader reader(filename);
    int evtNb = 0;
    reader.getTree()->SetBranchAddress("evtNb", &evtNb);
    while(reader.readNext())
    {
      auto const & nbLines = Colib::getTerminalRows();
      if (reader.getCursor() % nbLines == 0) {
        std::getline(std::cin, user_input);
        if (user_input == "q") break;
      }
      print(evtNb, reader.getHit());
    }
  }
}