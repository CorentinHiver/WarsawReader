#include "../CaenLib/CaenRawReader.hpp"
#include "../CaenLib/RootReader.hpp"
#include "TSystem.h"

#include <iostream>
#include <readline/readline.h>
#include <readline/history.h>


void caenIntegrityCheck(std::string filename = "")
{
  if (filename == "")
  {
    char* input;
    input = readline("enter the name of the .caendat file to inspect: ");
    add_history(input); // TODO : add a / instead of a space when tab completion on a folder....
    filename = input;
    filename = removeBlankSpace(filename);
  }

  CaenRawReader::skipData();

  std::string user_input;
  File file(filename);
       if (file.extension() == "caendat")
  {
    CaenRawReader reader(filename);

    int i = 0;
    
    auto checkBoard = [&]() -> bool 
    {
      ++i;
      try
      {
        return reader.readBoardAggregate();
      }
      catch (CaenRawReader::ErrorEof const & errorEof)
      {
        auto & board = reader.getBoard();
        print("Early end of file : board", i, "with size", board.size, "B and actual size", board.read_size, "B");
        return false;
      }
    };

    while(checkBoard()) continue;
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

int main(int argc, char** argv)
{
  if (argc == 1) error("Please give a file name");
  else caenIntegrityCheck(argv[1]);
}

// /home/corentin/data/tests/run_0102_06-06-2025_12h09m04s/eagleRU_i2628_0102_0000.caendat
// g++ -o check UserFriendly/caenIntegrityCheck.C -Wall -Wextra `root-config --cflags` `root-config --glibs` -g -lreadline -lhistory