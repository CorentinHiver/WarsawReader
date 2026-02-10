#include "../CaenLib/CaenRawReader.hpp"
#include "../CaenLib/RootReader.hpp"
#include "TSystem.h"

#include <iostream>
#include <readline/readline.h>
#include <readline/history.h>


void caendatIntegrityCheck(std::string filename = "")
{
  if (filename == "")
  {
    char* input;
    input = readline("enter the name of the .caendat file to inspect: ");
    add_history(input); // TODO : add a / instead of a space when tab completion on a folder....
    filename = input;
    filename = removeBlankSpace(filename);
  }

  CaenRawReader1725::skipData();

  std::string user_input;
  File file(filename);
  if (file.extension() == "caendat")
  {
    CaenRawReader1725 reader(filename);

    int i = 0;
    
    while(true)
    {
      ++i;
      try
      {
        if (!reader.readBoardAggregate()) break;
      }
      catch (CaenRawReader1725::ErrorEof const & errorEof)
      {
        auto & board = reader.getBoard();
        print("Early end of file :", i, "th board with size", board.size, "octets and actual size", board.read_size, "octets");
        break;
      }
      catch (Caen1725::CheckBinMissed const & checkBinMissed)
      {
        auto & board = reader.getBoard();
        print("Bin missmatch :", i, "th board with size", board.size, "B and actual size", board.read_size, "B");
        break;
      }
    }
  }

  else // Not a .caendat
  {
    error("Can't read anything else but .caendat files");
  }
}

int main(int argc, char** argv)
{
  if (argc == 1) error("Please give a file name");
  else caendatIntegrityCheck(argv[1]);
}

// /home/corentin/data/tests/run_0102_06-06-2025_12h09m04s/eagleRU_i2628_0102_0000.caendat
// g++ -o check UserFriendly/caenIntegrityCheck.C -Wall -Wextra `root-config --cflags` `root-config --glibs` -g -lreadline -lhistory