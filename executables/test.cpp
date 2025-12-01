#include "../CaenLib/Caen1725.hpp"

int main()
{
  Digitizer board;
  board.connect();
  board.start();
  return 0;
}