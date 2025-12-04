// #include "/home/corentin/Installation/CaenLibs/CAENDigitizer-v2.19.0/include/CAENDigitizer.h"
#include "../CaenLib/Caen1725.hpp"

int main()
{
  
  std::cout << Digitizer::listConnected().size() << " digitizers" << std::endl;
  for (auto const & e : Digitizer::listConnected()) std::cout << e << std::endl;

  // Digitizer board;
  // board.connect(0, 0);
  // board.start();
  return 0;
}