#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include "utils.hpp"

namespace CaenDataReader
{
  struct Sample
  {
    uint16_t sample = 0;
    bool     DP1    = false;
    bool     T      = false;
  
    friend std::ostream& operator<<(std::ostream& out, Sample const & sample) {
      out << 
        " sample " << sample.sample           <<
        " DP1 "    << nicer_bool(sample.DP1)  <<
        " T "      << nicer_bool(sample.T);
      return out;
    }
  };
};


#endif //SAMPLE_HPP