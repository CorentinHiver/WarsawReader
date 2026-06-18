#pragma once
#include "Hit.hpp"

namespace Caen1725
{
  struct Coeffs {double intercept, slope; double quadratic = 0; int m_order = 1;};
  class Calibration
  {
  public:
    Calibration() noexcept = default;
    Calibration(std::string const & filename) {load(filename);}
  
    void load(std::string const & filename)
    {
      std::ifstream file(filename);
      if (!file.is_open()) Colib::throw_error("In calibration::load(filename) ", filename, " not found !!");
      std::string line;
      while(getline(file, line))
      {
        std::istringstream iss;
        Label label = -1; iss >> label;
        auto & coeff = m_coeffs[label];
        iss >> coeff.intercept >> coeff.slope >> coeff.intercept >> coeff.quadratic;
      }
    }
  private:
    std::unordered_map<Label, Coeffs> m_coeffs;
  };
}