#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <unordered_map>
#include <sstream>
#include <thread>

struct CFDParameters
{
  double fraction = 0.5;
  int shift = 1;
  int nbBaseline = 10;
};

template<class T>
class CFDParametersMap
{
public:
  std::unordered_map<T, CFDParameters> map;
  CFDParametersMap (std::unordered_map<T, CFDParameters> const & _map):
    map(_map) {}

  void load(std::string const & filename)
  {
    std::ifstream paramFile(filename);
    std::string line;
    while(std::getline(paramFile, line))
    {
      // Get the parameter of each board or each detector
      std::istringstream iss(line);
      T label; iss >> label;
      auto & params = map[label];
      iss >> params.fraction >> params.shift;
    }
  }
};

class CFD
{
protected:
  /////////////////////////////////////////
  // Some helper functions (from Colib): //
  /////////////////////////////////////////

  template <typename T>
  std::tuple<T, size_t> minimum_and_index(std::vector<T> const & vector)
  {
    auto const & min_it = std::min_element(std::begin(vector), std::end(vector));
    auto const & min_index = static_cast<size_t>(std::distance(std::begin(vector), std::min_element(std::begin(vector), std::end(vector))));
    return std::make_tuple(*min_it, min_index);
  }

  // Type name
  using size_t = std::size_t;
  using CFD_t = std::vector<double>;
  
  // Type check
  template <typename T, typename std::enable_if<std::is_floating_point<T>::value, bool>::type = true>
  inline static constexpr bool is_floating() noexcept { return true;}
  template <typename T, typename std::enable_if<!std::is_floating_point<T>::value, bool>::type = true>
  inline static constexpr bool is_floating() noexcept { return false;}

  // Console colors
  struct Color
  {
    static constexpr const char* RED   = "\u001b[31m";
    static constexpr const char* RESET = "\u001b[0m" ;
  };
  
public:
  /// @brief Default constructor
  CFD() noexcept = default; 
  
  /// @brief Constructs a CFD object from a trace
  template<class T>
  CFD(std::vector<T> const & _trace, int shift, double fraction, size_t nbSamplesBaseline = 10)
  {
    generate(_trace, shift, fraction, nbSamplesBaseline);
  }

  /** @brief Constructs a CFD object from a trace
  * @details 
  * Fills an internal vector with a trace after the m_baseline is determined 
  * as the mean value of the first nbSamplesBaseline samples, 
  * and subtracted to all the samples of the trace. Call generate(shift, fraction)
  * to actually generate the CFD of the trace 
  **/
  template<class T = double>
  CFD(std::vector<T> const & _trace, size_t nbSamplesBaseline) 
  {
    setTrace(_trace, nbSamplesBaseline);
  }

  template<class T = double>
  CFD& setTrace(std::vector<T> const & _trace, size_t nbSamplesBaseline)
  {    
    if (_trace.empty()) return *this;

    m_baseline = 0;
    for (size_t sample_i = 0; sample_i<nbSamplesBaseline; ++sample_i) m_baseline += _trace[sample_i];
    m_baseline /= nbSamplesBaseline;
    
    trace.reserve(_trace.size());

    for (auto const & sample : _trace)
    {
      if constexpr (is_floating<T>()) trace.push_back(sample - m_baseline);
      else trace.push_back(adc_to_double(sample) - m_baseline);
    }
    return *this;
  }

  template<class T>
  void setBaseline(std::vector<T> const & _trace, int nbSamplesBaseline)
  {
    if (_trace.empty()) return;
    auto const nbl = std::min(static_cast<size_t>(nbSamplesBaseline), _trace.size()); // Number points in for baseline
    m_baseline = std::accumulate(_trace.begin(), _trace.begin()+nbl, 0) / nbl;
  }

  template<class T>
  void generate(std::vector<T> const & _trace, int shift, double fraction)
  {
    if (_trace.empty() || shift < 0 || _trace.size() <= static_cast<size_t>(shift) ) return;

    cfd.clear();
    cfd.reserve(_trace.size()-shift);
    for (size_t bin = shift; bin < _trace.size(); ++bin)
    {
      double delayed = static_cast<double>(_trace[bin - shift]) - m_baseline;
      double sample  = static_cast<double>(_trace[bin]        ) - m_baseline;
      cfd.push_back(delayed - (fraction * sample));
    }
  }

  template<class T>
  void generate(std::vector<T> const & _trace, int shift, double fraction, size_t nbSamplesBaseline)
  {
    setBaseline(_trace, nbSamplesBaseline);
    return generate(_trace, shift, fraction);
  }

  

  // template<class T>
  // void generate(int shift, double fraction)
  // {
  //   cfd.clear();
  //   cfd.reserve(trace.size());
  //   for (size_t bin = 5*shift; bin < (trace.size() - shift); ++bin)
  //   {
  //     auto const & value = trace[bin - shift] - fraction * trace[bin];
  //     cfd.push_back(value);
  //   }
  // }

  /// @brief Calculates the last zero crossing before the calculated cfd signal goes above the given threshold
  double findZero(double threshold)
  {
    if (threshold < 0) std::cout << Color::RED << "in CFD::findZero(threshold) : threshold < 0 !" << Color::RESET << std::endl;
    for (size_t bin_i = 0; bin_i < cfd.size(); ++bin_i){  // Loop through the cfd values
      if (threshold < cfd[bin_i]){                        // The cfd value crossed the threshold
        for (size_t bin_j = bin_i; bin_j>0; --bin_j){     // Looping back for looking for the zero crossing
          if (cfd[bin_j] < 0) return interpolate0(bin_j); // Zero crossing found, return the interpolated zero crossing between samples before and after
        }
        return noZero; // The 0 crossing happened before the first sample, so impossible to determine it
      }
    }
    return noSignal; // The signal never crosses the threshold -> the signal is too small, the cfd parameters are wrong, or the threshold is too large
  }

  /// @brief Finds the last zero crossing before the cfd trace reaches its maximum
  double findZero()
  {
    auto maximum_it = std::max_element(std::begin(cfd), std::end(cfd));
    if (*maximum_it < 0) return noSignal;            // If never crosses zero, returns noSignal
    auto const & max_bin = std::distance(std::begin(cfd), maximum_it); // Get the maximum bin number

    for (size_t bin_j = max_bin; bin_j>0; --bin_j)    // Looping back to look for the zero crossing
      if (cfd[bin_j] < 0) return interpolate0(bin_j); // Zero crossing found, return the interpolated zero crossing between samples before and after
    return noZero; // The 0 crossing happened before the first sample, so impossible to determine it
  }
  
  CFD_t trace;
  CFD_t cfd;

  auto size() const {return cfd.size();}

  // Static variables :
  constexpr static double noZero   = 1e-100;
  constexpr static double noSignal = 1e-101;

  /////////////////////////
  // Parameters handling //
  /////////////////////////
  
protected:

  double m_baseline = 0;
  
  inline double interpolate0(size_t bin) const noexcept 
  {
    auto const & cfd_0 = cfd[bin  ];
    auto const & cfd_1 = cfd[bin+1];
    if( cfd_0 == cfd_1) return 0;
    else return (bin -        cfd_0 /
                         (cfd_1 - cfd_0));
  }
  
public:
};