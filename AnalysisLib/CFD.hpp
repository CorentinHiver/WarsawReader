#ifndef CFD_HPP
#define CFD_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <unordered_map>
#include <sstream>

class CFD
{
protected:
  // Some helper functions (from libCo):
  using size_t = std::size_t;

  template <typename T, typename std::enable_if<std::is_floating_point<T>::value, bool>::type = true>
  inline static constexpr bool is_floating() noexcept { return true;}
  template <typename T, typename std::enable_if<!std::is_floating_point<T>::value, bool>::type = true>
  inline static constexpr bool is_floating() noexcept { return false;}

  virtual inline double fast_uniform() noexcept {
    static thread_local std::minstd_rand generator(std::random_device{}());
    static thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
  }
  static constexpr const char* RED   = "\u001b[31m";
  static constexpr const char* RESET = "\u001b[0m" ;
  
  using CFD_t = std::vector<double>;

public:

  CFD() noexcept = default;

  template<class T = double>
  CFD(std::vector<T> const & _trace, size_t nb_samples_baseline = 1) :
    m_size(_trace.size())
  {
    if (_trace.empty()) return;

    double baseline = _trace[0];
    if (nb_samples_baseline > 1)
    {
      for (size_t sample_i = 1; sample_i<nb_samples_baseline; ++sample_i) baseline += _trace[sample_i];
      baseline /= nb_samples_baseline;
    }
    
    trace.reserve(m_size);

    for (auto const & sample : _trace) 
    {
      if constexpr (is_floating<T>()) {
        trace.push_back(sample - baseline);
      }
      else {
        trace.push_back(static_cast<double>(sample  + fast_uniform()) - baseline);
      }
    }
  }

  CFD& operator=(CFD const & other)
  {
    this -> trace = other.trace;
    this -> cfd = other.cfd;
    this -> m_size = other.m_size;
    return *this;
  }

  template<class T = double>
  CFD& setTrace(std::vector<T> const & _trace)
  {
    *this = _trace;
    return *this;
  }

  template<class T = double>
  CFD& operator=(std::vector<T> const & _trace)
  {
    if (_trace.empty()) return *this;
    m_size = _trace.size();

    double baseline = 0;
    // for (size_t sample_i = 0; sample_i<nb_samples; ++sample_i) baseline += _trace[sample_i];
    // baseline /= nb_samples;
    
    trace.reserve(m_size);

    for (auto const & sample : _trace) 
    {
      if constexpr (is_floating<T>()) {
        trace.push_back(sample - baseline);
      }
      else {
        trace.push_back(static_cast<double>(sample  + fast_uniform()) - baseline);
      }
    }
    return *this;
  }

  template<class T = double>
  CFD(std::vector<T> const & _trace, int const & shift, double const & fraction) : CFD(_trace)
  {
    this -> calculate(shift, fraction);
  }

  virtual void calculate(size_t const & shift, double const & fraction)
  {
    if (fraction>1.) {std::cout << RED << "in CFD(trace, shift, fraction): fraction>1 !!" << RESET << std::endl; return;}

    cfd.clear();

    if (m_size < 2*shift) {std::cout << RED << "in CFD(trace, shift, fraction): m_size = " << m_size << " < 2*shift = " << 2*shift << " !!" << RESET << std::endl; return;}
    
    cfd.reserve(m_size - 2*shift);
    // for (size_t bin = 0; bin<shift; ++bin) cfd.push_back(0);
    for (size_t bin = 2*shift; bin<m_size - shift; ++bin){
      auto const & value = fraction * trace[bin] - trace[bin - shift] ;
      cfd.push_back(value);
    }
  }

  double findZero(double const & threshold)
  {
    if (threshold>0) std::cout << RED << "in CFD::findZero(threshold) : threshold > 0 !" << RESET << std::endl;
    for (size_t bin_i = 0; bin_i < cfd.size(); ++bin_i){  // Loop through the cfd values
      if (cfd[bin_i] < threshold){                        // The cfd value crossed the threshold
        for (size_t bin_j = bin_i; bin_j>0; --bin_j){     // Looping back for looking for the zero crossing
          if (cfd[bin_j] > 0) return interpolate0(bin_j); // Zero crossing found, return the interpolated zero crossing between samples before and after
        }
        return noZero; // The 0 crossing happened before the first sample, so impossible to determine it
      }
    }
    return noSignal; // The signal never crosses the threshold -> the signal is too small, the cfd parameters are wrong, or the threshold is too large
  }
  
  CFD_t trace;
  CFD_t cfd;

  // Static variables :
  constexpr static double noZero = 0.;
  constexpr static double noSignal = -1.;

  // TODO Static methods to fill the parameters vectors

  static void loadParameters(std::string filename)
  {
    std::ifstream paramFile(filename);
    std::string line;
    std::getline(paramFile, line);

    // Get the header. Anything can be written in it, but at least BOARD or LABELS
    // in order to know if the parameters are board wide and the label the BOARD_ID,
    // or detector per detector with the label the global label (BOARD_ID*16 + Channel_ID*2 + subchannel_ID)
         if (line.find("BOARDS") != std::string::npos) Param::sType = Param::BOARD;
    else if (line.find("LABELS") != std::string::npos) Param::sType = Param::LABEL;
    else std::cout << CFD::RED << "CFD::loadParameters " << filename << " : format issue. Should begin with LABELS or BOARDS" << RESET << std::endl; 

    while(std::getline(paramFile, line))
    {
      // Get the parameter of each board or each detector
      std::istringstream iss(line);
      int label; iss >> label;
      double tmp_d;
      iss >> tmp_d; sShifts    .emplace(label, tmp_d);
      iss >> tmp_d; sThresholds.emplace(label, tmp_d);
      iss >> tmp_d; sFractions .emplace(label, tmp_d);
    }
  }
  
protected:
  double interpolate0(size_t const & bin) const noexcept 
  {
    auto const & cfd_0 = cfd[bin];
    auto const & cfd_1 = cfd[bin+1];
    if( cfd_0 == cfd_1) return 0;
    else return bin - cfd_0 / (cfd_1 - cfd_0);
  }
  size_t m_size = 0;

  
public:
  // TODO Parameters :

  using Shifts     = std::unordered_map<int, int   >;
  using Thresholds = std::unordered_map<int, double>;
  using Fractions  = std::unordered_map<int, double>;

  static inline Shifts     sShifts     = {};
  static inline Thresholds sThresholds = {};
  static inline Fractions  sFractions  = {};
  
  class Param
  {
  public:
    enum Type {BOARD, LABEL, UNDEFINED};
    static inline int sType = UNDEFINED;
    static inline constexpr void setType(int type) noexcept {sType = type;}
  };
  
  
};



#endif //CFD_HPP
