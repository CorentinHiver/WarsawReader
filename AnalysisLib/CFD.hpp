#ifndef CFD_HPP
#define CFD_HPP

#include "../LibCo/libCo.hpp"

class CFD
{
  using CFD_t = std::vector<double>;
public:

  template<class T = double>
  CFD(std::vector<T> const & _trace, size_t nb_samples = 40) :
    m_size(_trace.size())
  {
    if (_trace.empty()) throw_error("In CFD::checkTrace() : trace is empty");

    double baseline = 0;
    for (size_t sample_i = 0; sample_i<nb_samples; ++sample_i) baseline += _trace[sample_i];
    baseline /= nb_samples;
    
    trace.reserve(m_size);

    for (auto const & sample : _trace) 
    {
      if constexpr (is_signed<T>()) {
        trace.push_back(sample - baseline);
      }
      else {
        trace.push_back(double_cast(sample  + randomCo::fast_uniform()) - baseline);
      }
    }
  }

  template<class T = double>
  CFD(std::vector<T> const & _trace, int const & shift, double const & fraction) : CFD(_trace)
  {
    this -> calculate(shift, fraction);
  }

  void calculate(size_t const & shift, double const & fraction)
  {
    if (fraction>1.) {error("in CFD(trace, shift, fraction): fraction>1 !!"); return;}

    cfd.clear();
    cfd.reserve(m_size - 2*shift);
    // for (size_t bin = 0; bin<shift; ++bin) cfd.push_back(0);
    for (size_t bin = 2*shift; bin<m_size - shift; ++bin){
      auto const & value = fraction * trace[bin] - trace[bin - shift] ;
      cfd.push_back(value);
    }
  }

  double findZero(double const & threshold)
  {
    if (threshold>0) error("in CFD::findZero(threshold) : threshold > 0 !");
    for (size_t bin_i = 0; bin_i<cfd.size(); ++bin_i){ // Loop through the cfd values
      if (cfd[bin_i] < threshold){ // The cfd value crossed the threshold
        for (size_t bin_j = bin_i; bin_j>0; --bin_j){ // Looping back for looking for the zero crossing
          if (cfd[bin_j] > 0) return interpolate0(bin_j); // Zero crossing found, return the interpolated zero crossing between samples before and after
        }
        return noZero; // The 0 crossing happened before the first sample, so impossible to determine it
      }
    }
    return noSignal; // The signal never crosses the threshold -> the signal is too small, the cfd parameters are wrong, or the threshold is too large
  }
  
  CFD_t trace;
  CFD_t cfd;
  constexpr static double noZero = 0.;
  constexpr static double noSignal = -1.;
  
private:
  double interpolate0(size_t const & bin) const noexcept 
  {
    auto const & cfd_0 = cfd[bin];
    auto const & cfd_1 = cfd[bin+1];
    if( cfd_0 == cfd_1) return 0;
    else return bin - cfd_0 / (cfd_1 - cfd_0);
  }
  size_t m_size = 0;
};

#endif //CFD_HPP
