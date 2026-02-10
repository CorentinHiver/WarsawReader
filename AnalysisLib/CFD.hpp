#ifndef CFD_HPP
#define CFD_HPP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <unordered_map>
#include <sstream>
#include <thread>

class CFD
{
protected:
  /////////////////////////////////////////
  // Some helper functions (from libCo): //
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

  // Random generation
  
  // virtual inline double random_fast_uniform() noexcept {
  //   static thread_local std::minstd_rand generator(std::random_device{}());
  //   static thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
  //   return distribution(generator);
  // }

  // inline double random_fast_uniform() noexcept {
  //   static thread_local std::mt19937_64 gen(123456789);
  //   return std::generate_canonical<double, 53>(gen);
  // }

  // inline double random_fast_uniform() noexcept 
  // {
  //   // 1. Use a simple LCG (Knuth's MMIX constants).
  //   // It has a state of just 8 bytes (fits in a register).
  //   static thread_local uint64_t state = []{
  //       // 2. Seed uniquely per thread!
  //       // Combine time and thread ID to ensure threads have different streams.
  //       auto seed = static_cast<uint64_t>(std::time(nullptr));
  //       seed ^= std::hash<std::thread::id>{}(std::this_thread::get_id());
  //       return seed;
  //   }();

  //   // 3. The LCG Algorithm: state = state * A + C
  //   state = state * 6364136223846793005ULL + 1442695040888963407ULL;

  //   // 4. Fast conversion to double [0, 1) without full division.
  //   // We take the top 53 bits (for double precision) and multiply by 2^-53.
  //   // 0x1.0p-53 is the hex literal for 1.0 / 2^53.
  //   return (state >> 11) * 0x1.0p-53;
  // }

  template<class T>
  inline double adc_to_double(T adc) noexcept
  {
      // // ------------------- Fast Lehmer64 (128-bit state) -------------------
      // static thread_local __uint128_t state1 = [] {
      //     uint64_t s = static_cast<uint64_t>(std::time(nullptr));
      //     s ^= std::hash<std::thread::id>{}(std::this_thread::get_id());
      //     return (__uint128_t)s << 64 | s;
      // }();

      // constexpr __uint128_t MUL = 0xda942042e4dd58b5ULL;   // passes BigCrush, very fast
      // state1 *= MUL;

      // uint64_t r1 = state1 >> 64;
      
      // constexpr double TO_DOUBLE = 0x1.0p-53; // Convert to [0,1) (top 53 bits)

      // return (static_cast<double>(adc) + (r1 >> 11) * TO_DOUBLE);
    return static_cast<double>(adc);
  }

  // Consol colors
  static constexpr const char* RED   = "\u001b[31m";
  static constexpr const char* RESET = "\u001b[0m" ;
  
public:
  /// @brief Default constructor
  CFD() noexcept = default; 
  
  template<class T = double>
  CFD(std::vector<T> && _trace, int shift, double fraction, size_t nb_samples_baseline = 10) :
    trace(std::move(_trace))
  {
    m_baseline = _trace[0];
    if (nb_samples_baseline > 1)
    {
      for (size_t sample_i = 1; sample_i<nb_samples_baseline; ++sample_i) m_baseline += _trace[sample_i];
      m_baseline /= nb_samples_baseline;
    }

    cfd.reserve(m_size);
    for (size_t bin = 5*shift; bin<m_size - shift; ++bin)
    {
      auto const & value = fraction * trace[bin] - trace[bin - shift] ;
      cfd.push_back(value);
    }
  }

  /** @brief Constructs a CFD object from a trace
  * @details 
  * Fills an internal vector with a trace after the m_baseline is determined 
  * as the mean value of the first nb_samples_baseline samples, 
  * and subtracted to all the samples of the trace. 
  * A random number in [0;1[ is added to integer values to smooth the values 
  **/
  template<class T = double>
  CFD(std::vector<T> const & _trace, size_t nb_samples_baseline = 10) :
    m_size(_trace.size())
  {
    setTrace(_trace, nb_samples_baseline);
  }

  template<class T = double>
  CFD& setTrace(std::vector<T> const & _trace, size_t nb_samples_baseline = 10)
  {    
    if (_trace.empty()) return *this;

    m_baseline = _trace[0];
    if (nb_samples_baseline > 1)
    {
      for (size_t sample_i = 1; sample_i<nb_samples_baseline; ++sample_i) m_baseline += _trace[sample_i];
      m_baseline /= nb_samples_baseline;
    }
    
    trace.reserve(m_size);

    for (auto const & sample : _trace)
    {
      if constexpr (is_floating<T>()) trace.push_back(sample - m_baseline);
      else trace.push_back(adc_to_double(sample) - m_baseline);
    }
    return *this;
  }

  template<class T = double>
  CFD(std::vector<T> const & _trace, int shift, double fraction, size_t nb_samples_baseline = 10) : CFD(_trace, nb_samples_baseline)
  {
    this -> calculate(shift, fraction);
  }

  virtual void calculate(size_t shift, double fraction)
  {
    if (fraction>1.) {std::cout << RED << "in CFD(trace, shift, fraction): fraction>1 !!" << RESET << std::endl; return;}

    cfd.clear();

    if (m_size < 2*shift) {std::cout << RED << "in CFD(trace, shift, fraction): m_size = " << m_size << " < 2*shift = " << 2*shift << " !!" << RESET << std::endl; return;}
    
    cfd.reserve(m_size);
    for (size_t bin = 5*shift; bin<m_size - shift; ++bin)
    {
      auto const & value = fraction * trace[bin] - trace[bin - shift] ;
      cfd.push_back(value);
    }
  }

  /// @brief Calculates the last zero crossing before the calculated cfd signal goes below the given threshold
  double findZero(double threshold)
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

  /// @brief Finds the last zero crossing before the cfd trace reaches its minimum
  double findZero()
  {
    auto minimum_it = std::min_element(std::begin(cfd), std::end(cfd));
    auto const & minimum = *minimum_it;
    if (0 < minimum) return noSignal;            // If never crosses zero, returns noSignal
    auto const & min_bin = std::distance(std::begin(cfd), minimum_it); // Get the minimum bin number

    for (size_t bin_j = min_bin; bin_j>0; --bin_j)    // Looping back to look for the zero crossing
      if (cfd[bin_j] > 0) return interpolate0(bin_j); // Zero crossing found, return the interpolated zero crossing between samples before and after
    return noZero; // The 0 crossing happened before the first sample, so impossible to determine it
  }
  
  CFD_t trace;
  CFD_t cfd;

  /////////////////////////
  // Parameters handling //
  /////////////////////////

  // Static variables :
  constexpr static double noZero   = 1e-100;
  constexpr static double noSignal = 1e-101;

  static void loadParameters(std::string const & filename)
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

  size_t m_size = 0;
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


    // if (0 < minimum(cfd)) return noSignal;            // If never crosses zero, returns noSignal
    // auto const & min_bin = minimum_index(cfd);        // Get the minimum bin number

    
  // template<class T = double>
  // CFD(std::vector<T> const & _trace, size_t shift, double fraction, size_t nb_samples_baseline = 10) : 
  //   m_size(_trace.size())
  // {
  //   m_baseline = _trace[0];
  //   if (nb_samples_baseline > 1)
  //   {
  //     for (size_t sample_i = 1; sample_i<nb_samples_baseline; ++sample_i) m_baseline += _trace[sample_i];
  //     m_baseline /= nb_samples_baseline;
  //   }

  //   cfd.clear();

  //   if (m_size < 2*shift) {std::cout << RED << "in CFD(trace, shift, fraction): m_size = " << m_size << " < 2*shift = " << 2*shift << " !!" << RESET << std::endl; return;}
    
  //   cfd.reserve(m_size);
  //   for (size_t bin = 5*shift; bin<m_size - shift; ++bin)
  //   {
  //     auto const & value = fraction * static_cast<float>(_trace[bin]+random_fast_uniform() - _trace[bin - shift]);
  //     cfd.push_back(value);
  //   }
  // }