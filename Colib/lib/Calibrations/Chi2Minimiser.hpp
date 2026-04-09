#pragma once
#include "TH2F.h"
#include "TH3F.h"

#include "../libCo.hpp"
#include "FlexibleHisto.hpp"

namespace Colib
{
  class Chi2Calculator
  {
  public:
    Chi2Calculator(TH1* reference) : m_reference(reference) {}

    ~Chi2Calculator() {delete m_histo_flex;}
    
    template<class THist>
    double operator()(THist* const test, FlexibleHisto const & histo)
    {
      auto testCal = histo(test);
      testCal -> SetName((test -> GetName() + std::string("calib : ") + mergeStrings(histo.getCalib().get(), "_")).c_str());
      return calculate(testCal);
    }
    
    /// @brief Calculates the chi2 between the reference histogram and a given histogram
    template<class THist> double calculate(THist* histo, bool normalize=false) const
    {
      if (normalize) normalise(histo);
      double sum_errors_squared = 0.0;
      auto bin_min = histo->FindBin(m_bounds.first);
      auto bin_max = histo -> GetNbinsX();

      if (0 < m_bounds.second)
      {
        auto const test_bin_max = histo->FindBin(m_bounds.second);
        if (test_bin_max < bin_max) bin_max = test_bin_max;
      }

      for (int bin = bin_min; bin<bin_max; ++bin) if (histo -> GetBinContent(bin)>0)
      {
        // Calculate the difference for this bin :
        auto const diff = m_reference -> GetBinContent(bin) - histo -> GetBinContent(bin);

        // Variance of the bin :
        double const weight = 1 / histo -> GetBinContent(bin); // V = sigma² = 1/N

        // Add the diff to the total squared diff of the spectra :
        sum_errors_squared += diff*diff*weight;

      }
      return sum_errors_squared/(bin_max - bin_min);
    }

    /// @brief Calculates the chi2 between the reference histogram and a given flexible histogram
    double calculate(FlexibleHisto const & histo_flex) const
    {
      auto histo = histo_flex.getHisto();

      double sum_errors_squared = 0.0;
      auto bin_min = histo->FindBin(m_bounds.first);
      auto bin_max = histo -> GetNbinsX();

      if (0 < m_bounds.second)
      {
        auto test_bin_max = histo->FindBin(m_bounds.second);
        if (test_bin_max < bin_max) bin_max = test_bin_max;
      }

      // std::vector<double> values;
      double constexpr tolerance = 0.005;
      double constexpr tolerance_2 = tolerance*tolerance;
      int nbBins = 0;
      for (int bin = bin_min; bin<bin_max; bin++) 
      {
        ++nbBins;
        auto const refValue = m_reference -> GetBinContent(bin);
        auto const testValue = histo_flex[bin];
        if (refValue < 5 || testValue < 5) continue;

        // Calculate the difference for this bin :
        auto const diff = refValue - testValue;
  
        // Variance of the bin :
        double const weight = 1/(testValue + testValue*testValue*tolerance_2); // 1/(V+tol²)
        // double const weight = 1/testValue; // 1/V = 1/sigma² = 1/N
  
        // Add the diff to the total squared diff of the spectra :
        sum_errors_squared += diff*diff*weight;
        // values.push_back(diff*diff*weight);
        // if (diff*diff*weight>1e8) print("Chi2 = ", diff*diff*weight); // work in progress
      }
      // print(*std::max_element(values.begin(), values.end()));
      // if (sum_errors_squared == 0) print("sum_errors_squared = 0");
      return sum_errors_squared/(nbBins);
    }

    double calculateForMinuit(double const * par)
    {
      double sum_errors_squared = 0.0;

      auto histo_flex = *m_histo_flex; // Aliasing
      histo_flex.setCalibAndScale({par[0], par[1], par[2]});

      auto const bins = histo_flex.getHisto() -> GetNbinsX();
  
      for (int bin = 0; bin<bins; bin++) if (histo_flex[bin]>0)
      {
        // Calculate the difference for this bin :
        auto const diff = m_reference -> GetBinContent(bin)-histo_flex[bin];
  
        // Variance of the bin :
        double const weight = 1/histo_flex[bin]; // V = sigma² = 1/N
  
        // Add the diff to the total squared diff of the spectra :
        sum_errors_squared += diff*diff*weight;
  
      }
      return sum_errors_squared/bins;
    }

    void setCalibForMinuit(FlexibleHisto * histo_flex)
    {
      m_histo_flex = histo_flex;
    }
    
    void setBounds(std::pair<int, int> bounds) {m_bounds = bounds;}
    void setBounds(int bound_min, int bound_max) {m_bounds = {bound_min, bound_max};}

    void normalise(TH1* testHisto) const
    {
      if (!m_reference) Colib::throw_error("in Chi2Minimizer::testHisto(test) : no reference !!");
      if (!testHisto) Colib::throw_error("in Chi2Minimizer::testHisto(test) : no test histogram !!");
      double integral_ref  = m_reference->Integral();
      double integral_test = testHisto->Integral();

      if (integral_test != 0) {
          double scale_factor = integral_ref / integral_test;
          testHisto->Scale(scale_factor);
      }
    }

  private:
    std::mutex m_mutex;
    FlexibleHisto * m_histo_flex = nullptr; // Used only for the minuit interface
    TH1* m_reference = nullptr;
    std::pair<int, int> m_bounds = {0,-1};
  };
  
  struct MinimiserVariable
  {
    double initGuess = 0.;
    double bound = 0.;
    double step = 0.;
    int nb_steps = 0;
    double min = 0; 
    double max = 0; 
  
    MinimiserVariable(double _initGuess, double _bound, double _nb_steps)
    {
      initGuess =_initGuess;
      bound =_bound;
      nb_steps =_nb_steps;
      initialize();
    }

    MinimiserVariable(std::initializer_list<double> init)
    {
      auto it = init.begin();
      initGuess = double_cast(*it++);
      bound = double_cast(*it++);
      nb_steps = int_cast(*it++);
      initialize();
    }
  
    MinimiserVariable& operator=(std::initializer_list<double> init)
    {
      auto it = init.begin();
      initGuess = double_cast(*it++);
      bound = double_cast(*it++);
      nb_steps = int_cast(*it++);
      initialize();
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& out, MinimiserVariable const & minvar)
    {
      out << 
            " initGuess " << minvar.initGuess <<
            " bound "     << minvar.bound     <<
            " step "      << minvar.step      <<
            " nb_steps "  << minvar.nb_steps  <<
            " min "       << minvar.min       <<
            " max "       << minvar.max       ;
      return out ;
    }
  
  private:

    void initialize() 
    {
      step = bound / nb_steps;
      min = initGuess - nb_steps * step;
      max = initGuess + nb_steps * step;
    }
  };
  
  class Chi2Minimiser
  {
  public:
    Chi2Minimiser(){}
    
    void calculate(Chi2Calculator & chi2Calc, FlexibleHisto const & testHisto)
    {
      auto chi2 = chi2Calc.calculate(testHisto);
      if (chi2 == 0) print("chi2 = 0");
      if (chi2 < m_min_chi2) 
      {
        m_min_chi2 = chi2;
        m_histo_flex = testHisto;
      }
      if (s_fill_histo && m_chi2map) 
      {
        auto parameters = testHisto.getCalib().get();
        m_chi2map->Fill(parameters[0], parameters[1], parameters[2], chi2);
      }
    }

    template<class THist>
    void minimise(Chi2Calculator & chi2Calc, THist* test, MinimiserVariable xParam, MinimiserVariable yParam, MinimiserVariable zParam, bool normalise = false)
    {
      if (normalise) chi2Calc.normalise(test);
      FlexibleHisto histo_flex(test);
      m_histo_flex.setHisto(test);
      if (m_bruteforce)
      {
        TString chi2map_name = "chi2map";
      #ifdef COMULTITHREADING
        chi2map_name += MTObject::getThreadIndex();
      #endif //COMULTITHREADING
        if (s_fill_histo && !m_chi2map) 
        {
          m_chi2map .reset(new TH3F(chi2map_name, "chi2map;x;y;z", 
                            xParam.nb_steps * 2, xParam.min, xParam.max * (1 + 1e-5), 
                            yParam.nb_steps * 2, yParam.min, yParam.max * (1 + 1e-5),
                            zParam.nb_steps * 2, zParam.min, zParam.max * (1 + 1e-5)
                          ));
          m_chi2map -> SetDirectory(nullptr);
        }

        for (int stepx = 0; stepx<xParam.nb_steps*2; ++stepx) {
          for (int stepy = 0; stepy<yParam.nb_steps*2; ++stepy) {
            for (int stepz = 0; stepz<zParam.nb_steps*2; ++stepz) {

              // 1. Set the calibration
              histo_flex.setCalibAndScale({
                xParam.min + stepx*xParam.step,
                yParam.min + stepy*yParam.step,
                zParam.min + stepz*zParam.step
              });
              
              // 2. Calculate the chi2 between the reference spectrum and the calibrated and scaled histogram
              this -> calculate(chi2Calc, histo_flex);
            }
          }
        }
      }
      else
      {
      #ifdef INCLUDE_MINUIT // Work in progress
        // chi2Calc.setCalibForMinuit(&calib);
        // ROOT::Math::Minimizer* minimizer = new ROOT::Minuit2::Minuit2Minimizer("Migrad") ;
        auto minimizer = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);
        // ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        
        // Set properties
        minimizer->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
        // minimizer->SetMaxIterations(10000);     // for GSL
        minimizer->SetTolerance(0.001);
        minimizer->SetPrintLevel(1);            // 0: silent, 1: default
        
  
        // Create Functor object
        auto func = [&chi2Calc](const double* x) {
          return chi2Calc.calculateForMinuit(x);
        };
        ROOT::Math::Functor f(func, 3);
  
        // Set function and initial values
        minimizer->SetFunction(f);
        minimizer->SetVariable(0, "param_x0", 0.0, xParam.step);
        minimizer->SetVariable(1, "param_x1", 1.0, yParam.step);
        minimizer->SetVariable(2, "param_x2", 1.0, zParam.step);
        
        // Perform minimization
        minimizer->Minimize();
  
        // Get results
        const double* xs = minimizer->X();
        m_histo_flex = {xs[0], xs[1], xs[2]};
        // std::cout << "Minimum at: x0 = " << xs[0] << ", x1 = " << xs[1] << ", x2 = " << xs[2] << std::endl;
        // std::cout << "Minimum value: f = " << minimizer->MinValue() << std::endl;
      #else 
        Colib::throw_error("Compile with -DINCLUDE_MINUIT");
      #endif //INCLUDE_MINUIT
      }
    }
  
    auto const & getCalib() const {return m_histo_flex.getCalib();}
    auto const & getHistoTest() const {return m_histo_flex;}
    auto const & getMinChi2() const {return m_min_chi2;}
    auto getChi2Map() const {return m_chi2map.get();}
  
    void brutefore(bool b) {m_bruteforce = b;}
    void multistages(int nb_stages)
    {
      m_nb_stages = nb_stages;
      m_multistages = true;
    }
  
    static void fillHisto(bool b) {s_fill_histo = b;}
  
  private:

    static bool s_fill_histo; 
    bool m_bruteforce = true;
    bool m_multistages = false;
    int m_nb_stages = 1;
    double m_min_chi2 = 1e100;
    
    FlexibleHisto m_histo_flex;
    std::unique_ptr<TH3F> m_chi2map;
  };
  
  bool Chi2Minimiser::s_fill_histo = true; 
}