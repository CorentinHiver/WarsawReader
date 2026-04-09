#pragma once

#include "libCo.hpp"
#include "TH2F.h"
#include "TH3F.h"
#include "../CaenLib/RootHit.hpp"
#include "TraceAnalysis.hpp"
#include "CFD.hpp"

namespace Caen1725
{
  struct CFDParameter
  {
    double first{};
    double last{};
    int nb_steps{};

    std::vector<double> values;
    void set() noexcept
    {
      values.clear();
      values.reserve(nb_steps);
      auto const delta = (last-first)/nb_steps;
      for (int i = 0; i<nb_steps; ++i) values.push_back(first+i*delta);
    }
    friend std::ostream& operator<< (std::ostream& out, CFDParameter const & param)
    {
      out << "first " << param.first << " last " << param.last << " nb_steps " << param.nb_steps;
      return out; 
    }
  };

  class CFDParameters
  {
  public:
    CFDParameter fractions;
    CFDParameter shifts;

    void set(std::string filename)
    {
      if (!Colib::fileExists(filename)) Colib::throw_error(filename+" not found !!");
      std::ifstream file(filename);
      std::string line;
      while(std::getline(file, line))
      {
        print(line);
        std::istringstream iss(line);
        std::string temp; 
        while(iss >> temp)
        {
          print(temp);
          if (temp == "shifts") iss >> shifts.first >> shifts.last >> shifts.nb_steps;
          else if (temp == "fractions") iss >> fractions.first >> fractions.last >> fractions.nb_steps;
        }
      }
      fractions.set();
      shifts.set();
    }
    friend std::ostream& operator<< (std::ostream& out, CFDParameters const & params)
    {
      out << "fractions " << params.fractions << "\nshifts " << params.shifts << "\n";
      return out;
    }
  };

  class OptimizerHistograms
  {
  public:
    OptimizerHistograms() noexcept = default;

    void init(std::string const & name, CFDParameters const & parameters)
    {
      auto const & frac = parameters.fractions;
      auto const & shift = parameters.shifts;
      print(parameters);
      resolution_histos.reset(new TH2F((name+"_resolution").c_str(), (name+"_resolution").c_str(), frac.nb_steps, frac.first, frac.last,  shift.nb_steps, shift.first, shift.last)); 
      dT_histos.reset(new TH3F((name+"_dT").c_str(), (name+"_dT").c_str(), frac.nb_steps, frac.first, frac.last,  shift.nb_steps, shift.first, shift.last, 4000,-2e6,2e6)); 
      resolution_histos -> SetDirectory(nullptr);
      dT_histos -> SetDirectory(nullptr);
    }

    void write(std::string rootFilename)
    {
      auto file = TFile::Open(rootFilename.c_str(), "recreate");
      if (!file) Colib::throw_error(rootFilename+" not created !! Is the path ok ?");
      resolution_histos -> Write();
      dT_histos -> Write();
      file->Close();
    }

    void write()
    {
      resolution_histos -> Write();
      dT_histos -> Write();
    }

    template<class... ARGS> void fill_dT  (ARGS &&... args) {dT_histos  ->Fill(std::forward<ARGS>(args)...);}
    template<class... ARGS> void fill_resolution(ARGS &&... args) {resolution_histos->Fill(std::forward<ARGS>(args)...);}
    
    std::unique_ptr<TH2F> resolution_histos;
    std::unique_ptr<TH3F> dT_histos;
  };

  class CFDOptimizer
  {
  public:
    CFDOptimizer(std::vector<int> labels) noexcept :
      m_nbDetectors(labels.size()),
      m_listLabels(labels) //eg {0,2,4,6,8}
    {
      size_t labelMax = *std::max_element(m_listLabels.begin(), m_listLabels.end());
      m_labelToDetectorIndex.resize(labelMax, -1);
      for (size_t det_i = 0; det_i<m_listLabels.size(); ++det_i) 
      {
        auto const & label = m_listLabels[det_i];
        m_labelToDetectorIndex[label] = det_i; // eg {0,-1,1,-1,2,-1,3,-1,4}
      }
    }

    void setParameters(CFDParameters const & parameters) 
    {
      m_parameters = parameters;
      m_histograms.resize(m_nbDetectors);
      for (size_t det_i = 0; det_i<m_nbDetectors; ++det_i) 
      {
        auto const & label = m_listLabels[det_i];
        m_histograms[det_i].init(std::to_string(label), m_parameters);
      }
    }

    void calculate_dT(Hit const & refHit, Label label, Timestamp time, Trace const & trace) 
    {
      static thread_local CFD cfd;
      for (auto const & shift : m_parameters.shifts.values) for (auto const & fraction : m_parameters.fractions.values)
      {
        if (m_labelToDetectorIndex.size() <= label) continue; // Detector non treated
        auto const & index = m_labelToDetectorIndex[label];
        if (index < 0) continue; // Detector non treated

        cfd.generate(trace, shift, fraction, 10);

        auto const zero = cfd.findZero();
        if (zero == CFD::noSignal || zero == CFD::noZero) continue;
        auto const time_cfd = time + zero*ticks_to_ps;

        m_histograms[index].fill_dT(fraction, shift, refHit.time - time_cfd);
      }
    }

    void calculateResolution(TH1* histo)
    {
      auto max = histo->GetMaximum();
      auto maxBin = histo->GetMaximumBin();
      auto leftBin = histo->FindFirstBinAbove(max*0.7);
      auto rightBin = histo->FindLastBinAbove(max*0.7);
      histo->GetXaxis()->SetRange(leftBin*2, rightBin*2);

      // Raw ROOT estimates in the range
      double mean_est = histo->GetMean();
      double sigma_est = histo->GetRMS();
      double height_est = histo->GetMaximum();

      // Set fit range
      double x_min = mean_est - 3 * sigma_est;
      double x_max = mean_est + 3 * sigma_est;

      // Instanciate the fit
      TF1* fit = new TF1("fit", "gaus", x_min, x_max);
      // Initialise the parameters
      fit->SetParameters(height_est, mean_est, sigma_est);
      // Fit
      histo->Fit(fit, "RQ");
      // Get this first estimate
      double final_mean  = gFit->GetParameter(1);
      double final_sigma = gFit->GetParameter(2);
    }

    void calculateResolutions()
    {
      for (auto & label : m_listLabels)
      {
        auto const & index = m_labelToDetectorIndex[label];
        auto & dTs = m_histograms[index].dT_histos;
        for (int binx = 1; binx<=dTs->GetNbinsX(); ++binx) for (int biny = 1; biny<=dTs->GetNbinsX(); ++biny)
        {
          std::string name = std::to_string(binx)+std::to_string(biny);
          auto dT = dTs->ProjectionZ(name.c_str(), binx, binx, biny, biny);
          m_histograms[index].resolution_histos->SetBinContent(binx, biny, calculate_resolution(dT));
        }
        m_histograms[index].resolution_histos;
      }
    }

    void write(std::string filename)
    {
      auto file = TFile::Open(filename.c_str(), "recreate");
      if (!file) Colib::throw_error(filename+" not created !! Is the path ok ?");
      for (auto & histos : m_histograms) histos.write();
      file->Close();
      print(filename, "written");
    }
    
  private:
    size_t m_nbDetectors = 0;
    std::vector<int> m_listLabels;
    std::vector<int> m_labelToDetectorIndex;
    CFDParameters m_parameters;
    std::vector<OptimizerHistograms> m_histograms;
  };
}