#pragma once

#include "../Colib/lib/libCo.hpp"
#include "../CaenLib/Hit.hpp"
#include "TraceAnalysis.hpp"
#include "CFD.hpp"

#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph2D.h"

namespace Caen1725
{
  class CFDMinimisationParameters
  {
  public:
    CFDMinimisationParameters() noexcept = default;
    CFDMinimisationParameters(std::string const & filename) noexcept {set(filename);}

    struct Fractions
    {
      std::vector<double> values; double first; double last; int nb_steps; double delta;
      void set()
      {
        nb_steps = 1 + (last-first)/delta;
        values.clear(); values.reserve(nb_steps);
        for (int i = 0; i<nb_steps; ++i) values.push_back(first+i*delta);
      }
    } fractions;

    struct Shifts
    {
      std::vector<int> values; int first; int last; int nb_steps;
      void set()
      {
        nb_steps = last - first + 1;
        values.clear();values.reserve(nb_steps);
        for (int i = 0; i<nb_steps; ++i) values.push_back(first+i);
      }
    } shifts;

    void set(std::string const & filename)
    {
      if (!Colib::fileExists(filename)) Colib::throw_error(filename+" not found !!");
      std::ifstream file(filename);
      std::string line;
      while(std::getline(file, line))
      {
        std::istringstream iss(line);
        std::string temp; 
        while(iss >> temp)
        {
          if (temp == "shifts") iss >> shifts.first >> shifts.last;
          else if (temp == "fractions") iss >> fractions.first >> fractions.last >> fractions.delta;
        }
      }
      
      fractions.set();
      shifts.set();
      print("fractions", fractions.values);
      print("shifts", shifts.values);
    }
    friend std::ostream& operator<< (std::ostream& out, CFDMinimisationParameters const & params)
    {
      out << "fractions " << params.fractions.first << " " << params.fractions.last << " " << params.fractions.nb_steps << "\n"
          << "shifts " << params.shifts.first << " " << params.shifts.last << " " << params.shifts.nb_steps << "\n";
      return out;
    }
  };

  class OptimizerHistograms
  {
  public:
    OptimizerHistograms() noexcept = default;

    void init(std::string const & name, CFDMinimisationParameters const & parameters)
    {
      auto const & fracv = parameters.fractions.values;
      auto const & shiftv = parameters.shifts.values;
      auto const & frac = parameters.fractions;
      auto const & shift = parameters.shifts;
      resolution_histos.reset(new TH2F(("resolution"+name).c_str(), ("resolution"+name+";fraction;shift").c_str(), 
        frac.nb_steps,frac.first-frac.delta/2,frac.last+frac.delta/2, shift.nb_steps,shift.first,shift.last+1)); 
      dT_histos.reset(new TH3F(("dT"+name).c_str(), ("dT"+name).c_str(), 
        frac.nb_steps, frac.first-frac.delta/2, frac.last+frac.delta/2,  shift.nb_steps,shift.first,shift.last+1, 4000,-2e6,2e6)); 
      resolution_histos -> SetDirectory(nullptr);
      dT_histos -> SetDirectory(nullptr);
    }

    void write(std::string rootFilename)
    {
      auto file = TFile::Open(rootFilename.c_str(), "recreate");
      if (!file) Colib::throw_error(rootFilename+" not created !! Is the path ok ?");
      dT_histos -> Write();
      resolution_histos -> Write();
      file->Close();
    }

    void write()
    {
      dT_histos -> Write();
      resolution_histos -> Write();
    }

    template<class... ARGS> void fill_dT  (double fraction, double shift, double dT) {dT_histos->Fill(fraction, shift, dT);}
    
    int pointNb = 0;
    std::unique_ptr<TH2F> resolution_histos;
    std::unique_ptr<TH3F> dT_histos;
  };

  class CFDOptimizer
  {
  public:
    CFDOptimizer(std::vector<int> const & labels) noexcept :
      m_nbDetectors(labels.size()),
      m_listLabels(labels) //eg {0,2,4,6,8}
    {
      size_t labelMax = *std::max_element(m_listLabels.begin(), m_listLabels.end());
      m_labelToDetectorIndex.resize(labelMax+1, -1);
      for (size_t det_i = 0; det_i<m_listLabels.size(); ++det_i) 
      {
        auto const & label = m_listLabels[det_i];
        m_labelToDetectorIndex[label] = det_i; // eg {0,-1,1,-1,2,-1,3,-1,4}
      }
    }

    void setParameters(CFDMinimisationParameters const & parameters) 
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

    double calculateResolution(TH1* histo)
    {
      if (!histo || histo->IsZombie() || histo->GetEntries() < 1) return 1e42;

      // Raw ROOT estimates in the range
      double max = histo->GetMaximum();
      double mean_est = histo->GetMean();
      double sigma_est = histo->GetStdDev();

      double bin_min = histo->FindFirstBinAbove(max*0.25);
      double bin_max = histo->FindLastBinAbove(max*0.25);

      if (bin_min < 0 || bin_max < 0) return 1e42;
      if (bin_min == bin_max) {bin_min-=1; bin_max+=1;}

      double x_min = histo->GetBinLowEdge(bin_min);
      double x_max = histo->GetBinLowEdge(bin_max);

      double FWHM = x_max - x_min;

      // return FWHM/1000.;

      // Instanciate the fit
      TF1* fit = new TF1("fit", "gaus", x_min, x_max);
      // Initialise the parameters
      fit->SetParameters(max, mean_est, FWHM/2.35);
      // Fit
      histo->Fit(fit, "RQ");
      // Get this first estimate
      double resolution = fit->GetParameter(2)*2.35/1000.;

      return resolution;
    }

    void calculateResolutions()
    {
      for (auto & label : m_listLabels)
      {
        auto const & index = m_labelToDetectorIndex[label];
        auto & dTs = m_histograms[index].dT_histos;
        for (int binx = 1; binx<=dTs->GetNbinsX(); ++binx) for (int biny = 1; biny<=dTs->GetNbinsY(); ++biny)
        {
          std::string name = std::to_string(binx)+std::to_string(biny);
          auto dT = dTs->ProjectionZ(name.c_str(), binx, binx, biny, biny);
          auto const resolution = calculateResolution(dT);
          m_histograms[index].resolution_histos->SetBinContent(binx, biny, resolution);
        }
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
    CFDMinimisationParameters m_parameters;
    std::vector<OptimizerHistograms> m_histograms;
  };
}