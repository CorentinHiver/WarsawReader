#pragma once

#include "../Colib/lib/libRoot.hpp"
#include "../CaenLib/Hit.hpp"
// #include "TraceAnalysis.hpp"
#include "CFD.hpp"

namespace Caen1725
{
  struct CFDMinimisationParameterI
  {
    std::vector<int> values; int first; int last; int nb_steps; int delta;
    void set()
    {
      nb_steps = last - first + 1;
      delta = 1;
      values.clear();values.reserve(nb_steps);
      for (int i = 0; i<nb_steps; ++i) values.push_back(first+i);
    }
  };
  struct CFDMinimisationParameterD
  {
    std::vector<double> values; double first; double last; int nb_steps; double delta;
    void set()
    {
      nb_steps = 1 + (last-first)/delta;
      values.clear(); values.reserve(nb_steps);
      for (int i = 0; i<nb_steps; ++i) values.push_back(first+i*delta);
    }
  };
  class CFDMinimisationParameters
  {
    // std::array<bool, 10000> labelLUT{{}};
    std::vector<Label> labels;
  public:
    CFDMinimisationParameters() noexcept = default;
    CFDMinimisationParameters(std::string const & filename) noexcept {set(filename);}

    CFDMinimisationParameterD fractions;
    CFDMinimisationParameterI shifts;

    int nbBaseline = 10;

    void set(std::string const & filename)
    {
      if (!Colib::fileExists(filename)) Colib::throw_error("cfd parameter file", filename, "not found !!");
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
          else if (temp == "nbBaseline") iss >> nbBaseline;
          else if (temp == "labels")
          {
            size_t label{};
            while(iss >> label) 
            {
              // if (labelLUT.size() < label) Colib::throw_error("Label", label, "too high (<", labelLUT.size(), ")");
              labels.push_back(label);
              // labelLUT[label] = true;
            }
          }
          else if (temp == "boards")
          {
            size_t board{};
            while(iss >> board) for (size_t channel = 0; channel<16; ++channel) 
            {
              auto label = 16*board+channel;
              // if (labelLUT.size() < label) Colib::throw_error("Label", label, "too high (<", labelLUT.size(), ")");
              labels.push_back(label);
              // labelLUT[label] = true;
            }
          }
        }
      }
      
      fractions.set();
      shifts.set();
      // print("fractions", fractions.values);
      // print("shifts", shifts.values);
    }
    auto const & getLabels() const noexcept {return labels;}
    friend std::ostream& operator<< (std::ostream& out, CFDMinimisationParameters const & params)
    {
      out << "fractions " << params.fractions.first << " " << params.fractions.last << " " << params.fractions.nb_steps << "\n"
          << "shifts " << params.shifts.first << " " << params.shifts.last << " " << params.shifts.nb_steps << "\n"
          << "nbBaseline " << params.nbBaseline;
      out << "labels ";
      for (auto const & label : params.labels) out << label << " ";
      return out;
    }
  };

  class OptimizerHistograms
  {
  public:
    OptimizerHistograms() noexcept = default;

    void init(std::string const & name, CFDMinimisationParameters const & parameters, int nb_dT_points = 4000, double min_dT = -2e6, double max_dT = 2e6)
    {
      auto const & frac = parameters.fractions;
      auto const & shift = parameters.shifts;
      resolution_histos.reset(new TH2F(("resolution"+name).c_str(), ("resolution"+name+";fraction;shift").c_str(), 
        frac.nb_steps,frac.first,frac.last+frac.delta, shift.nb_steps,shift.first,shift.last+shift.delta)); 
      dT_histos.reset(new TH3F(("dT"+name).c_str(), ("dT"+name).c_str(), 
        frac.nb_steps,frac.first,frac.last+frac.delta,  shift.nb_steps,shift.first,shift.last+shift.delta, nb_dT_points,min_dT,max_dT)); 
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
    template<class... ARGS> void fillbin_dT  (int fraction_bin, int shift_bin, double dT) 
    {
      auto dT_bin = dT_histos->GetZaxis()->FindBin(dT);
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6, 8, 0)
      if (0 < dT_bin && dT_bin < dT_histos->GetNbinsZ()) dT_histos->AddBinContent(dT_histos->GetBin(fraction_bin+1, shift_bin+1, dT_bin+1));
    #else
      if (0 < dT_bin && dT_bin < dT_histos->GetNbinsZ()) dT_histos->AddBinContent(fraction_bin+1, shift_bin+1, dT_bin+1);
    #endif
    }

    std::unique_ptr<TH1D> get_dT(int fraction_bin, int shift_bin) const
    {
      auto dTName = dT_histos->GetName() + std::to_string(fraction_bin) + std::to_string(shift_bin);  
      return std::make_unique<TH1D>(*(dT_histos->ProjectionZ(dTName.c_str(), fraction_bin+1, shift_bin+1)));
    }
    
    std::unique_ptr<TH2F> resolution_histos;
    std::unique_ptr<TH3F> dT_histos;
  };

  class CFDOptimizer
  {
    std::array<int, 10000> labelLUT{{}};
  public:
    CFDOptimizer() noexcept = default;
    CFDOptimizer(CFDMinimisationParameters const & parameters):
      m_parameters(parameters),
      m_nbDetectors(m_parameters.getLabels().size())
    {
      size_t labelMax = *std::max_element(m_parameters.getLabels().begin(), m_parameters.getLabels().end());
      if (labelLUT.size() < labelMax) Colib::throw_error("Label", labelMax, "too high (>", labelLUT.size(), ")");
      m_histograms.resize(m_nbDetectors);
      labelLUT.fill(-1);
      for (size_t det_i = 0; det_i<m_parameters.getLabels().size(); ++det_i) 
      {
        auto const & label = m_parameters.getLabels()[det_i];
        // m_labelToDetectorIndex[label] = det_i;
        labelLUT[label] = det_i;
        m_histograms[det_i].init(std::to_string(label), m_parameters);
      }
      gaus.reset( new TF1("gaus", "gaus"));
      gaus_and_bkgd.reset( new TF1("gaus_and_bkgd", "gaus(0)+pol1(3)"));
    }

    void calculate_dT(Timestamp timeRef, Timestamp time, Label label, Trace const & trace) 
    {
      if (trace.empty() || labelLUT[label] < 0) return;
      auto const & index = labelLUT[label];

      static thread_local CFD cfd;
      cfd.setBaseline(trace, m_parameters.nbBaseline);
      for (int shift_i = 0; shift_i<m_parameters.shifts.nb_steps; ++shift_i) 
        for (int fraction_i = 0; fraction_i<m_parameters.fractions.nb_steps; ++fraction_i) 
      {

        cfd.generate(trace, m_parameters.shifts.values[shift_i], m_parameters.fractions.values[fraction_i]);

        auto const zero = cfd.findZero();
        if (zero == CFD::noSignal || zero == CFD::noZero) continue;
        auto const time_cfd = time + zero*ticks_to_ps;

        m_histograms[index].fillbin_dT(fraction_i, shift_i, timeRef - time_cfd);
        // m_histograms[index].fill_dT(fraction, shift, refHit.time - time_cfd);
      }
    }

    double calculateResolution(TH1* histo)
    {
      if (!histo || histo->IsZombie() || histo->GetEntries() < 1) return 1e42;

      double max = histo->GetMaximum();

      // Full Width at Quarter Maximum (FWQM) :
      double T = max * 0.25;
      double bin_min = histo->FindFirstBinAbove(T);
      double bin_max = histo->FindLastBinAbove(T);

      if (bin_min < 1 || bin_max < 1 || histo->GetNbinsX() <= bin_min || histo->GetNbinsX() <= bin_max) return 1e42;
      if (bin_min == bin_max) {bin_min-=1; bin_max+=1;}

      double x1 = histo->GetBinCenter(bin_min - 1);
      double y1 = histo->GetBinContent(bin_min - 1);
      double x2 = histo->GetBinCenter(bin_min);
      double y2 = histo->GetBinContent(bin_min);
      double x_min_interp = x1 + (T - y1) * (x2 - x1) / (y2 - y1);

      double x3 = histo->GetBinCenter(bin_max);
      double y3 = histo->GetBinContent(bin_max);
      double x4 = histo->GetBinCenter(bin_max + 1);
      double y4 = histo->GetBinContent(bin_max + 1);
      double x_max_interp = x3 + (T - y3) * (x4 - x3) / (y4 - y3);

      double FWQM = x_max_interp - x_min_interp;
      
      // Full Width at Half Maximum (FWHM) :

      T = max * 0.25;
      bin_min = histo->FindFirstBinAbove(T);
      bin_max = histo->FindLastBinAbove(T);

      if (bin_min < 1 || bin_max < 1 || histo->GetNbinsX() <= bin_min || histo->GetNbinsX() <= bin_max) return 1e42;
      if (bin_min == bin_max) {bin_min-=1; bin_max+=1;}

      x1 = histo->GetBinCenter(bin_min - 1);
      y1 = histo->GetBinContent(bin_min - 1);
      x2 = histo->GetBinCenter(bin_min);
      y2 = histo->GetBinContent(bin_min);
      x_min_interp = x1 + (T - y1) * (x2 - x1) / (y2 - y1);

      x3 = histo->GetBinCenter(bin_max);
      y3 = histo->GetBinContent(bin_max);
      x4 = histo->GetBinCenter(bin_max + 1);
      y4 = histo->GetBinContent(bin_max + 1);

      x_max_interp = x3 + (T - y3) * (x4 - x3) / (y4 - y3);

      double FWHM = x_max_interp - x_min_interp;
      double FWHM_from_FWQM = FWQM/sqrt(2);
      return (FWHM_from_FWQM+FWHM) / 2000;

      // Here lies previous attemps :
      {
        // Check if the difference is less than 10%:
        // if ( (FWHM - FWHM_from_FWQM) / ((FWHM + FWHM_from_FWQM) / 2) < 0.1) 

        // We measured FWHM and FWQM. If the peak is gaussian, FWQM = sqrt(2)*FWHM.
        // In case of bad CFD parameters, the peak won't be gaussian and a gaussian
        // gaussian fit is very likely to fail. Therefore, if FWHM is too different from 
        // sqrt(2)*FWQM then we don't try to fit it and FWQM/sqrt(2) is actually more
        // likely to be a good measurement of the FWHM.
        
        // auto max_binX = histo -> GetBinLowEdge(histo -> GetMaximumBin());
        // gaus->SetRange(max_binX-FWQM, max_binX+FWQM);
        // histo->GetXaxis()->SetRangeUser(max_binX-FWQM, max_binX+FWQM);

        // auto mean_est = histo -> GetMean();
        // gaus->SetParameters(max, mean_est, FWHM_from_FWQM);
        // // Fit
        // histo->Fit(gaus.get(), "RQ");
        // // Get this first estimate
        // double sigma = gaus->GetParameter(2);

        // // if (20 < FWHM) return FWHM;

        // // Initialise the parameters
        // // gaus_and_bkgd->SetRange(mean_est-FWHM, mean_est+FWHM);
        // // gaus_and_bkgd->SetParameters(max, mean_est, sigma, 0, 1);
        // // // Fit
        // // histo->Fit(gaus_and_bkgd.get(), "RQ");
        // // // Get this first estimate
        // // sigma = gaus_and_bkgd->GetParameter(2);

        // return Colib::sigtofwhm(sigma)/1000.;
      }
    }

    void calculateResolutions()
    {
      for (auto & label : m_parameters.getLabels())
      {
        // auto const & index = m_labelToDetectorIndex[label];
        auto const & index = labelLUT[label];
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

    void findMinima(std::string filename)
    {
      auto file = TFile::Open(filename.c_str(), "read");
      if (!file) Colib::throw_error(filename+" not found !!");

      {
        auto histos = Colib::file_get_map_of<TH2F>();
        for (auto const & [name, histo] : histos)
        {
          std::string label_str = name;
          Colib::remove(label_str, "resolution");
          int label = std::stoi(label_str);
          Colib::Smooth(histo, 2);
          auto const globalMin = histo->GetMinimumBin();
          int xbin, ybin, zbin;
          histo->GetBinXYZ(globalMin, xbin, ybin, zbin);
          m_minimaBin.emplace(label, std::array<int   , 2>({xbin, ybin}));
          m_minima   .emplace(label, std::array<double, 2>({histo->GetXaxis()->GetBinLowEdge(xbin), histo->GetYaxis()->GetBinLowEdge(ybin)}));
        }
      }

      auto dTs_names = Colib::file_get_names_of<TH3F>();
      for (auto const & name : dTs_names)
      {
        auto dTs = file->Get<TH3F>(name.c_str());
        std::string label_str = name;
        Colib::remove(label_str, "dT");
        int label = std::stoi(label_str);

        auto const & [fraction_bin, shift_bin] = m_minimaBin.at(label);
        auto dTName = "best_" + std::string(dTs->GetName()) + "_" + std::to_string(fraction_bin) + std::to_string(shift_bin);  
        m_best_dT.emplace(label, dTs->ProjectionZ(dTName.c_str(), fraction_bin+1, fraction_bin+1, shift_bin+1, shift_bin+1));
        m_best_dT.at(label)->SetDirectory(nullptr);
      }

      file->Close();
    }

    void writeRoot(std::string filename)
    {
      auto file = TFile::Open(filename.c_str(), "recreate");
      if (!file) Colib::throw_error(filename+" not created !! Is the path ok ?");
      for (auto & histos : m_histograms) histos.write();
      file->Close();
      print(filename, "written");
    }

    void write_dT(std::string filename)
    {
      if (m_minima.empty()) Colib::throw_error("Can't write dT because no minima have been found !!");
      std::ofstream file(filename);
      for (auto const & [label, min] : m_minima) 
      {
        file << label << " ";
        for (auto const & param : min) file << param << " ";
        file << "\n";
      }
      file.close();
      print(filename, "written");

      auto rootFilename = Colib::setExtension(filename, "root");
      auto resolutionsFilename = Colib::setExtension(filename, "resolutions");
      auto tfile = TFile::Open(rootFilename.c_str(), "recreate");
      if (!tfile) Colib::throw_error(rootFilename+" not created !! Is the path ok ?");
      std::ofstream rfile(resolutionsFilename);
      for (auto & [label, histo] : m_best_dT) 
      {
        if (!histo) {error(label, "not written"); continue;}
        auto const R = Colib::resolution(histo)/1000;
        rfile << label << " " << R << "\n"; 
        histo->Write();
      }
      tfile->Close();
      rfile.close();
      print(rootFilename, "written");
      print(resolutionsFilename, "written");
    }

    // void writeBest_dT(std::string filename)
    // {
    //   if (m_best_dT.empty()) Colib::throw_error("Can't write best dT because no  !!");
    //   auto file = TFile::Open(filename.c_str(), "recreate");
    //   if (!file) Colib::throw_error(filename+" not created !! Is the path ok ?");
    //   for (auto const & [label, histo] : m_best_dT)
    //   {
    //     histo -> Write();
    //   }
    // }
    
  private:
    CFDMinimisationParameters m_parameters;
    size_t m_nbDetectors = 0;
    std::vector<OptimizerHistograms> m_histograms;
    std::map<int, std::array<double, 2>> m_minima;
    std::map<int, std::array<int, 2>> m_minimaBin;
    std::map<int, TH1D*> m_best_dT;

    std::unique_ptr<TF1> gaus; 
    std::unique_ptr<TF1> gaus_and_bkgd;

  };
}