#ifndef CAENDAT_HPP
#define CAENDAT_HPP

#include "utils.hpp"

namespace Caen1725
{
  class CaenDat
  {
  public:
    CaenDat() {}
    CaenDat(std::string const & filename) : m_file(filename)
    {
      setup();
    }
    void setup()
    {
      if (Colib::extension(m_file) != "caendat") error(m_file, "not a .caendat file !!");
      auto temp = Colib::split(Colib::getShortname(m_file), '_');
      try {
        m_fileNumber = std::stoi(temp.back());
        m_runNumber = std::stoi(temp[temp.size()-2]);
      }
      catch (std::invalid_argument const & error){
        Colib::throw_error(Colib::removePath(m_file)+" : not a valid name (should be runName_boardVersion_runNumber_fileNumber.caendat)");
      }
      m_ok = true;
    }

    constexpr operator bool() const {return m_ok;}
  
    void read(std::string const & filename)
    {
      m_file = filename;
      setup();
    }
  
    CaenDat& operator=(std::string const & filename)
    {
      this->read(filename);
      return *this;
    }
  
    auto filename() const {return Colib::removePath(m_file);}
    auto const & fileNumber() const {return m_fileNumber;}
    auto const & runNumber() const {return m_runNumber;}
    auto path() const {return Colib::getPath(m_file);}
  
  private:
    std::string m_file;
    int m_fileNumber = -1;
    int m_runNumber = -1;
    bool m_ok = false;
  };
};

#endif //CAENDAT_HPP