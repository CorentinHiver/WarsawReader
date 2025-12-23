#ifndef CAENDAT_HPP
#define CAENDAT_HPP

#include "utils.hpp"

namespace CaenDataReader1725
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
      auto temp = Colib::split(m_file.filename().shortName(), '_');
      try {
        m_fileNumber = std::stoi(temp.back());
        m_runNumber = std::stoi(temp[temp.size()-2]);
      }
      catch (std::invalid_argument const & error){
        Colib::throw_error(m_file.filename().string()+" : not a valid name (should be runName_boardVersion_runNumber_fileNumber.caendat)");
      }
    }
  
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
  
    auto const & filename() {return m_file.filename();}
    auto const & fileNumber() {return m_fileNumber;}
    auto const & runNumber() {return m_runNumber;}
    auto const & path() {return m_file.path();}
  
  private:
    Colib::File m_file;
    int m_fileNumber = -1;
    int m_runNumber = -1;
  };
};

#endif //CAENDAT_HPP