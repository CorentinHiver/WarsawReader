#ifndef CAENREADER_HPP
#define CAENREADER_HPP

#include "CaenDat.hpp"
#include "AgavaHeader.hpp"
#include "BoardAggregate.hpp"

namespace CaenDataReader1725
{
  class CaenReaderBase
  {
  public:
    CaenReaderBase(std::string const & filename)
    {
      this->open(filename);
    }

    virtual void makePureVirtual(bool const & isVirtual = true) = 0; // To make this class pure virtual

    virtual ~CaenReaderBase() {if (p_datafile.is_open()) p_datafile.close();}

    virtual void open(std::string const & filename)
    {
      m_filename = filename  ;
      caenFile   = m_filename;
      p_datafile.open(m_filename, std::ios::binary | std::ios::in);
      if (!p_datafile.is_open())
      {
        error(std::strerror(errno));
        Colib::throw_error("Could not open file '"+ m_filename+"'");
        return;
      }

      print("Reading file", m_filename);

      // Print agava header informations for the first file of each run :
      if (caenFile.fileNumber() == 0) 
      {
        AgavaHeader ah(p_datafile);
      }
    }

    virtual bool eof() {return p_datafile.eof();}

  protected:    
    std::ifstream p_datafile;
    std::string   m_filename;
    CaenDat caenFile;
  };
};
#endif //CAENREADER_HPP
