#ifndef CAENREADER_HPP
#define CAENREADER_HPP

#include "CaenDat.hpp"
#include "AgavaHeader.hpp"
#include "BoardAggregate.hpp"

namespace Caen1725
{
  class CaenReaderBase
  {
  public:
    CaenReaderBase(std::string const & filename)
    {
      open(filename);
    }

    virtual ~CaenReaderBase() {if (p_datafile.is_open()) p_datafile.close();}

    virtual void open(std::string const & filename)
    {
      caenFile = CaenDat(m_filename = filename);
      p_datafile.open(m_filename, std::ios::binary | std::ios::in);
      if (!p_datafile.is_open())
      {
        error(std::strerror(errno));
        Colib::throw_error("Could not open file '"+ m_filename+"'");
        return;
      }

      print("Reading file", m_filename);

      // Print agava header informations for the first file of each run :
      if (s_hasAgava && caenFile.fileNumber() == 0) 
      {
        AgavaHeader ah(p_datafile);
      }
    }

    virtual bool eof() {return p_datafile.eof();}

    bool readHit() {++m_nb_hits; return false;}
    auto const & nbHits() const {return m_nb_hits;}

    static void setAgavaHeader(bool b = true) noexcept {s_hasAgava = b;}

  protected:    
    virtual void makePureVirtual(bool isVirtual = true) = 0; // To make this class pure virtual
    
    size_t        m_nb_hits = 0;
    std::ifstream p_datafile;
    std::string   m_filename;
    CaenDat       caenFile;

    inline static bool s_hasAgava = true; 
  };
};
#endif //CAENREADER_HPP
