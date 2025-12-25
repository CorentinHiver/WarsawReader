#ifndef CAENROOTREADER_HPP
#define CAENROOTREADER_HPP

#include "RootHit.hpp"
#include "CaenReaderBase.hpp"
#include "TROOT.h"

namespace CaenDataReader1725
{
  /**
 * @brief Reads a .caendat file written by a 1725 board and provides a hit interface to write to a root file
 */
  class RootInterface : public CaenReaderBase
  {
  public:
    /// @brief ROOT I/O handled by the user.
    RootInterface(std::string const & caenFilename, bool handle_traces = true) : 
      CaenReaderBase(caenFilename)
    {
      handleTraces(handle_traces);
    }

    /// @brief DEV - Resets the reader. The user needs to handle closing and/or writting the ROOT interface (TTrees, TFiles...)
    void reset()
    {
      m_nb_hits = 0;
    }

    /// @brief Internal ROOT tree. ROOT file I/O handled by the user.
    RootInterface(std::string const & caenFilename, TTree * tree, bool handle_traces = true) : 
      CaenReaderBase(caenFilename), m_tree(tree)
    {
      handleTraces(handle_traces);
      if (m_tree) m_rootHit.writeTo(m_tree);
      else Colib::throw_error("in RootInterface::RootInterface(std::string caenFilename, TTree * tree) : tree is nullptr");
    }

    /// @brief Internal ROOT tree and file. Internal I/O management, call RootInterface::Write() to write the TTree in the TFile and clean everything.
    RootInterface(std::string const & caenFilename, std::string const & rootFilename, TTree * tree, bool handle_traces = true) : 
      CaenReaderBase(caenFilename), m_tree(tree)
    {
      m_file = TFile::Open(rootFilename.c_str(), "RECREATE");
      handleTraces(handle_traces);
      if (m_tree) m_rootHit.writeTo(m_tree);
      else Colib::throw_error("in RootInterface::RootInterface(std::string caenFilename, TTree * tree) : tree is nullptr");
    }

    void Write()
    {
      if (!m_file) Colib::throw_error("in RootInterface::Write() : No output file !!");
      if (!m_tree) Colib::throw_error("in RootInterface::Write() : No output tree !!");
      m_file -> cd();
      m_tree -> Write();
      m_file -> Close();
      delete m_file;
    }

    void setTreeInMemory(bool b) {if (b) gROOT->cd();}

    ~RootInterface() 
    {
      if (m_file && m_file->IsZombie() && m_file -> IsOpen()) m_file -> Close();
    }

    void setBoardReadTrace(std::vector<bool> boards)
    {
      m_boardsReadTraces = boards;
    }

    void writeTo(TTree * tree)
    {
      m_rootHit.writeTo(tree);
    }

    bool readHit() override
    {
      std::ifstream& data = CaenReaderBase::p_datafile; // Simple aliasing

      if (read_board_header)
      { // New board aggregate, reading its header
        m_board.clear();
        if (!m_board.readHeader(data)) return false;
        read_board_header = false;
      }

      if (read_channel_header)
      { // New channel aggregate, reading its header
        m_channel.clear();
        m_board.readChannelHeader(m_channel, data); // Crucial to use BoardAggregate::readChannelHeader and not ChannelAggregate::readHeader because it loads the mask and increments the size_read of the board (crucial for m_board.hasMoreChannels())
        read_channel_header = false;
      }

      if (m_channel.hasMoreEvents())
      { // There is at least one more caen event (=detector hit) to read
        ++m_nb_hits;
        m_caenEvent.clear();
        m_rootHit.clear();
        m_rootHit.readCaenEvent(data, m_board, m_channel, m_caenEvent, m_boardsReadTraces);
      }
      else
      { // At the end of each channel aggregate, different aggregate headers need proper handling: 
        read_channel_header = true; // Need to read the next channel aggregate header
        read_board_header = !m_board.hasMoreChannels(); // If this was last channel aggregate of the board aggregate, need to read the next board aggregate header
        return this -> readHit(); // Reads the next hit. Iteration : in principle, never called more than twice in a row
      }
      return true;
    }

    virtual void handleTraces(bool b) {m_rootHit.handleTraces(b);}

    auto & getHit() {return m_rootHit;}

  private:
    void makePureVirtual(bool isVirtual = false) override {print(isVirtual);}; // To make this class real (printed to get rid of the warning)
    
    RootHit m_rootHit;
    TTree*  m_tree = nullptr;
    TFile*  m_file = nullptr;
    
    BoardAggregate    m_board     ;
    ChannelAggregate  m_channel   ;
    CaenEvent         m_caenEvent ;

    bool read_board_header   = true;
    bool read_channel_header = true;

    std::vector<bool> m_boardsReadTraces;
  };
};

/// @brief Alias to class CaenDataReader1725::RootInterface : reads a .caendat file written by a 1725 board and provides an interface to write in a root files
using Caen1725RootInterface = CaenDataReader1725::RootInterface;
#endif //CAENROOTREADER_HPP
