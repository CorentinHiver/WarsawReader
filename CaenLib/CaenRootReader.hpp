#ifndef CAENROOTREADER_HPP
#define CAENROOTREADER_HPP

#include "RootHit.hpp"
#include "CaenReaderBase.hpp"


namespace CaenDataReader1725
{
  /**
 * @brief Reads a .caendat file written by a 1725 board and provides an interface to write in a root files
 */
  class RootReader : public CaenReaderBase
  {
  public:
    RootReader(std::string const & filename, bool handle_traces = true) : 
      CaenReaderBase(filename)
    {
      m_rootHit.handleTraces(handle_traces);
    }

    /// @brief DEV - Resets the reader. The user needs to handle closing and/or writting the ROOT interface (TTrees, TFiles...)
    void reset()
    {
      m_nb_hits = 0;
    }

    void makePureVirtual(bool const & isVirtual = false) override {print(isVirtual);}; // To make this class real (printed to get rid of the warning)

    RootReader(std::string const & filename, TTree * tree, bool handle_traces = true) : 
      CaenReaderBase(filename), m_tree(tree)
    {
      m_rootHit.handleTraces(handle_traces);
      if (m_tree) m_rootHit.writeTo(m_tree);
      else Colib::throw_error("in RootReader::RootReader(std::string filename, TTree * tree) : tree is nullptr");
    }

    void writeTo(TTree * tree)
    {
      m_rootHit.writeTo(tree);
    }

    bool readHit() override
    {
      if (read_board_header)
      { // We hit a new board aggregate, reading the header
        m_board.clear();
        if (!m_board.readHeader(CaenReaderBase::p_datafile)) return false;
        read_board_header = false;
      }

      if (read_channel_header)
      { // We hit a new channel aggregate, reading the header
        m_channel.clear();
        m_board.readChannelHeader(m_channel, CaenReaderBase::p_datafile);
        read_channel_header = false;
      }

      if (m_channel.hasMoreEvents())
      { // There is at least one more caen event (=detector hit) to read
        ++m_nb_hits;
        m_caenEvent.clear();
        m_rootHit.readCaenEvent(CaenReaderBase::p_datafile, m_board, m_channel, m_caenEvent);
      }
      else
      { // At the end of each channel aggregate, different aggregate headers need proper handling: 
        read_channel_header = true; // Need to read the next channel aggregate header
        read_board_header = !m_board.hasMoreChannels(); // If this was last channel aggregate of the board aggregate, need to read the next board aggregate header
        return this -> readHit(); // Reads the next hit. Iteration : in principle, never called more than twice in a row
      }
      return true;
    }

    /**
     * @brief 
     * 
     * @param b 
     * @param resetConnection 
     */
    void handleTraces(bool const & b, bool resetConnection = true) 
    {
      m_rootHit.handleTraces(b);
      if (resetConnection && m_rootHit.writting && m_tree) writeTo(m_tree); // If the hit is already connected to a tree, then reset the connection to remove or add the trace 
    }

    auto & getHit() {return m_rootHit;}

  private:

    RootHit m_rootHit;
    TTree* m_tree = nullptr;
    
    BoardAggregate    m_board     ;
    ChannelAggregate  m_channel   ;
    CaenEvent         m_caenEvent ;

    bool read_board_header   = true;
    bool read_channel_header = true;
  };
};

/// @brief Alias to class CaenDataReader1725::RootReader : reads a .caendat file written by a 1725 board and provides an interface to write in a root files
using CaenRootReader1725 = CaenDataReader1725::RootReader;
#endif //CAENROOTREADER_HPP
