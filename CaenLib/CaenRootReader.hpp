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
    RootReader(std::string const & filename) : 
      CaenReaderBase(filename)
    {
    }

    void makePureVirtual(bool const & isVirtual = false) override {print(isVirtual);}; // To make this class real (printed to get rid of the warning)

    RootReader(std::string const & filename, TTree * tree) : 
      CaenReaderBase(filename)
    {
      if (tree) this->writeTo(tree);
    }

    void writeTo(TTree * tree)
    {
      m_rootHit.writeTo(tree);
    }

    bool readHit()
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
      { // There is at least one more event to read
        ++m_nb_hits;
        m_caenEvent.clear();
        m_rootHit.readCaenEvent(CaenReaderBase::p_datafile, m_board, m_channel, m_caenEvent);
      }
      else
      { // At the end of the channel aggregate, need to handle the different aggregate headers : 
        read_channel_header = true; // Need to read the next channel aggregate header
        read_board_header = !m_board.hasMoreChannels(); // If this was last channel aggregate of the board aggregate, need to read the next board aggregate header
        return this -> readHit(); // Reads the next hit. Iteration : should never be called more than twice in a row
      }
      return true;
    }

    void handleTraces(bool const & b) {m_rootHit.handle_traces = b;} 

    auto & getHit() {return m_rootHit;}
    auto const & nbHits() const {return m_nb_hits;}

  private:
    int m_nb_hits = 0;

    RootHit m_rootHit;
    
    BoardAggregate m_board;
    ChannelAggregate m_channel;
    CaenEvent m_caenEvent;

    bool read_board_header = true;
    bool read_channel_header = true;
  };
};

/// @brief Alias to class CaenDataReader1725::RootReader : reads a .caendat file written by a 1725 board and provides an interface to write in a root files
using CaenRootReader1725 = CaenDataReader1725::RootReader;
#endif //CAENROOTREADER_HPP
