#ifndef CAENRAWREADER_HPP
#define CAENRAWREADER_HPP

#include "CaenReaderBase.hpp"

namespace Caen1725
{
  class RawReader : public CaenReaderBase
  {
  public:
    RawReader(std::string const & filename) : CaenReaderBase(filename)
    {}

    struct ErrorEof
    {
      public: ErrorEof(){}
    } static errorEof;

    /**
     * @brief Main function, reads an entire board aggregate
     * IS NOT compatible with skipAll() == true
     * 
     * @return true while the end of the file is not reached
     */
    bool readBoardAggregate()
    {
      if (!readBoardHeader()) return false;

        while(boardNewChannel()) 
          while ((sSkipData) ? skipEvent() : readEvent())
            if (CaenReaderBase::p_datafile.eof()) throw errorEof;

        return !CaenReaderBase::p_datafile.eof();
    }

    /**
     * @brief Does exactly the same as RawReader::readBoardAggregate, 
     * but explicits each step instead of calling functions.
     * IS compatible with skipAll() == true
     */
    bool readBoardAggregatePlain()
    {
      //////////////////////////
      // READ BOARD AGGREGATE //
      //////////////////////////

      // 1. Size and Check bin :
      read_data(CaenReaderBase::p_datafile, &tmp_u32); // Read the whole word
      debug("size, check_bin:", std::bitset<32>(tmp_u32));

      if (CaenReaderBase::p_datafile.eof()) return false;

      if (getBitField (tmp_u32, 31, 28) != m_board.check_bin_ref) // Check whether the check bin is correct
      {
        error("Check bin missed : ", std::bitset<32>(tmp_u32)); 
        return false;
      }

      m_board.clear();

      m_board.read_size += sizeof(tmp_u32); // Accounting for the previous call to read_data

      m_board.size      = getBitField (tmp_u32, 27    ) * sizeof(tmp_u32); // To get the number of octets
      m_board.check_bin = getBitField (tmp_u32, 31, 28);

      if (sSkipData) skip(CaenReaderBase::p_datafile, m_board.header_size - sizeof(tmp_u32));
      else 
      {
        // 2. Miscellaneous :
        read_data(CaenReaderBase::p_datafile, &tmp_u32, m_board.read_size); // Read the whole word
        debug("word n°2", std::bitset<32>(tmp_u32));
        m_board.DUAL_CHANNEL_MASK = getBitField (tmp_u32, 7     );
        m_board.PATTERN           = getBitField (tmp_u32, 22,  8);
        m_board.BF                = getBit      (tmp_u32, 26    );
        m_board.BOARD_ID          = getBitField (tmp_u32, 31, 27);
  
        // Handle the channel mask that allows to know the channel ID
        m_board.maskHelper = std::bitset<8>(m_board.DUAL_CHANNEL_MASK);
  
        // 3. Counter :
        read_data(CaenReaderBase::p_datafile, &m_board.COUNTER, m_board.read_size);
        debug("COUNTER", std::bitset<32>(m_board.COUNTER));
        m_board.COUNTER &= mask(22);
  
        // 4. Timestamp : 
        read_data(CaenReaderBase::p_datafile, &m_board.TIME_TAG, m_board.read_size);
        debug("TIME_TAG", std::bitset<32>(m_board.TIME_TAG));
      }    

      while (m_board.hasMoreChannels())
      {
        ////////////////////////////
        // READ CHANNEL AGGREGATE //
        ////////////////////////////

        m_board.channels.emplace_back();

        auto & channel = m_board.channels.back();

        channel.handle_traces = m_board.handle_traces;
        
        channel.ID = m_board.maskHelper.getID();

        // 1. Size and FI
        read_data(CaenReaderBase::p_datafile, &tmp_u32, m_board.read_size);
        debug("size, FI:", std::bitset<32>(tmp_u32));
        channel.size = getBitField(tmp_u32, 30) * sizeof(uint32_t);
        channel.FI   = getBit     (tmp_u32, 31);

        if (channel.FI)
        {
          // 2. Format
          read_data(CaenReaderBase::p_datafile, &tmp_u32, m_board.read_size);
          debug("Format:", std::bitset<32>(tmp_u32));

          channel.NUM_SAMPLES = getBitField (tmp_u32, 15) * 8; // The data stored is NUM_SAMPLES/8 (see doc)
          channel.DP          = getBitField (tmp_u32, 19, 16);
          channel.AP2         = getBitField (tmp_u32, 21, 20);
          channel.AP1         = getBitField (tmp_u32, 23, 22);
          channel.EX          = getBitField (tmp_u32, 26, 24);

          channel.ES = getBit (tmp_u32, 27);
          channel.E2 = getBit (tmp_u32, 28);
          channel.ET = getBit (tmp_u32, 29);
          channel.EE = getBit (tmp_u32, 30);
          channel.DT = getBit (tmp_u32, 31);
        }

        channel.read_size += 2 * sizeof(tmp_u32); // Take into account the size of the header

        while (channel.hasMoreEvents())
        {
          auto const pos_begin_event = CaenReaderBase::p_datafile.tellg();

          if (sSkipData) skip(CaenReaderBase::p_datafile, sizeof(tmp_u32) + channel.NUM_SAMPLES * sizeof(channel.NUM_SAMPLES) + 2*sizeof(tmp_u32));

          /////////////////////
          // READ CAEN EVENT //
          /////////////////////

          // CaenEvent event;
          channel.events.emplace_back();
          auto & event = channel.events.back();

          read_data(CaenReaderBase::p_datafile, &tmp_u32);
          debug("TRIGGER_TIME_TAG, CH:", std::bitset<32>(tmp_u32));
          
          event.TRIGGER_TIME_TAG = getBitField(tmp_u32, 30);
          event.CH               = getBit     (tmp_u32, 31);
          
          // Looping through the NUM_SAMPLES samples of the waveform of each event :
          if (channel.handle_traces) {
            debug("Trace with :", channel.NUM_SAMPLES, "samples");
            event.samples.resize(channel.NUM_SAMPLES);
            for (auto & sample : event.samples)
            {
              read_data(CaenReaderBase::p_datafile, &tmp_u16);
              sample.sample  = getBitField(tmp_u16, 13);
              sample.DP1     = getBit     (tmp_u16, 14);
              sample.T       = getBit     (tmp_u16, 15);
            }
          }
          else {
            auto const & size_to_skip = channel.NUM_SAMPLES * sizeof(channel.NUM_SAMPLES);
            skip(CaenReaderBase::p_datafile, size_to_skip);
          }
      
          // if (channel.E2)  // TODO: do we need to read EXTRAS2 if it is disabled, i.e. if channel.E2=false ?
          read_data(CaenReaderBase::p_datafile, &event.EXTRAS2);
          debug("EXTRAS2", std::bitset<32>(event.EXTRAS2));
      
          event.EX = channel.EX;
          event.extra = Extra2(event.EXTRAS2, event.TRIGGER_TIME_TAG, event.EX);
          
          read_data(CaenReaderBase::p_datafile, &tmp_u32);
      
          debug("ENERGY, PU, EXTRAS:", std::bitset<32>(tmp_u32));
      
          event.ENERGY = getBitField (tmp_u32, 14    );
          event.PU     = getBit      (tmp_u32, 15    );
          event.EXTRAS = getBitField (tmp_u32, 25, 16);
          
          auto const & event_length = int_cast(CaenReaderBase::p_datafile.tellg() - pos_begin_event);
          channel.read_size  += event_length;
          m_board.read_size  += event_length;

          ++m_nb_hits;
        }
      }
      return true;
    }

    auto const & getBoard() const {return m_board;}
    auto       & getBoard()       {return m_board;}
    
    void handleTraces(bool b = true) {m_board.handle_traces = b;} 
    static void skipData(bool b = true) {sSkipData = b;}
    
    bool readBoardHeader() {return m_board.readHeader(CaenReaderBase::p_datafile);}
    bool boardNewChannel() {return m_board.newChannel(CaenReaderBase::p_datafile);}
    bool readEvent()       {return m_board.readEvent (CaenReaderBase::p_datafile);}    
    bool skipEvent()       {return m_board.skipEvent (CaenReaderBase::p_datafile);}    

  private:
  
    void makePureVirtual(bool isVirtual = false) override {print(isVirtual);}; // To make this class real (print to get rid of the warning)
    static inline bool sSkipData = false;
    BoardAggregate m_board;
  };
};

using CaenRawReader1725 = Caen1725::RawReader;

#endif //CAENRAWREADER_HPP