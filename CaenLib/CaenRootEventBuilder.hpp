#ifndef CAENROOTEVENTBUILDER_HPP
#define CAENROOTEVENTBUILDER_HPP

#include "RootHit.hpp"

namespace Caen1725
{
  class EventBuilder
  {
    public:
    
    /// @brief Holds the indices of coincident hits in the hit buffer
    using EventId = std::vector<size_t>;

    EventBuilder(size_t reserved_buffer_size = (50000ul)) noexcept
    {
      m_hit_buffer.reserve(reserved_buffer_size);
    }

    /// @brief Forbidden to copy a RootHit
    bool fill_buffer(RootHit const & hit) noexcept = delete;
    // {
    //   if (m_hit_buffer.size() < m_hit_buffer.capacity())
    //   {
    //     m_hit_buffer.emplace_back(hit);
    //     return true;
    //   }
    //   else return false;
    // }

    bool fill_buffer(RootHit && hit) noexcept
    {
      if (m_hit_buffer.size() < m_hit_buffer.capacity())
      {
        m_hit_buffer.emplace_back(std::move(hit));
        return true;
      }
      else return false;
    }

    void align()
    {
      Colib::linspace(m_ordered_index, m_hit_buffer.size());
      if (m_buildOnTimestamp) std::sort(m_ordered_index.begin(), m_ordered_index.end(), [this](size_t i, size_t j){
        return m_hit_buffer[j].timestamp > m_hit_buffer[i].timestamp;
      });
      else                    std::sort(m_ordered_index.begin(), m_ordered_index.end(), [this](size_t i, size_t j){
        return m_hit_buffer[j].time      > m_hit_buffer[i].time     ;
      });
      m_aligned = true;
    }

    void fast_event_building(Long64_t time_window) noexcept
    {
      // 1. Initialize the event buffer
      if (!m_aligned) this -> align();
      EventId event;
      event.emplace_back(m_ordered_index[0]); 
      // 2. Loop through the hits buffer
      for (size_t loop_i = 1; loop_i < m_hit_buffer.size(); ++loop_i)
      {
        auto const & hit_i      =  m_ordered_index[loop_i       ];
        auto const & hit        =  m_hit_buffer   [hit_i        ];
        auto const & first_hit  =  m_hit_buffer   [event.front()];
        // 3. Add new hits until one is out of time window with the first hit of the event
        // moment referred to as "closing the event"
        if ((m_buildOnTimestamp) ? (static_cast<Long64_t>(hit.timestamp - first_hit.timestamp) < time_window)
                                 : (static_cast<Long64_t>(hit.time      - first_hit.time     ) < time_window)) 
        {
          event.emplace_back(hit_i);  
          continue;
        }
        // -- Piece of code only executed when the event is full -- //
        // 4. Fill the event buffer
        m_event_buffer.emplace_back(event);
        // 5. Prepare next event : 
        event.clear();
        event.emplace_back(hit_i); // Save the current hit
      }
    }

    void clear() noexcept
    {
      m_hit_buffer   .clear();
      m_event_buffer .clear();
      m_ordered_index.clear();
      m_aligned = false;
    }

    auto const & operator[](size_t i) const {return m_hit_buffer[i];}
    auto & operator[](size_t i) {return m_hit_buffer[i];}

    auto const & getHitBuffer() const {return m_hit_buffer;}

    auto begin()       {return m_event_buffer.begin();}
    auto end  ()       {return m_event_buffer.end  ();}
    auto begin() const {return m_event_buffer.begin();}
    auto end  () const {return m_event_buffer.end  ();}

    void buildOnTimestamp(bool b) {m_buildOnTimestamp = b;}

  private:

    std::vector<size_t>  m_ordered_index;
    std::vector<RootHit> m_hit_buffer   ;
    std::vector<EventId> m_event_buffer ;
    bool m_aligned = false;

    // Options
    bool m_buildOnTimestamp = false;
  };
};

using Caen1725EventBuilder = Caen1725::EventBuilder;

#endif //CAENROOTEVENTBUILDER_HPP
