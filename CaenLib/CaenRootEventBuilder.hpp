#ifndef CAENROOTEVENTBUILDER_HPP
#define CAENROOTEVENTBUILDER_HPP

#include "RootHit.hpp"

namespace CaenDataReader
{
  class RootEventBuilder
  {
    /// @brief Holds the indices of coincident hits in the hit buffer
    using Event = std::vector<size_t>;

  public:
    RootEventBuilder()
    {
      m_hit_buffer.reserve(m_reserved_buffer_size);
    }

    RootEventBuilder(size_t const & reserved_buffer_size) : 
      m_reserved_buffer_size(reserved_buffer_size)
    {
      m_hit_buffer.reserve(m_reserved_buffer_size);
    }

    bool fill_buffer(RootHit const & hit)
    {
      if (m_hit_buffer.size() < m_reserved_buffer_size)
      {
        m_hit_buffer.emplace_back(hit);
        return true;
      }
      else return false;
    }

    bool fill_buffer(RootHit && hit)
    {
      if (m_hit_buffer.size() < m_reserved_buffer_size)
      {
        m_hit_buffer.emplace_back(std::move(hit));
        return true;
      }
      else return false;
    }

    void align()
    {
      Colib::linspace(m_ordered_index, m_hit_buffer.size());
      std::sort(m_ordered_index.begin(), m_ordered_index.end(), [this](size_t const & i, size_t const & j)
      {
        return m_hit_buffer[j] > m_hit_buffer[i];
      });
      m_aligned = true;
    }

    void fast_event_building(Long64_t const & time_window) noexcept
    {
      
      // 1. Initialize the event buffer
          
      if (!m_aligned) this -> align();
      Event event;
      event.emplace_back(m_ordered_index[0]); 

      // 2. Loop through the hits buffer

      for (size_t loop_i = 1; loop_i < m_hit_buffer.size(); ++loop_i)
      {
        auto const & hit_i      =  m_ordered_index[loop_i       ];
        auto const & hit        =  m_hit_buffer   [hit_i        ];
        auto const & first_hit  =  m_hit_buffer   [event.front()];

        // 3. Add new hits until one is out of time window with the first hit of the event
        // moment referred to as "closing the event"

        if (Long64_t(hit.timestamp - first_hit.timestamp) < time_window) 
        {
          event.emplace_back(hit_i);
          continue;
        }

        // -- Piece of code only executed when the event is closed -- //
        
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

    auto const & operator[](size_t const & i) const {return m_hit_buffer[i];}
    auto & operator[](size_t const & i) {return m_hit_buffer[i];}

    auto const & getHitBuffer() const {return m_hit_buffer;}

    auto begin()       {return m_event_buffer.begin();}
    auto end  ()       {return m_event_buffer.end  ();}
    auto begin() const {return m_event_buffer.begin();}
    auto end  () const {return m_event_buffer.end  ();}

  private:

    size_t m_reserved_buffer_size = 5000ul;

    std::vector<size_t> m_ordered_index;
    std::vector<RootHit> m_hit_buffer;
    std::vector<Event> m_event_buffer;
    bool m_aligned = false;
  };
};

using CaenRootEventBuilder = CaenDataReader::RootEventBuilder;

#endif //CAENROOTEVENTBUILDER_HPP
