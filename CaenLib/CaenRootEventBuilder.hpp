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
      hit_buffer.reserve(m_reserved_buffer_size);
    }

    RootEventBuilder(size_t const & reserved_buffer_size) : 
      m_reserved_buffer_size(reserved_buffer_size)
    {
      hit_buffer.reserve(m_reserved_buffer_size);
    }

    bool fill_buffer(RootHit const & hit)
    {
      if (hit_buffer.size() < m_reserved_buffer_size)
      {
        hit_buffer.emplace_back(hit);
        return true;
      }
      else return false;
    }

    void fast_event_building(Long64_t const & time_window)
    {
      std::vector<size_t> ordered_index;
      bubble_sort(hit_buffer, ordered_index);
      Event event;
      event.emplace_back(ordered_index[0]);

      for (size_t loop_i = 1; loop_i < hit_buffer.size(); ++loop_i)
      {
        auto const & hit_i = ordered_index[loop_i];
        auto const & hit = hit_buffer[hit_i];
        auto const & first_hit = hit_buffer[event.front()];

        if (Long64_t(hit.timestamp - first_hit.timestamp) < time_window) 
        {
          event.emplace_back(hit_i);
          continue;
        }
        event_buffer.emplace_back(event);

        ////////////////////
        // Event Analysis //
        ////////////////////

        // Prepare next event : 
        event.clear();
        event.emplace_back(hit_i);
      }
    }

    void clear(){
      hit_buffer  .clear();
      event_buffer.clear();
    }

    auto const & operator[](size_t const & i) const {return hit_buffer[i];}

    auto begin()       {return event_buffer.begin();}
    auto end  ()       {return event_buffer.end  ();}
    auto begin() const {return event_buffer.begin();}
    auto end  () const {return event_buffer.end  ();}

  private:
    std::vector<RootHit> hit_buffer;
    std::vector<Event> event_buffer;
    size_t m_reserved_buffer_size = 5000ul;
  };
};

using CaenRootEventBuilder = CaenDataReader::RootEventBuilder;

#endif //CAENROOTEVENTBUILDER_HPP
