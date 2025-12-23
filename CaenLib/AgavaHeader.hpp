#ifndef AGAVAHEADER_HPP
#define AGAVAHEADER_HPP

#include "utils.hpp"

namespace CaenDataReader1725
{
  class AgavaHeader
  {
  public:
    AgavaHeader(std::ifstream & datafile){
      std::cout << "Reading Agava header" << std::endl;
      read_data(datafile, &tmp_word); // Word 1
      m_nb_words = size_cast(getBitField(tmp_word, 27, 0));
      read_data(datafile, &tmp_word); // Word 2
      // Some parameters
      read_data(datafile, &tmp_word); // Word 3
      m_nb_cards = size_cast(tmp_word);
      std::cout << m_nb_cards << " cards" << std::endl;
      read_data(datafile, &tmp_word); // Word 4
      // Some other parameters
      read_data(datafile, &tmp_word); // Word 5
      // Some other parameters

      // Looping through the m_nb_cards cards parameters of the agva header:
      for (size_t word_i = 5; word_i<m_nb_words; ++word_i){
          read_data(datafile, &tmp_word);
          // Some other parameters
      }
    }

  private:
    uint32_t tmp_word;
    size_t m_nb_words = 0;
    size_t m_nb_cards = 0;
  };
}


#endif //AGAVAHEADER_HPP
