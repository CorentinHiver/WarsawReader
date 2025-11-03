#ifndef ROOTREADER_HPP
#define ROOTREADER_HPP

#include "RootHit.hpp"
#include "TFile.h"

/**
 * @brief Reads a .root files
 */

class RootReader
{
public:
  RootReader(int hits_nb = 0) : m_size(hits_nb) {};
  RootReader(std::string const & filename, int hits_nb = 0) : RootReader(hits_nb)
  {
    m_file = TFile::Open(filename.c_str(), "READ");
    if (!m_file)  {error("RootReader::RootReader(std::string filename) : Can't read m_file" + filename); return;}
    m_tree = m_hit.readFrom(m_file, "HIL");
    m_tree->SetBranchAddress("evtNb", &m_evtNb);
    if (!m_tree) return;
    if (m_size == 0) m_size = m_tree->GetEntries();
  }

  RootReader(TFile * m_file, int hits_nb = 0) : RootReader(hits_nb)
  {
    if (!m_file) {error("RootReader::RootReader(TFile * m_file) : Can't read the given TFile"); return;}
    m_tree = m_hit.readFrom(m_file, "HIL");
    if (!m_tree) return;
    m_size = m_tree->GetEntries();
    m_tree->SetBranchAddress("evtNb", &m_evtNb);
    if (m_size == 0) m_size = m_tree->GetEntries();
  }

  RootReader(TTree * tree, int hits_nb = 0) : RootReader(hits_nb)
  {
    m_tree = tree;
    if (!m_tree || m_tree->IsZombie()) {error("RootReader::RootReader(TTree * m_tree) : Can't read the given TTree"); return;}
    m_file = dynamic_cast<TFile*>(m_tree->GetDirectory());
    m_hit.readFrom(m_tree);
    m_tree->SetBranchAddress("evtNb", &m_evtNb);
    if (m_size == 0) m_size = m_tree->GetEntries();
  }

  bool readNextHit()
  {
    if (m_cursor < m_size)
    {
      m_tree->GetEntry(m_cursor++);
      return true;
    }
    return false;
  }

  bool readNextEvent()
  {
    if (m_finished) return false;
    if (readNextHit()) 
    {
      if (m_oldEvt < m_evtNb) m_oldEvt = m_evtNb;
      else m_event.push_back(m_hit);
    }
    return (m_finished = true);
  }

  // Setters :

  void resetCursor() {m_cursor = 0;}
  
  // Getters :
  
  auto const & getHit()   const {return m_hit  ;}
  auto       & getHit()         {return m_hit  ;}
  auto const & getEvent() const {return m_event;}
  auto       & getEvent()       {return m_event;}

  auto getTree() {return m_tree;}
  auto getFile() {return m_file;}

  auto const & getCursor() const {return m_cursor;}

  void Scan(bool stop = 0)
  {
    std::string user_input;
    while(this -> readNextHit())
    {
      auto const & nbLines = Colib::getTerminalRows()-2;
      if (!stop && (m_cursor) % nbLines == 0) {
        println("Press enter for continuing, enter q for stopping (or Ctrl+C of course) ");
        std::getline(std::cin, user_input);
        if (user_input == "q") break;
      }
      print(m_evtNb, m_hit);
    }
  }

private:

  RootCaenHit m_hit;
  RootCaenEvent m_event;

  TFile *m_file = nullptr;
  TTree *m_tree = nullptr;
  size_t m_cursor = 0;
  size_t m_size = 0;

  int m_oldEvt = 0;
  int m_evtNb = 0;
  bool m_finished = false;
};

#endif //ROOTREADER_HPP
