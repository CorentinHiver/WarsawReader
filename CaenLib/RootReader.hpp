#ifndef ROOTREADER_HPP
#define ROOTREADER_HPP

#include "RootEvent.hpp"
#include "TFile.h"

/**
 * @brief Reads a .root files
 */

class RootReader
{
public:
  RootReader(int hits_nb = 0) : m_size(hits_nb) {};
  
  RootReader(TTree * tree, int hits_nb = 0) : RootReader(hits_nb)
  {
    m_tree = tree;
    connectTree(tree);
  }

  RootReader(TFile * file, int hits_nb = 0) : RootReader(hits_nb)
  {
    connectFile(file);
  }

  RootReader(std::string const & filename, int hits_nb = 0) : RootReader(hits_nb)
  {
    connectFile(filename);
  }

  TTree* connectTree(TTree* tree)
  {
    if (!tree) {error("in connectTree(TTree* tree) : tree is nullptr"); return nullptr;}
    m_tree = m_event.readFrom(tree);
    if (!m_grouped)
    {
      m_tree->SetBranchAddress("evtNb", &m_evtNb);
      m_tree->SetBranchAddress("mult", &m_evtMult);
    }
    m_size = m_tree->GetEntries();
    print(m_size);
    return tree;
  }

  TTree* connectFile(TFile * file)
  {
    if (!file) {error("in connectFile(TFile * file) : file is nullptr"); return nullptr;}
    m_file = file;
    auto const & listTrees = file_get_map_of<TTree>(m_file);
         if (Colib::key_found(listTrees, std::string("HIL")     )) {m_grouped = true ; return connectTree(listTrees.at("HIL"     ));}
    else if (Colib::key_found(listTrees, std::string("HILplain"))) {m_grouped = false; return connectTree(listTrees.at("HILplain"));}
    else return nullptr;
  }

  TTree* connectFile(std::string const & filename)
  {
    print(filename);
    m_file = TFile::Open(filename.c_str(), "READ");
    if (!m_file)  {error("RootReader::RootReader(std::string filename) : Can't read m_file" + filename); return nullptr;}
    return connectFile(m_file);
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
    if (m_grouped)
    {
      return readNextHit(); // Because we're in group mode, a TTree entry is not a single hit but already an event
    }
    else
    {
      if (m_finished) return false;
      if (readNextHit()) 
      {
        if (m_oldEvt < m_evtNb) m_oldEvt = m_evtNb;
        else m_event.push_back(m_hit);
      }
      return (m_finished = true);
    }
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
    if (m_grouped) while(this -> readNextEvent())
    {
      auto const & nbLines = Colib::getTerminalRows()-2;
      if (!stop && (m_cursor) % nbLines == 0) {
        println("Press enter for continuing, enter q for stopping (or Ctrl+C of course) ");
        std::getline(std::cin, user_input);
        if (user_input == "q") break;
      }
      print(m_event);
    }
    else while(this -> readNextHit())
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

  Caen1725RootHit m_hit;
  Caen1725RootEvent m_event;

  TFile *m_file = nullptr;
  TTree *m_tree = nullptr;
  size_t m_cursor = 0;
  size_t m_size = 0;
  bool m_grouped = true;

  size_t m_oldEvt = 0;
  size_t m_evtNb = 0;
  int m_evtMult = 0;
  bool m_finished = false;
};

#endif //ROOTREADER_HPP
