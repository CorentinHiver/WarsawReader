#pragma once
#include "RootEvent.hpp"
#include "TFile.h"

namespace Caen1725
{
  class RootReader
  {
  public:
    RootReader(size_t hits_nb = Colib::max<int>()) : m_size(hits_nb) {};
    
    RootReader(TTree * tree, size_t hits_nb = Colib::max<int>()) : RootReader(hits_nb) 
      {connectTree(tree);}  
    RootReader(TFile * file, size_t hits_nb = Colib::max<int>()) : RootReader(hits_nb) 
      {connectFile(file);}  
    RootReader(std::string const & filename, size_t hits_nb = Colib::max<int>()) : RootReader(hits_nb) 
      {connectFile(filename);}

    ~RootReader() {if (m_file && !m_file->IsZombie()) m_file->Close();}
  
    TTree* connectTree(TTree* tree)
    {
      if (!tree) {error("in connectTree(TTree* tree) : tree is nullptr"); return nullptr;}
      m_tree = m_event.readFrom(tree);
      if (m_plain)
      {
        m_hit.readFrom(tree);
        m_tree->SetBranchAddress("eventID", &m_eventID);
        m_tree->SetBranchAddress("mult"   , &m_evtMult);
      }
      m_size = std::min(m_size, static_cast<size_t>(m_tree->GetEntries()));
      print(m_size);
      return tree;
    }
  
    TTree* connectFile(TFile * file)
    {
      if (!file) {error("in connectFile(TFile * file) : file is nullptr"); return nullptr;}
      m_file = file;
      auto const & listTrees = file_get_map_of<TTree>(m_file);
           if (Colib::key_found(listTrees, std::string("HIL")     )) 
      {
        m_plain = false ; 
        return connectTree(listTrees.at("HIL"));
      }
      else if (Colib::key_found(listTrees, std::string("HILplain"))) 
      {
        m_plain = true; 
        return connectTree(listTrees.at("HILplain"));
      }
      else return nullptr;
    }
  
    TTree* connectFile(std::string const & filename)
    {
      m_file = TFile::Open(filename.c_str(), "READ");
      if (!m_file)  {error("RootReader::RootReader(std::string filename) : Can't read m_file" + filename); return nullptr;}
      return connectFile(m_file);
    }
        
    /// @brief Advanced users only. Read next event in grouped mode, or next hit in plain mode
    /// @return false if last entry
    bool readNextEntry()
    {
      if (m_cursor < m_size)
      {
        m_tree->GetEntry(m_cursor++);
        return true;
      }
      return false; 
    }

    /// @brief Read next event
    /// @return false if last event
    bool readNextEvent()
    {
      // In plain mode, a TTree entry is a single hit. This loop is therefor used to reconstruct the full event
      if (m_plain)
      {
        m_event.clear();
        bool continuing = false;
        while((continuing = readNextEntry()) && m_event.size() < size_cast(m_evtMult))
          m_event.push_back(m_hit);
        return continuing;
      }
      // In event mode, a TTree entry is already an event
      else return readNextEntry(); 
    }
  
    // Setters :
  
    /// @brief Starts reading the tree from the beginning
    void resetCursor() {m_cursor = 0;}
    
    // Getters :
    
    /// @brief Get a read-only access to the current hit.
    auto const & getHit()   const {return m_hit  ;}
    /// @brief Get an access to the current hit.
    auto       & getHit()         {return m_hit  ;}
    /// @brief Get a read-only access to the current event.
    auto const & getEvent() const {return m_event;}
    /// @brief Get an access to the current event.
    auto       & getEvent()       {return m_event;}
  
    /// @brief Advanced users only. Get a pointer to the TTree.
    auto getTree() {return m_tree;}
    /// @brief Advanced users only. Get a pointer to the TFile.
    auto getFile() {return m_file;}
  
    /// @brief Get a read-only access to the cursor.
    auto const & getCursor() const {return m_cursor;}
  
    /// @brief 
    void Scan(bool scanHits = false)
    {
      std::string user_input;
      int nbLinesWritten = 0;
      if (m_plain) while((scanHits) ? this->readNextEntry() : this -> readNextEvent())
      {
        if (Colib::Terminal::getRows()-2 < ++nbLinesWritten) 
        {
          nbLinesWritten = 0;
          println("Press enter for continuing, enter q for stopping (or Ctrl+C of course) ");
          std::getline(std::cin, user_input);
          if (user_input == "q") break;
        }
        if (scanHits) print("EventID", m_eventID, "eventMult", m_evtMult, m_hit);
        else print("EventID", m_eventID, "eventMult", m_evtMult, m_event);
      }
      else while(this -> readNextEvent())
      {
        if (Colib::Terminal::getRows()-2 < ++nbLinesWritten) 
        {
          nbLinesWritten = 0;
          println("Press enter for continuing, enter q for stopping (or Ctrl+C of course) ");
          std::getline(std::cin, user_input);
          if (user_input == "q") break;
        }
        if (scanHits) print(m_event);
        else for (size_t hit_i = 0; hit_i<m_event.size(); ++hit_i) print(m_hit);
      }
    }
  
  private:
  
    RootHit m_hit;
    RootEvent m_event;
  
    TFile *m_file = nullptr;
    TTree *m_tree = nullptr;
    size_t m_cursor = 0;
    size_t m_size = 0;
    bool m_plain = false;
  
    Caen1725::EventID   m_eventID = 0;
    Caen1725::EventMult m_evtMult = 0;
  };
}