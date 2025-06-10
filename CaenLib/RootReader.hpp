#ifndef ROOTREADER_HPP
#define ROOTREADER_HPP


#include "RootHit.hpp"
#include "TFile.h"

class RootReader
{
public:
  RootReader() noexcept = default;
  RootReader(std::string const & filename)
  {
    m_file = TFile::Open(filename.c_str(), "READ");
    if (!m_file)  {error("RootReader::RootReader(std::string filename) : Can't read m_file" + filename); return;}
    m_tree = m_hit.readFrom(m_file, "HIL");
    if (!m_tree) return;
    m_size = m_tree->GetEntries();
  }

  RootReader(TFile * m_file)
  {
    if (!m_file) {error("RootReader::RootReader(TFile * m_file) : Can't read the given TFile"); return;}
    m_tree = m_hit.readFrom(m_file, "HIL");
    if (!m_tree) return;
    m_size = m_tree->GetEntries();
  }

  RootReader(TTree * tree)
  {
    m_tree = tree;
    if (!m_tree || m_tree->IsZombie()) {error("RootReader::RootReader(TTree * m_tree) : Can't read the given TTree"); return;}
    m_file = dynamic_cast<TFile*>(m_tree->GetDirectory());
    m_hit.readFrom(m_tree);
    m_size = m_tree->GetEntries();
  }

  bool readNext()
  {
    if (m_cursor < m_size)
    {
      m_tree->GetEntry(m_cursor++);
      return true;
    }
    return false;
  }

  void resetCursor() {m_cursor = 0;}
  
  auto const & getHit() const {return m_hit;}
  auto & getHit() {return m_hit;}
  auto getTree() {return m_tree;}
  auto getFile() {return m_file;}
  auto const & getCursor() const {return m_cursor;}

  void Scan(int nb_hits = -1)
  {
    // TODO
  }

private:

  RootCaenHit m_hit;
  TFile *m_file = nullptr;
  TTree *m_tree = nullptr;
  size_t m_cursor = 0;
  size_t m_size = 0;
};



#endif //ROOTREADER_HPP
