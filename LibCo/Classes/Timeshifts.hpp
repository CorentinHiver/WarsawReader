#ifndef TIMESHIFTS_HPP
#define TIMESHIFTS_HPP

// #include "../Classes/Hit.hpp"
#include "../libCo.hpp"

class Timeshifts
{
  using Label = int;
  
public:

  using Timeshift_t = int64_t;
  using Timeshifts_t = std::vector<Timeshift_t>;

  Timeshifts() = default;
  Timeshifts(std::string const & filename) : m_filename(filename) {this -> load(m_filename);}
  Timeshifts(Timeshifts const & other) : m_timeshifts(other.m_timeshifts) {}
  
  Timeshifts& operator=(Timeshifts const & other)
  {
    m_timeshifts = other.m_timeshifts;
    m_ok = other.m_ok;
    return *this;
  }

  /// @brief Use this method to load timeshifts from a .dT file
  bool load(std::string const & filename);

  auto const & operator[] (int const & i) const {return m_timeshifts[i];}
#ifdef HIT_HPP // If Hit.hpp is included
  void operator() (Hit & hit) const {hit.stamp += m_timeshifts[hit.label];}
#endif //HIT_HPP

  Timeshifts_t const & get() const {return m_timeshifts;}
  Timeshifts_t const & data() const {return m_timeshifts;}
  auto const & get(int const & i) const {return m_timeshifts[i];}
  operator bool() const & {return m_ok;}

  void write(std::string const & fullpath, std::string const & name);

  auto size() const {return m_timeshifts.size();}

private:
  std::string m_filename;
  Timeshifts_t m_timeshifts;
  bool m_ok = false;
  int m_nb_detectors = 0;
  Colib::Path m_outPath;
  
  public:
  class NotFoundError
  {
  public:
    NotFoundError(std::string const & filename) : m_filename(filename) {error("Timeshiftor::NotFoundError", filename);}
    auto const & getFilename() const {return m_filename;}
  private:
    std::string m_filename;
  };

};

bool Timeshifts::load(std::string const & filename)
{
  std::ifstream inputFile(filename, std::ifstream::in);
  if (!inputFile.good()) {throw NotFoundError(filename);}
  else if (Colib::fileIsEmpty(inputFile)) {print("TIMESHIFT FILE", filename, "EMPTY !");return false;}
  std::string line = ""; // Reading buffer
  Label label = 0; // Reading buffer

  // ----------------------------------------------------- //
  // First extract the maximum label
  Label size = 0;
  // Infer the number of detectors from the higher label in calbration file
  while (getline(inputFile, line))
  { 
    std::istringstream iss(line);
    iss >> label;
    if (size<label) size = label;
  }
  size++; // The size of the vector must be label_max+1
  // if (detectors && size<detectors.size()) size = detectors.size(); // Ensure there is no mismatch with the detectors module
  inputFile.clear(); 
  inputFile.seekg(0, inputFile.beg);

  // ----------------------------------------------------- //
  // Now fill the vector
  m_timeshifts.resize(size, 0);
  Timeshift_t shift = 0;
  while (getline(inputFile, line))
  { // Then fill the array
    m_nb_detectors++;
    std::istringstream iss(line);
    iss >> label >> shift;
    m_timeshifts[label] = shift;
    shift = 0;
  }
  inputFile.close();
  print("Timeshifts extracted from", filename);
  return (m_ok = true);
}

void Timeshifts::write(std::string const & fullpath, std::string const & name)
{
  m_outPath = Colib::Path (fullpath, true);
  if (!m_outPath) {m_ok = false; return;}

  Colib::File outData (m_outPath+name);
  outData.setExtension(".dT");

  std::ofstream outTimeshiftsFile(outData, std::ios::out);
  
  for (size_t label = 0; label<m_timeshifts.size(); label++)
  {
    auto const & dT = m_timeshifts[label];
    if (dT != 0) outTimeshiftsFile << label << "\t" << dT << std::endl;
  }

  outTimeshiftsFile.close();

  print("Timeshifts data written to", outData);
}

#endif //TIMESHIFTS_HPP
