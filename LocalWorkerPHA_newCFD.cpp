#include "LocalWorkerPHA.h"

#include <cstdio>
#include <cstring>
#include <algorithm>

#include "rubuilder/utils/ConfReader.h"

// Content of CFD.hpp :
#ifndef CFD_HPP	
#define CFD_HPP

enum tick_lengths{four, ten};
constexpr auto tick_length = ten;
constexpr bool verboseCFD = true;

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <unordered_map>
#include <sstream>   
#include <map>
class CFD
{
protected:
  // Some helper functions (from libCo):
  using size_t = std::size_t;

  template <typename T, typename std::enable_if<std::is_floating_point<T>::value, bool>::type = true>
  inline static constexpr bool is_floating() noexcept { return true;}
  template <typename T, typename std::enable_if<!std::is_floating_point<T>::value, bool>::type = true>
  inline static constexpr bool is_floating() noexcept { return false;} 
  virtual inline double fast_uniform() noexcept {
    static thread_local std::minstd_rand generator(std::random_device{}());
    static thread_local std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
  }
  static constexpr const char* RED   = "\u001b[31m";
  static constexpr const char* RESET = "\u001b[0m" ;
  
  using CFD_t = std::vector<double>;

public:

  CFD() noexcept = default;

  CFD(std::string setting_file)
  {
    loadParameters(setting_file);
  }

  template<class T = double>
  CFD(std::vector<T> const & _trace, size_t nb_samples_baseline = 10) :
    m_size(_trace.size())
  {
    if (_trace.empty()) return;

    double baseline = _trace[0]; 
    if (nb_samples_baseline > 1)
    {
      for (size_t sample_i = 1; sample_i<nb_samples_baseline; ++sample_i) baseline += _trace[sample_i];
      baseline /= nb_samples_baseline;
    }
    
    trace.reserve(m_size);

    for (auto const & sample : _trace) 
    {
      if (is_floating<T>()) {
        trace.push_back(sample - baseline);
      }
      else { 
        trace.push_back(static_cast<double>(sample  + fast_uniform()) - baseline);
      }
    }
  }

  CFD& operator=(CFD const & other)
  {
    this -> trace = other.trace;
    this -> cfd = other.cfd;
    this -> m_size = other.m_size;
    return *this;
  }

  template<class T = double>
  void setTrace(std::vector<T> const & _trace)
  {
    trace.clear();
    if (_trace.empty()) return;

    double baseline = _trace[0]; 
    int nb_samples_baseline=10; 
    for (size_t sample_i = 1; sample_i<nb_samples_baseline; ++sample_i) baseline += _trace[sample_i];
    baseline /= nb_samples_baseline; 
  
    trace.reserve(m_size);

    for (auto const & sample : _trace) 
    {
      if (is_floating<T>()) {
        trace.push_back(sample - baseline);
      }
      else { 
        trace.push_back(static_cast<double>(sample  + fast_uniform()) - baseline);
      }
    }
    //this->trace = _trace;
    //return *this;
  }
  
  

  template<class T = double>
  CFD& operator=(std::vector<T> const & _trace)
  {
    if (_trace.empty()) return *this;
    m_size = _trace.size();

    double baseline = 0;
    // for (size_t sample_i = 0; sample_i<nb_samples; ++sample_i) baseline += _trace[sample_i];
    // baseline /= nb_samples;
    
    trace.reserve(m_size);

    for (auto const & sample : _trace) 
    {
      if (is_floating<T>()) {
        trace.push_back(sample - baseline);
      }
      else {
        trace.push_back(static_cast<double>(sample  + fast_uniform()) - baseline);
      }
    }
    return *this;
  }

  template<class T = double,typename Shift_Type=int,  typename Fraction_Type=double>
  CFD(std::vector<T> const & _trace, Shift_Type shift, Fraction_Type fraction) : CFD(_trace)
  {
    this -> calculate(shift, fraction);
  }
  template<class T = double,typename Shift_Type=int,  typename Fraction_Type=double>
  CFD(std::vector<T> const & _trace, Shift_Type shift, Fraction_Type fraction, size_t nb_samples_baseline) : CFD(_trace, nb_samples_baseline)
  {
    this -> calculate(shift, fraction);
  } 

  virtual void calculate(size_t shift, double fraction)
  {
    if (fraction>1.) {std::cout << RED << "in CFD(trace, shift, fraction): fraction>1 !!" << RESET << std::endl; return;} 
    cfd.clear(); 
    if (m_size < 2*shift) {std::cout << RED << "in CFD(trace, shift, fraction): m_size = " << m_size << " < 2*shift = " << 2*shift << " !!" << RESET << std::endl; return;} 
    cfd.reserve(m_size - 2*shift);  
    for (size_t bin = 2*shift; bin<m_size - shift; ++bin)\
    {
      auto const & value = fraction * trace[bin] - trace[bin - shift] ;
      cfd.push_back(value);
    }
  }

  template<typename Shift_Type=int,typename Fraction_Type=double >
  void calculate(Shift_Type const & shift, Fraction_Type const & fraction)
  {
    if (fraction>1.) {std::cout << RED << "in CFD(trace, shift, fraction): fraction>1 !!" << RESET << std::endl; return;}
    cfd.clear();
    if (m_size < 2*shift) {std::cout << RED << "in CFD(trace, shift, fraction): m_size = " << m_size << " < 2*shift = " << 2*shift << " !!" << RESET << std::endl; return;}
    cfd.reserve( static_cast<int>(m_size - std::floor( 2*shift)));
    
    // for (size_t bin = std::floor(2*shift); bin<m_size - shift; ++bin){
    //   auto const & value = std::is_integral_v<Shift_Type> ?  fraction * trace[bin] - trace[bin - shift]: fraction*trace[bin]-interpolate_next_value(  &trace[static_cast<int>(std::floor(bin-shift))], bin-shift- std::floor(bin-shift) ) ;
    //   cfd.push_back(value);
    // }
    for (size_t bin = std::floor(2*shift); bin<m_size - shift; ++bin){
      auto const & value = std::is_integral<Shift_Type>::value ?  fraction * trace[bin] - trace[bin - shift]: fraction*trace[bin]-interpolate_next_value(  &trace[static_cast<int>(std::floor(bin-shift))], bin-shift- std::floor(bin-shift) ) ;
      cfd.push_back(value);
    }
  } 

  virtual void calculate(int glabel)
  {
    //std::cout<<"glabel="<< glabel<<"\n";
    if(sType==LABEL)
    {
      if(key_found(sShifts,glabel))
    {
      calculate(sShifts[glabel],sFractions[glabel]);
    }
    else
    {
      std::cout<<"label glabel "<< glabel<< "has not had any settings loaded \n";
      calculate(0,0);
    }
    }
    else
    {
       if(sType==BOARD)
      {
        if(key_found(sShifts,glabel))
        {
          calculate(sShifts[glabel],sFractions[glabel]);
        }
        else
        {
          std::cout<<"label glabel "<< glabel<< "has not had any settings loaded \n";
          calculate(0,0);
        }
        calculate(sShifts[glabel],sFractions[glabel]);
      }
      else
      {
        std::cout<<"no settings have been loaded";
        calculate(0,0);
      }
    }
  
  }
  //=added
    template<typename T>
  double interpolate_next_value(const T*  value_0, double desired_index)
  {
    if(desired_index>=1)
    {
      std::cout<<"error, gap >= 1";
      return 0;
    }
    const T* value_1=value_0+1;
    return (*value_0) + (*value_1-*value_0)*desired_index;
  } 
  
   


  /// @brief Returns the zero crossing before the signal first crosses the threshold (unit : sample length)
  double findZero(double const & threshold)
  {
    if (threshold>0) std::cout << RED << "in CFD::findZero(threshold) : threshold > 0 !" << RESET << std::endl;
    for (size_t bin_i = 0; bin_i < cfd.size(); ++bin_i){  // Loop through the cfd values
      if (cfd[bin_i] < threshold){                        // The cfd value crossed the threshold
        for (size_t bin_j = bin_i; bin_j>0; --bin_j){     // Looping back for looking for the zero crossing
          if (cfd[bin_j] > 0) return interpolate0(bin_j); // Zero crossing found, return the interpolated zero crossing between samples before and after
        }
        return noZero; // The 0 crossing happened before the first sample, so impossible to determine it
      }
    }
    return noSignal; // The signal never crosses the threshold -> the signal is too small, the cfd parameters are wrong, or the threshold is too large
  }

  /// @brief Returns the last zero crossing before the trace reaches its minimum (unit : sample length)
  double findZero()
  {
    //1. Find the minimum of the trace :
    double minimum=cfd[0]; size_t min_index;
    for (size_t bin_i = 0; bin_i < cfd.size(); ++bin_i)
    {  // Loop through the cfd values
      if (cfd[bin_i] < minimum)
      {
        min_index=bin_i;
        minimum=cfd[bin_i];
      }
    }

    if(minimum>0) return noZero; // Trace is never below 0 -> can't determine 0 crossing (for obvious reason)

    // 2. Find the last zero crossing before the trace reaches its minimum :
    for(size_t bin_i=min_index;bin_i>=0;--bin_i)
    {
      if(cfd[bin_i]>0) return interpolate0(bin_i);
    }
    return noSignal; // Should never return this now.
  } 
  CFD_t trace;
  CFD_t cfd;

  // Static variables :
  constexpr static double noZero = 0.;
  constexpr static double noSignal = -1.;
 

  
  
protected:  
  double interpolate0(size_t const & bin) const noexcept 
  {
    auto const & cfd_0 = cfd[bin];
    auto const & cfd_1 = cfd[bin+1];
    if( cfd_0 == cfd_1) return 0; 
    else return bin - cfd_0 / (cfd_1 - cfd_0);
  }
 
  size_t m_size = 0;

  
public:
  
  enum Type {BOARD, LABEL, UNDEFINED};
  int sType = UNDEFINED;
  constexpr void setType(int type) noexcept {sType = type;}

  using Shifts     = std::unordered_map<int, int   >;
  //using Thresholds = std::unordered_map<int, double>;
  using Fractions  = std::unordered_map<int, double>;

  Shifts     sShifts     = {};
  //static inline Thresholds sThresholds = {};
  Fractions  sFractions  = {};
  
  template<typename K, typename V> 
  static inline bool key_found(std::unordered_map<K,V> const & map, K const & key)
  {
    typename std::unordered_map<K, V>::const_iterator it = map.find(key);
    return it != map.end();
  }


  void loadParameters(std::string filename)
  {
    std::ifstream paramFile(filename);
    std::string line;
    std::getline(paramFile, line);

    // Get the header. Anything can be written in it, but at least BOARD or LABELS
    // in order to know if the parameters are board wide and the label the BOARD_ID,
    // or detector per detector with the label the global label (BOARD_ID*16 + Channel_ID*2 + subchannel_ID) 
        if (line.find("BOARDS") != std::string::npos) sType = BOARD;
    else if (line.find("LABELS") != std::string::npos) sType =LABEL;
    else std::cout << CFD::RED << "CFD::loadParameters " << filename << " : format issue. Should begin with LABELS or BOARDS" << RESET << std::endl; 
    //=  should we replace the data?
    bool replace=false;
    while(std::getline(paramFile, line))
    {
      // Get the parameter of each board or each detector
      std::istringstream iss(line);
      int label; iss >> label;
    
      if( key_found(sShifts,label))
      {
      double tmp_d;
          iss >> tmp_d; sShifts    [label]= tmp_d;
          //iss >> tmp_d; sThresholds  [label]= tmp_d;
          iss >> tmp_d; sFractions   [label]= tmp_d;
      }
      else
      {
      double tmp_d;
          iss >> tmp_d; sShifts    .emplace(label, tmp_d);
          //iss >> tmp_d; sThresholds.emplace(label, tmp_d);
          iss >> tmp_d; sFractions .emplace(label, tmp_d);
      }
    }
  }

  //= added
  //= note:this will (obviously) overwrite whatever's already there so it's recommended to read the file first
  void saveParameters(std::string filename)
  {
    std::ofstream paramFile(filename);
    if(sType==BOARD){
    paramFile<< "BOARDS \n";
    }
      
    if(sType ==LABEL)
    {
    paramFile<< "LABELS";
    }   
    
    for (auto pair: sShifts)
    {
      int key=pair.first;
      paramFile<< key <<"\t"<<sShifts[key] << "\t"<<sFractions[key]<<"\n";
      //paramFile<< key <<"\t"<<sShifts[key] <<"\t" << sThresholds[key]<<"\t"<<sFractions[key]<<"\n";
      
      //keys.push_back(kv.first);
    }
    paramFile.close();
  }
  //= added
  void para_print_test()
  {
    if(sType=BOARD)
    {
      printf("Per board params \n");
    }
    if(sType == LABEL)
    {
      printf("Per Glabel params \n");
    }
    //std::vector<int> keys; keys.reserve(sShifts.size());
    for (auto pair: sShifts)
    {
      int key=pair.first;
      
      std::cout<<"label:" << key <<"\t shift:"<<sShifts[key] <<"\t fraction:"<<sFractions[key]<<"\n";
      //std::cout<<"label:" << key <<"\t shift:"<<sShifts[key] <<"\t threshold:" << sThresholds[key]<<"\t fraction:"<<sFractions[key]<<"\n";
      
      //keys.push_back(kv.first);
    }
    std::cout<<"End";
  }
};
#endif //CFD_HPP

const int   specCal     =  4096;    // adc range
const int   specRaw     = 16384;    // adc range
const int   matPackChan =    10;

using namespace std;

LocalWorkerPHA::LocalWorkerPHA(CycleServer& _cServer, int& _tcpPortno, int boardID,
			       int numdoms, int numsubs, int numsamp)
  : LocalWorker(_cServer, _tcpPortno, boardID)
{
  gDaughterClass = "LocalWorkerPHA";
  fOutMagicData = nullptr;
  calibration.clear();
  fGlobalTSOffset = 0L;
  fEnergyGain = 1.f;
  HasHeader   = false;
  FirstBlock  = true;
  outEvt = new caenEventPHA_t;
  tstCount = 0;
  fDomID = 0;
  tstoff = 0;
  v1730MaxChans=16;
  fOutMagicData = new uint32_t[16];
  fLastTS =  0;
  aTrace = nullptr;
  wspace = nullptr;

  waveAnalyzer = new MWD;
  waveAnalyzer->Reset(254, true);      // user-provided working buffers
  waveAnalyzer->nDcwin     = 254 / 5;     // the first 20% --> check if this is valid
  waveAnalyzer->nSmooth    = 5;
  waveAnalyzer->nWidthT    = 8;
  waveAnalyzer->nDelayCFD  = 20;
  waveAnalyzer->fFractCDF  = 0.25f;
  waveAnalyzer->fThreshTFA = 20;
  waveAnalyzer->fThreshCFD = 2;
  waveAnalyzer->nMinWidthP = 2;
  waveAnalyzer->nMinWidthN = 6;

  // cfdCalc=  CFD(boardID+"_settings.txt");
  

  memset(TimeShaper,  0, sizeof(TimeShaper));
  for (int ii = 0; ii < 16; ii++) {
    if (TimeShaper[ii]) {
      delete TimeShaper[ii];
      TimeShaper[ii] = nullptr;
    }
  }
}

// todo: add the delete []
LocalWorkerPHA::~LocalWorkerPHA()
{
  cServer.Finish();
  if(aTrace != nullptr)
    delete [] aTrace;
  if(wspace != nullptr)
    delete [] wspace;
}

int LocalWorkerPHA::GetParameters(bool doList, bool partial)
{
  float       cTstamp(0);
  std::string histPrefix;
  uint32_t    calBoard(0), calChan(0);
  float       calThresh(0), calOffs(0), calGain(1.f);
  uint16_t    boardID, boardFW, boardCH, boardTS, boardSampling;
  float       intercept, slope;
  std::string boardName;
  std::string graphServ;
  uint32_t    graphPort;
  int16_t     channelID;
  float       tsOffs;
  std::string tmpOutData;
  uint16_t    preTrig;
  std::string confFile;
  bool        AnaRiseTime;
  int         rdbPort;
  std::string rdbHost;
  hGroup.destroy();
  rubuilder::utils::ConfReader conf;
  conf.Add("SaveDataDir",      "where to write data and spectra"
	   , &fOutPath);
  conf.Add("SpecPrefix",       "prefix to names of spectra",
	   &boardID, &histPrefix); 
  conf.Add("Board",      "Board description: BoardID, Name, Firmware, number of channel, nsPerTimestamp, ns per sample",
	   &boardID, &boardName, &boardFW, &boardCH, &boardTS, &boardSampling);
  conf.Add("TSOffset",   "Timestamp offset for board X chanel Y",
	   &boardID, &channelID, &tsOffs);
  conf.Add("GraphiteServer", "Graphite Server to send the rates of the boards: Hostname and port",
	   &graphServ, &graphPort);
  conf.Add("OutDataKey", "Output ADF frame key (0xFA0201A2=SPIDER, 0xFA0201A3=DANTE, ...",
	   &boardID, &tmpOutData);
  conf.Add("PreTrigger",       "Pre-trigger time in ns (expected to be common for the full board)",
	   &boardID, &preTrig);
  conf.Add("ConfFilePSA",      "Configuration file for the PSA",
	   &boardID, &confFile);
  conf.Add("AnaRiseTime",      "Enable analysis of the short traces for the particle discrimination",
	   &boardID, &AnaRiseTime);
  conf.Add( "RedisServer"    ,"REDIS Server host and port",
	    &rdbHost, &rdbPort);
  conf.Add( "CalShort"    ,"Calibration coefficient for short intergral: board channel intercept slope",
	    &boardID, &channelID, &intercept, &slope);
  conf.Add( "CalLong"    ,"Calibration coefficient for long intergral: board channel intercept slope",
	    &boardID, &channelID, &intercept, &slope);
  // refine behaviour
  
  if (doList) { conf.Show(gMotherClass); return 0; }
  
  if(partial) {
    // disable "non-calibration" keywords
    conf.Disable("Domains");
    conf.Disable("ServerPort");
    conf.Disable("SaveDataDir");
    conf.Disable("SpecPrefix");
    conf.Disable("NumBoards");
  }
  
  int rv = conf.Read(gMotherClass, fConfigFile);
  
  if (rv != 0)
    return rv;  // error
  
  int which, times;
  doRise = false;
  if(conf.Seen("AnaRiseTime")) doRise = true;
  which = conf.Find("Board");
  times = conf.Times(which);
  for (int nn = 0; nn < times; nn++) {
    conf.Restore(which, nn);
    if(boardID!=fBoardID)continue;
    filterName << "LF_"<< boardName;
    ratesMonitor.setName(boardName);
    fNsPerSample  = boardSampling;
    fNsPerTStamp  = boardTS;
    nbChannel     = boardCH;
    
    for(int ch = 0 ; ch < nbChannel ; ++ch){
      lastTS[ch]  =0;
      nbEvents[ch]=0;
      nbGoodEvents[ch]=0;
      nbIdleEvents[ch]=0;
      nbTraces[ch]=1000;
    }
  }  

  which = conf.Find("SpecPrefix");
  times = conf.Times(which);
  for (int nn = 0 ; nn < times ; ++nn){
    conf.Restore(which,nn);
    if(boardID!=fBoardID)continue;
    fHistPfx = histPrefix;
    specAmpl = new MultiHist<uint32_t>(2,nbChannel, 16384);     // [0]-fromShaper [1]-fromTrace
    specAmpl->setFileName(fOutPath+fHistPfx+"?Ampl.spec");
    specAmpl->setComment("raw energy");
    hGroup.add(specAmpl);
  }  

  if(conf.Seen("GraphiteServer")){
    ratesMonitor.setHost(graphServ,graphPort);
  }

  which = conf.Find("OutDataKey");
  times = conf.Times(which);
  for (int nn = 0; nn < times; nn++) {
    conf.Restore(which, nn);
    if(boardID!=fBoardID)continue;
    outDataKey = OutDataKey(tmpOutData) ;
  }    
 
  which = conf.Find("TSOffset");
  times = conf.Times(which);
  for (int nn = 0; nn < times; nn++) {
    conf.Restore(which, nn);
    if(boardID!=fBoardID) continue;
    if(channelID==-1){
      for(int ch = 0 ; ch < nbChannel ; ++ch){
	tsOffset[ch]=tsOffs;
      }
    }else{
      tsOffset[channelID]=tsOffs;
    }
  }

  fPromptPeriod = max(0, fPromptPeriod);
  cServer.SetHistGroup(&hGroup);

  return rv;
}

int LocalWorkerPHA::ProcessConfig(const char* fname)
{
  cout << gMotherClass + "::process_config() called with file = " << fname << endl;
  
  string confString;
  fConfigFile = fname;
  
  if (!fConfigFile.empty()) {          // use defaults, if process_config was not called
    if (!fileExists(fConfigFile)) {
      return 102;              // fatal error because configuration file not present
    }
    
    // VME/AGAVA
    int error_code = GetParameters();
    if (error_code) {
      return error_code;
    }
    
    fScaleTime = 0;
    if (fScaleDown > 0)
      fScaleTime = (tScaleMask + 1) / fScaleDown;
  }
  
  cServer.SetCommandFile(fConfigFile + ".live");
#ifdef GAL_OFFLINE
  cServer.AcceptExit(); // enable ECommand::kExit
#endif
  cServer.Start(filterName.str().c_str(), max(fPromptPeriod, 0));
  
  return 0;
}

int LocalWorkerPHA::ProcessStart(uint32_t run)
{
  fRunNumber = run;
  fCountBlocks = 0;
  cServer.Start(filterName.str().c_str(), max(fPromptPeriod, 0));
  return 0;
}

int LocalWorkerPHA::ProcessStop()
{
  cServer.Save(true);
  
  cout << "Number of events in board " << fBoardID << " of this run " << cServer.GetCounts() << endl;
  printStatistics();
   
  return 0;
}

int LocalWorkerPHA::ProcessBlock(const void *input_buffer,  int input_size,  int packet_ID,
				 void *output_buffer, int output_size, int &used_output)
{
  uint32_t* inpBuffer   = (uint32_t*)input_buffer;
  uint32_t* outBuffer   = (uint32_t*)output_buffer;
  uint32_t  inpBuffSize = input_size  / sizeof(uint32_t);

  //  std::cout << "Receiving data in  " << __FILE__ << std::endl;
  int error_code = 0;
  uint32_t posBuffer =0;
  std::bitset<32> decodeField;
  int outByteSize = used_output;
  int countEvents = 0;
  uint32_t outEvtLen = (sizeof(agataKey) + sizeof(subDataPHA_t));
  // Aggregate buffer has a header containing different information
  uint32_t aggLength         = inpBuffer[posBuffer++]&0x0FFFFFFF;
  decodeField                = inpBuffer[posBuffer++];
  uint32_t board             = project_range<27,31>(decodeField).to_ulong();
  bool     boardFailFlag     = decodeField.test(26);
  std::bitset<8> channelMask = project_range<0,7>(decodeField).to_ulong();
  uint32_t aggregateCounter  = inpBuffer[posBuffer++]&0x7FFFFF;
  uint32_t aggregateTimeTag  = inpBuffer[posBuffer++];
  uint8_t couples[8]={0,0,0,0,0,0,0,0};
  // we need to reconstruct the channel number from the couple id for the 
  // 16 channels boards
  int index = 0 ;
  for(uint16_t bit = 0 ; bit < 8 ; bit++){ 
    if(channelMask.test(bit)){
      couples[index]=bit;
      index++;
    }
  }

  for(int ch=0; ch < nbChannel ; ch++)
  {  
    if(lastTS[ch]==0) lastTS[ch]=fGlobalTSOffset;
  }

  uint32_t coupleAggregateSize = 0 ;
  uint64_t tstamp_4ns          = 0 ;
  uint64_t tstamp_10ns         = 0 ;
  uint32_t fineTS              = 0 ;
  int      offsetTS            = 0 ;
  uint32_t flags               = 0 ;
  bool     pur                 = false;
  bool     satu                = false;
  bool     lost                = false;
  bool     trap_sat            = false;
  bool     coinc               = false;
  uint32_t chanNum             = 0 ;
  float    cfd                 = 0.;
  // here starts the real decoding of the data for the channels 
  for(uint8_t cpl = 0 ; cpl < channelMask.count() ; cpl++)
  {
    coupleAggregateSize = inpBuffer[posBuffer++]&0x3FFFFF;
    dataForm.SetDataFormat(inpBuffer[posBuffer++]);
    //dataForm.printDataFormat();
    const uint16_t evtSize = dataForm.evtSize();
    if(evtSize==0){
      dataForm.printDataFormat();
      break;
    }

    for(uint16_t evt = 0 ; evt < (coupleAggregateSize-2)/evtSize;evt++)
    {
      if(posBuffer==inpBuffSize) break;
      float trig = 0.;
      cfd = 0.;
      decodeField = inpBuffer[posBuffer++];
      tstamp  = project_range<0,30>(decodeField).to_ulong();
      chanNum = static_cast<uint32_t>(project_range<31,31>(decodeField).to_ulong()
				      +couples[cpl]*2);
      
      if(dataForm.enaTrace())
      {
        std::valarray<float> samples;
        samples.resize(dataForm.numSamples());
        for(uint16_t s = 0 ; s < dataForm.numSamples()/2 ; ++s){
          decodeField = inpBuffer[posBuffer++];
          samples[s*2]   = project_range< 0,13>(decodeField).to_ulong() + (rand()%10000)/10000.-0.5;
          samples[s*2+1] = project_range<16,29>(decodeField).to_ulong() + (rand()%10000)/10000.-0.5;
        }

        { // Corentin and Piotr modifications :
          if((fBoardID<2 && chanNum%2==0))// To target only eagle HPGe
          { 
            std::vector<double> traceVec(std::begin(samples), std::end(samples));
            static constexpr size_t nbSamplesBaseline = 20;
            CFD cfdModule(traceVec, nbSamplesBaseline);
            // 1. With default parameters (To remove)
            static constexpr int    shift = 9  ;
            static constexpr double fract = 0.7;
            cfdModule.calculate(shift, fract);
            // 2. With parameter file (TBD)
            // auto const & global_label = board*16+chanNum;
            // cfdModule.calculate(global_label);
            // cfd = cfdModule.findZero()*fNsPerSample*1000; // in ps
            cfd = cfdModule.findZero(); // in number of ticks
            if (verboseCFD && cfd == CFD::noSignal) std::cerr << "label " << board*16+chanNum << " : no signal" << std::endl;
            if (verboseCFD && cfd == CFD::noZero  ) std::cerr << "label " << board*16+chanNum << " : no zero"   << std::endl;
          }
        }
        
        /* Corentin and Piotr comnmenting out
        //IP_start      //only for eagle, try:
        if((fBoardID<0 && chanNum%2==0)){
          //std::cout<<"AAAAAAAAAAAAAAAAA::::   "<<fBoardID<<std::endl;

          float baseline = static_cast<std::valarray<float>>(samples[std::slice(0,60,1)]).sum()/60.;
          samples-=baseline;
          if(samples.sum() < 0){
            samples *=-1;
          }
          
          if(aTrace == nullptr)
            aTrace = new float[samples.size()];
          if(wspace == nullptr)
            wspace = new float[3 * samples.size()];
          
          std::copy(begin(samples),end(samples),aTrace);
          float otrig = 0 ;
          //float *wspace = new float[3 * samples.size()];
          trig = waveAnalyzer->TriggerCFD(aTrace, 12, wspace, otrig);//was 5
          if (TimeShaper[chanNum] == nullptr && nbTraces[chanNum]>0) {
            std::ostringstream on;
            on << fOutPath + fHistPfx << "__1000-4-"
              << samples.size() << "-F__timing-"
              << setfill('0') << setw(2) << chanNum << ".dat";
            TimeShaper[chanNum] = new Append<float>(on.str(), 4*samples.size());
          }
          if(nbTraces[chanNum]>0){
            float *fw0 = wspace;
            float *fw1 = fw0 + samples.size();
            float *fw2 = fw1 + samples.size();
            TimeShaper[chanNum]->add(aTrace,samples.size());
            TimeShaper[chanNum]->add(fw0,samples.size());
            TimeShaper[chanNum]->add(fw1,samples.size());
            TimeShaper[chanNum]->add(fw2,samples.size());
          nbTraces[chanNum]--;
          }
          //	  if(chanNum%2==0){
            cfd = otrig*fNsPerSample*1000;
            // to have it in ps units
            //}
          //std::cout << "trig = " << trig << " otrig = " << otrig << std::endl; 
        
        }*/
        //IP_end
      }

      if(dataForm.enaExtras())
      { // decoding of EXTRAS2
        decodeField = inpBuffer[posBuffer++];
        switch(dataForm.confExtras()) 
        {
          case 0: // Extended Time Stamp [31:16] ; baseline*4 [15:0]
            tstamp = (static_cast<uint64_t>(project_range<16,31>(decodeField).to_ulong())<<31) 
              | (uint64_t)tstamp;
            break;
          case 1: // Reserved??
            break;
          case 2: // Extended Time stamp [31:16] ; Reserved [15:10] ; Fine Time Stamp [9:0]
            fineTS = static_cast<uint16_t>(project_range< 0, 9>(decodeField).to_ulong());
            outEvt->theData.cfd = fNsPerTStamp/1.024* fineTS; // to have it in ps units
            
            tstamp = (static_cast<uint64_t>(project_range<16,31>(decodeField).to_ulong())<<31) 
              | (uint64_t)tstamp;
            break;
          case 3:// Reserved
            break;
          case 4: // Lost Trigger Counter [31:16] ; Total Trigger [15:0]
            break;
          case 5: // Positive zero crossing [31:16] ; Negative zero crossing [15:0]
            break;
          case 7: // Reserved
            break;
          default:
            break;
        }
      }
      
      decodeField = inpBuffer[posBuffer++];
      // the last word contains energy[0-14], PU [15] and extra[16:25]
      satu     = decodeField.test( 4+16);
      lost     = decodeField.test( 5+16);
      coinc    = decodeField.test( 8+16);
      pur      = decodeField.test( 9+16);
      trap_sat = decodeField.test(10+16);

      //	statsData.time[cpl] = tstamp*2/pow( 10, 9 );   // must be in ns

      if(tsOffset.find(chanNum) != tsOffset.end()) offsetTS = tsOffset[chanNum];
      else                                         offsetTS = 0;
      outEvt->theKey.size     = outEvtLen              ; 
      outEvt->theKey.key      = outDataKey             ;                      // default is data:ranc1

      /* Corentin and Piotr commenting out :

      tstamp_10ns             = tstamp * fNsPerTStamp / 10 ;
      cfd                    += (tstamp*fNsPerTStamp - tstamp_10ns*10)*1000; // in ps no to loose precision when using uint16_t
      tstamp_10ns            += offsetTS + fGlobalTSOffset;
      	std::cout << __func__ << ":" << fGlobalTSOffset << std::endl;

      GJ - also here has to be - eagle only:
           if(dataForm.enaTrace()){
      if(dataForm.enaTrace() && fBoardID<0 && chanNum%2==0 )
      {
        if(cfd>10000){
          int cfdOffset = std::floor(cfd)/10000.;
          tstamp_10ns         += cfdOffset;
          cfd = cfd - cfdOffset*10000.;
        }
      }
      else
      {
        cfd += outEvt->theData.cfd ;
        if(cfd>10000)
        {
          int cfdOffset = std::floor(cfd)/10000.;
          tstamp_10ns         += cfdOffset;
          cfd = cfd - cfdOffset*10000.;
        }
      }

      */
      
      { // Corentin and Piotr modifications :
        //1. Replacing the original 10 ns recalculated ticks by the raw 4 ns ticks
        if (tick_length == four)
        {
          tstamp_4ns = tstamp;
          tstamp_4ns += static_cast<uint64_t>((offsetTS + fGlobalTSOffset)/fNsPerSample*10.); // Recalculate the global offsets from the I suppose 10ns ticks to the 4 ns
        }
        else if (tick_length == ten)
        { // Keeping 10ns-long ticks
          tstamp_10ns = tstamp * fNsPerSample / 10.;
          tstamp_10ns += offsetTS + fGlobalTSOffset;
        }

        //2. Modify the cfd correction to only be a correction in the [0;fNsPerSample] ns range (typically [0,4]ns) :
        if((fBoardID<2 && chanNum%2==0))// To target only eagle HPGe
        {
          if (tick_length == four)
          {
            auto const & floor_4ns = std::floor(cfd);
            uint64_t nb_4ns_ticks = static_cast<uint64_t>(floor_4ns); // Get the integer part of the cfd correction
            tstamp_4ns += nb_4ns_ticks; // Correct the absolute timestamp with the integer number of ticks 
            cfd -= floor_4ns;     // Adjusts the cfd correction to only be the fractionnal part
            cfd *= fNsPerSample * 1000; // Convert the cfd correction to ps for respecting the output format
          }
          else if (tick_length == ten)
          {
            cfd *= fNsPerSample / 10.; // Convert the cfd correction in 10ns-long ticks for compatibility
            auto const & floor_10ns = std::floor(cfd);
            uint64_t nb_10ns_ticks = static_cast<uint64_t>(floor_10ns); // Get the integer part of the cfd correction
            tstamp_10ns += nb_10ns_ticks; // Correct the absolute timestamp with the integer number of ticks 
            cfd -= floor_10ns;       // Adjusts the cfd correction to only be the fractionnal part
            cfd *= 10000;                  // Convert the cfd correction to ps for respecting the output format
          }
        }
      }

      outEvt->theData.cfd     = static_cast<uint16_t>(cfd);
      if (tick_length == four)
      {
        outEvt->theKey.ts_0     = uint32_t(tstamp_4ns & 0xFFFFFFFF);         // tstLow
        outEvt->theKey.ts_1     = uint32_t((tstamp_4ns >> 32) & 0x0000FFFF); // tstHigh
      }
      else if (tick_length == ten)
      {
        outEvt->theKey.ts_0     = uint32_t(tstamp_10ns & 0xFFFFFFFF);         // tstLow
        outEvt->theKey.ts_1     = uint32_t((tstamp_10ns >> 32) & 0x0000FFFF); // tstHigh
      }
      outEvt->theData.energy  = static_cast<uint16_t>(project_range<0,15>(decodeField).to_ulong());

      /// DEBUG:
      //      if((fBoardID<3 && chanNum%2==0)){
      //	std::cout<<"LF data to be adf: b "<<fBoardID<<" ch "<<chanNum<<" -- energy  "<<outEvt->theData.energy<<std::endl;
      //      }      

      if(HCHECK(specAmpl))
      {
        if(!coinc) specAmpl->IncrQ(0,chanNum,outEvt->theData.energy/4,0); // Energy is on 16 bit.
        else       specAmpl->IncrQ(1,chanNum,outEvt->theData.energy/4,0); // Energy is on 16 bit.
      }
      
      // Adding PUR to the flag

      decodeField = (board << 8 | chanNum);
     
      if(ratesMonitor.values.find(chanNum) == ratesMonitor.values.end()) ratesMonitor.values[chanNum].Init();
      ratesMonitor.values[chanNum].time = tstamp*fNsPerTStamp;
      ++ratesMonitor.values[chanNum].totalEvents ;
      if( pur )      { ++ratesMonitor.values[chanNum].pileEvents ; decodeField.set(16); }
      if( trap_sat ) { decodeField.set(17); }
      if( satu )     { ++ratesMonitor.values[chanNum].satuEvents ; decodeField.set(18); }
      if( lost )     { ratesMonitor.values[chanNum].lostEvents += 1024; decodeField.set(19);}
      if( boardFailFlag ) {decodeField.set(20); }
	
      // N is set from bits[17:16] of register 0x1nA0 (default value 1024)
      // Adding flags to the evt num

      outEvt->theKey.evnum    = decodeField.to_ulong();
      
      outEvt->theKey.evnum   |= boardFailFlag;
      ++nbEvents[chanNum];
      countEvents++;
      if(outEvt->theData.energy<2)continue;
      ++nbGoodEvents[chanNum];
      lastTS[chanNum] = tstamp_10ns;
      if (fLastTS < tstamp_10ns) fLastTS = tstamp_10ns;
      int outEvtSize = SetOutput(outBuffer,output_size-outByteSize,outEvt);
      if(outEvtSize<0)
      { //output buffer full
        std::cout << "Returning to LF" << std::endl;
        return outEvtSize;
      }
      outByteSize += outEvtSize ;
      outBuffer   += outEvtSize / sizeof(uint32_t) ;
    }
  }
   
  if(error_code) return error_code;
  
  used_output = outByteSize;
  fCountBlocks++;
  
  cServer.Exec(tstamp_10ns, countEvents);
  cServer.Prompt(0, countEvents, uint32_t(outByteSize)*sizeof(uint32_t));
  
  if (cServer.CheckCommand()) 
  {
    if (cServer.CheckCommand(CycleServer::ECommand::kExit))      return -1000;
    if (cServer.CheckCommand(CycleServer::ECommand::kConfigure)) GetParameters(false, true);
  }
  return 0;  
}

int LocalWorkerPHA::SetOutput(uint32_t *outBuff, int availableBytes, const caenEventPHA_t * pE)
{
  // agataKey + numSubs * subData_t
  int outEvLenByte = sizeof(agataKey) + sizeof(subDataPHA_t); // max needed size

  
  agataKey* psaKey = (agataKey *)outBuff;
  outEvLenByte = sizeof(agataKey);                      // recalculate effective value
  
  // hand-made agata key
  psaKey->size  = pE->theKey.size;
  psaKey->key   = pE->theKey.key;
  psaKey->ts_0  = pE->theKey.ts_0;
  psaKey->ts_1  = pE->theKey.ts_1;
  psaKey->evnum = pE->theKey.evnum;
  int evtSize = (sizeof(agataKey) + sizeof(subDataPHA_t));
  if ( evtSize > availableBytes) {
    return availableBytes - outEvLenByte;
  }  
  
  outBuff += sizeof(agataKey) / sizeof(uint32_t);       // immediatly after agataKey
  
  // the full size as defined in subData_t
  subDataPHA_t pS = pE->theData;
  memcpy(outBuff, &pS, sizeof(subDataPHA_t));
  outEvLenByte += sizeof(subDataPHA_t);                  // recalculate effective value
  outBuff      += sizeof(subDataPHA_t) / sizeof(uint32_t);
  psaKey->size  = outEvLenByte;
  return outEvLenByte;
}

int LocalWorkerPHA::GenerateIdles(uint64_t largest_ts,
				  void *output_buffer, int output_size, int &used_output)
{
  uint32_t* outBuffer   = (uint32_t*)output_buffer;
  int outByteSize = used_output;
  uint32_t outEvtLen = (sizeof(agataKey) + sizeof(subDataPHA_t));  
  uint64_t newTS =0 ;
  
  for(int ch=0; ch < nbChannel ; ch++){  
    // if difference larger than 50ms for individual channels 
    // we keep 50 ms to be larger than the autoflush time - hoping that it is set to the minimum value possible
    if(lastTS[ch]<fGlobalTSOffset) lastTS[ch]=fGlobalTSOffset;
    if ( largest_ts - lastTS[ch] < 100000 ) continue; 
    ++nbIdleEvents[ch];
    //at this point we generate an idle event using the eventIdx to identify them       
    newTS = lastTS[ch] + 10000; 
    outEvt->theKey.size     = outEvtLen;
    outEvt->theKey.key      = outDataKey;
    outEvt->theKey.evnum    = 0x80000000 | ((fBoardID&0xFF) << 8) | (ch&0xff) ; // 32 bit set to high as idle flag
    outEvt->theKey.ts_0     = uint32_t(newTS & 0xFFFFFFFF);         // tstLow
    outEvt->theKey.ts_1     = uint32_t((newTS >> 32) & 0x0000FFFF); // tstHigh
    outEvt->theData.energy  = 0;
    outEvt->theData.cfd     = 0;
    int outEvtSize = SetOutput(outBuffer,output_size-outByteSize,outEvt);
    if(outEvtSize<0){ //output buffer full
      return outEvtSize;
    }
    outByteSize += outEvtSize ;
    outBuffer   += outEvtSize / sizeof(uint32_t) ;
    lastTS[ch] = newTS;
  }
  used_output = outByteSize;
  return 0;
}
