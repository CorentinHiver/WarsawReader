#pragma once

#define Nuball2

#ifdef CoMT
  std::mutex Nuball2mutex;
#endif // CoMT

#include <iostream>
#include <charconv> // Required for from_chars
#include <string_view>
#include <system_error>

#include "units.hh"
#include "libCo.hpp"
#include "RtypesCore.h"

//////////////////
/// Data types ///
//////////////////

using Label     = ushort;    // Label (ushort)
using ADC       = int;       // ADC (int)
using NRJ       = Energy_t;  // Energy in keV (float)
using Timestamp = ULong64_t; // Timestamp in ps (absolute)
using Time_ns   = float;     // Time in ns (relative) !deprecated! 
using Pileup    = bool;      // Pileup bit (bool) !unused!
using Index     = uint8_t;   // Used in analysis structures (Clovers, Paris...). Be careful to the max value of 255 !
// using Time      = int64_t;    // Time in ps (relative) Already declared in units.hh
// using Time      = int64_t;    // Time in ps (relative) Already declared in units.hh

////////////////////
/// Data vectors ///
////////////////////

using Label_vec   = std::vector<Label  >; // Vector of Label (ushort)
using ADC_vec     = std::vector<ADC    >; // Vector of ADC (int)
using Energy_vec  = std::vector<NRJ    >; // Vector of Energy in keV (float)
using Time_vec    = std::vector<Time   >; // Vector of Time in ps (relative)
using Time_ns_vec = std::vector<Time_ns>; // Vector of Time in ns (relative) !deprecated!  
using Pileup_vec  = std::vector<Pileup >; // Vector of Pileup bit (bool) !unused!

//////////////////
/// Data casts ///
//////////////////

/// @brief Casts a number into unsigned Label
template<typename T,  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
constexpr inline Label Index_cast(T t) {return static_cast<Index>(t);}

/// @brief Casts a number into unsigned Label
template<typename T,  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
constexpr inline Label Label_cast(T t) {return static_cast<Label>(t);}

/// @brief Casts a number into unsigned ADC
template<typename T,  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
constexpr inline ADC ADC_cast(T t) {return static_cast<ADC>(t);}

/// @brief Casts a number into unsigned Time
template<typename T,  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
constexpr inline Time Time_cast(T t) {return static_cast<Time>(t);}

/// @brief Casts a number into unsigned NRJ
template<typename T,  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
constexpr inline NRJ NRJ_cast(T t) {return static_cast<NRJ>(t);}

/// @brief Casts a number into unsigned Time_ns
/// @deprecated
template<typename T,  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
constexpr inline Time_ns Time_ns_cast(T t) {return static_cast<Time_ns>(t);}

/// @brief Casts a number into unsigned Timestamp
template<typename T,  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
constexpr inline Timestamp Timestamp_cast(T t) {return static_cast<Timestamp>(t);}

/// @brief Casts a number into unsigned Timestamp
template<typename T,  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
constexpr inline ULong64_t ULong64_cast(T t) {return static_cast<ULong64_t>(t);}

/// @brief Casts a number into unsigned Timestamp
template<typename T,  typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
constexpr inline Long64_t Long64_cast(T t) {return static_cast<Long64_t>(t);}

/////////////////
// FILE NAMING //
/////////////////

inline std::string runName     (std::string runPath) {return Colib::removeExtension(Colib::removePath(Colib::removeLast(runPath, "/")));}
inline std::string runNumberStr(std::string runName) {return Colib::split(runName, '_')[1];}
inline int         runNumber   (std::string runName) 
{
  auto nbStr = runNumberStr(runName);
  if (Colib::isNumber(nbStr))return std::stoi(nbStr);
  else return -1;
}

int rootRunNumber(std::string fullpath)
{
  print(Colib::removeAll(Colib::removeAll(Colib::removePath(fullpath), "run_"),"_rf.root"));
  return std::stoi(Colib::removeAll(Colib::removeAll(Colib::removePath(fullpath), "run_"),"_rf.root"));
}

////////////
// Values //
////////////

constexpr inline auto LabelMax = std::numeric_limits<Label>::max();

namespace NSI136
{
  static constexpr size_t LUT_s = 1000; // Look-up-table size

  static const auto U_runs  = Colib::computeList<int> (130, [](int run_number) {return 74 < run_number && run_number < 123;});
  static const auto Th_runs = Colib::computeList<int> (130, [](int run_number) {return 13 < run_number && run_number <  75;});
  
  static const std::unordered_map<std::string, std::vector<int>> both_runs {{"U", U_runs},{"Th", Th_runs}};

  // CLOVERS //

  static constexpr auto cloverIndex  = Colib::LUT<LUT_s> ([](Label label) -> Index {return Index_cast((label-23)/6);});
  static constexpr auto crystalIndex = Colib::LUT<LUT_s> ([](Label label) -> Index {return Index_cast((label-23)%6);});
  static constexpr auto crystalIndexGe = Colib::LUT<LUT_s> ([](Label label) -> Index {return crystalIndex[label]-2;});

  static constexpr auto isClover     = Colib::LUT<LUT_s> ([](Label label) -> bool  {return 22 < label && label < 168;});
  static constexpr auto isGe         = Colib::LUT<LUT_s> ([](Label label) -> bool  {return isClover[label] && crystalIndex[label] > 1;});
  static constexpr auto isBGO        = Colib::LUT<LUT_s> ([](Label label) -> bool  {return isClover[label] && crystalIndex[label] < 2;});
  static constexpr auto isR2         = Colib::LUT<LUT_s> ([](Label label) -> bool  {return isClover[label] && cloverIndex [label] > 12;});
  static constexpr auto isR3         = Colib::LUT<LUT_s> ([](Label label) -> bool  {return isClover[label] && cloverIndex [label] < 13;});
  
  static constexpr auto CloverRingNumber = Colib::LUT<LUT_s> ([](Label label) -> Index {if (isR2[label]) return 2; else if (isR3[label]) return 3; else return 0;});
  static constexpr auto CloverRingIndex  = Colib::LUT<LUT_s> ([](Label label) -> Index {return (isClover[label]) ? Index_cast(cloverIndex[label] % 12) : -1;});
  
  //   RF   //

  static constexpr auto isRF = Colib::LUT<LUT_s> ([](Label label) -> bool {return label == 251;});
  
  //   REF   //

  static constexpr auto isRef = Colib::LUT<LUT_s> ([](Label label) -> bool {return label == 252;});

  //  PARIS  //

  // static constexpr auto isBR1   = Colib::LUT<LUT_s> ([](Label label) -> bool  {return 200 < label && label < 209;});
  static constexpr auto isBR2   = Colib::LUT<LUT_s> ([](Label label) -> bool  {return 300 < label && label < 317;});
  static constexpr auto isBR3   = Colib::LUT<LUT_s> ([](Label label) -> bool  {return 400 < label && label < 413;});
  static constexpr auto isFR1   = Colib::LUT<LUT_s> ([](Label label) -> bool  {return 500 < label && label < 509;});
  static constexpr auto isFR2   = Colib::LUT<LUT_s> ([](Label label) -> bool  {return 600 < label && label < 617;});
  static constexpr auto isFR3   = Colib::LUT<LUT_s> ([](Label label) -> bool  {return 700 < label && label < 713;});
  static constexpr auto isBack  = Colib::LUT<LUT_s> ([](Label label) -> bool  {return Colib::lut_OR(label, isBR2, isBR3);});
  static constexpr auto isFront = Colib::LUT<LUT_s> ([](Label label) -> bool  {return Colib::lut_OR(label, isFR1, isFR2, isFR3);});
  static constexpr auto isParis = Colib::LUT<LUT_s> ([](Label label) -> bool  {return Colib::lut_OR(label, isBR2, isBR3, isFR1, isFR2, isFR3);});

  //  Dssd  //
  
  static constexpr auto isSector1 = Colib::LUT<LUT_s> ([](Label label) {return 799 < label && label < 816;});
  static constexpr auto isSector2 = Colib::LUT<LUT_s> ([](Label label) {return 819 < label && label < 836;});
  static constexpr auto isRing    = Colib::LUT<LUT_s> ([](Label label) {return 839 < label && label < 856;});
  static constexpr auto isSector  = Colib::LUT<LUT_s> ([](Label label) {return Colib::lut_OR(label, isSector1, isSector2);});
  static constexpr auto isDssd    = Colib::LUT<LUT_s> ([](Label label) {return Colib::lut_OR(label, isRing, isSector);});
  
  ////////////
  // Labels //
  ////////////

  static const auto labelsGe     = Colib::computeList<Label> (LUT_s, [](Label label){return (isGe    [label]);});
  static const auto labelsBGO    = Colib::computeList<Label> (LUT_s, [](Label label){return (isBGO   [label]);});
  static const auto labelsClover = Colib::computeList<Label> (LUT_s, [](Label label){return (isClover[label]);});
  static const auto labelsParis  = Colib::computeList<Label> (LUT_s, [](Label label){return (isParis [label]);});
  static const auto labelsBR     = Colib::computeList<Label> (LUT_s, [](Label label){return (isBack  [label]);});
  static const auto labelsFR     = Colib::computeList<Label> (LUT_s, [](Label label){return (isFront [label]);});
  static const auto labelsDssd   = Colib::computeList<Label> (LUT_s, [](Label label){return (isDssd  [label]);});

  static const auto labelsAll = Colib::append(labelsClover, Label_vec({251, 252}), labelsParis, labelsDssd);

  ///////////
  // Index //
  ///////////
  
  static constexpr auto GeIndex      = Colib::LUT<LUT_s> ([](Label label) -> Label {return (isGe[label ]) ? crystalIndex[label]-2+4*cloverIndex[label] : LabelMax;});
  static constexpr auto BGOIndex     = Colib::LUT<LUT_s> ([](Label label) -> Label {return (isBGO[label]) ? crystalIndex[label]+2*cloverIndex[label] : LabelMax;});
  static constexpr auto DssdIndex    = Colib::LUT<LUT_s> ([](Label label) -> Label 
  {
    if (!isDssd[label]) return LabelMax;
         if (isSector1[label]) return label - 800;
    else if (isSector2[label]) return label - 804;
    else if (isRing   [label]) return label - 808;
    else return LabelMax;
  });

  static constexpr auto ParisRingIndex = Colib::LUT<LUT_s> ([](Label const label) -> Label  
  {
    if (!isParis[label]) return LabelMax;
         if (isBR2[label]) return label - 301;
    else if (isBR3[label]) return label - 401;
    else if (isFR1[label]) return label - 501;
    else if (isFR2[label]) return label - 601;
    else if (isFR3[label]) return label - 701;
    else return LabelMax;
  });
  static constexpr auto ParisIndex = Colib::LUT<LUT_s> ([](Label const label) -> Label  
  {
      if (!isParis[label]) return LabelMax;
         if (isBR2[label]) return ParisRingIndex[label];
    else if (isBR3[label]) return ParisRingIndex[label] + Colib::lutEntries(isBR2);
    else if (isFR1[label]) return ParisRingIndex[label] + Colib::lutsEntries_sum(isBR2, isBR3);
    else if (isFR2[label]) return ParisRingIndex[label] + Colib::lutsEntries_sum(isBR2, isBR3, isFR1);
    else if (isFR3[label]) return ParisRingIndex[label] + Colib::lutsEntries_sum(isBR2, isBR3, isFR1, isFR2);
    else return LabelMax;
  });
  static constexpr auto ParisClusterIndex = Colib::LUT<LUT_s> ([](Label const label) -> Label  
  {
     if (!isParis[label]) return LabelMax;
    else if (isBack [label]) return ParisIndex[label];
    else if (isFront[label]) return ParisIndex[label] - Colib::lutEntries(isBack);
    else return LabelMax;
  });

  static constexpr auto ParisRingNumber = Colib::LUT<LUT_s> ([](Label const label) -> Label  {
    if (isFR1[label]) return 1;
    else if (isBR2[label] || isFR2[label]) return 2;
    else if (isBR3[label] || isFR3[label]) return 3;
    else return 0;
  });

  /////////////////////
  // Detectors types //
  /////////////////////

  namespace Detector
  {
    enum Types {Ge, BGO, Ref, RF, Paris, Dssd, DEFAULT};
    constexpr size_t nbTypes = DEFAULT;
    constexpr std::array<std::string, nbTypes> Names = {"Ge", "BGO", "Ref", "RF", "Paris", "Dssd"};
  }

  static constexpr auto detectorType = Colib::LUT<NSI136::LUT_s> ([](Label label){
    if (isGe[label]) return Detector::Ge;
    else if (isBGO[label]) return Detector::BGO;
    else if (isRef[label]) return Detector::Ref;
    else if (isRF[label]) return Detector::RF;
    else if (isParis[label]) return Detector::Paris;
    else if (isDssd[label]) return Detector::Dssd;
    else return Detector::DEFAULT;
  });

  static constexpr auto detectorTypeName = Colib::LUT<NSI136::LUT_s> ([](Label label) -> std::string { 
    auto const & type = detectorType[label];
    if (type == Detector::DEFAULT) return "NaD"; // return NaD (Not a Detector)
    else return Detector::Names[type];
  });

  /////////////////////
  // Detectors names //
  /////////////////////

  static constexpr std::array<std::string, 4> GeLeafs{"red", "green", "black", "blue"};
  static const auto CloverRingNumberStr = Colib::LUT<LUT_s> ([](Label label) -> std::string {return "R"+std::to_string(CloverRingNumber[label]);});
  static const auto CloverRingIndexStr  = Colib::LUT<LUT_s> ([](Label label) -> std::string {return "A"+std::to_string(CloverRingIndex [label]);});
  static const auto CloverName = Colib::LUT<LUT_s> ([](Label label) -> std::string {return CloverRingNumberStr[label]+CloverRingIndexStr[label];});
  static const auto GeName     = Colib::LUT<LUT_s> ([](Label label) -> std::string {return (isGe[label]) ? CloverName[label] + "_" + GeLeafs[crystalIndex[label]-2] : "NaD";});
  static const auto BGOName    = Colib::LUT<LUT_s> ([](Label label) -> std::string {return CloverName[label] + "_BGO" + std::to_string(crystalIndex[label]);});
  
  static const auto ParisWallStr = Colib::LUT<LUT_s> ([](Label label) -> std::string {
    if (isFront[label]) return "F"; 
    else if (isBack[label]) return "B"; 
    else return "ERROR";
  });
  static const auto ParisRingIndexName  = Colib::LUT<LUT_s> ([](Label label) -> std::string {return "R"+std::to_string(ParisRingIndex [label]);});
  static const auto ParisRingNumberName = Colib::LUT<LUT_s> ([](Label label) -> std::string {return "D"+std::to_string(ParisRingNumber[label]);});
  static const auto ParisName           = Colib::LUT<LUT_s> ([](Label label) -> std::string {return "Paris_" + ParisWallStr[label] + ParisRingNumberName[label] + ParisRingIndexName[label];});
  
  static const auto DssdStripName = Colib::LUT<LUT_s> ([](Label label) -> std::string {
    if (isSector1[label]) return "S1"; 
    else if (isSector2[label]) return "S2"; 
    else if (isRing[label]) return "R"; 
    else return "ERROR";
  });
  static const auto DssdName = Colib::LUT<LUT_s> ([](Label label) -> std::string {return "Dssd_" + DssdStripName[label] + "_" + std::to_string(DssdIndex[label]);});

  static const auto detectorName = Colib::LUT<NSI136::LUT_s> ([](Label label) -> std::string { 
    switch (detectorType[label])
    {
      case Detector::Ge    : return GeName [label];
      case Detector::BGO   : return BGOName[label];
      case Detector::Ref   : return "R1A9_Ref_LaBr3";
      case Detector::RF    : return "RF";
      case Detector::Paris : return ParisName[label];
      case Detector::Dssd  : return DssdName [label];
      default: return "NaD"; // return NaD (Not a Detector)
    }
  });

  static constexpr Label maxLabel = Label_cast(LUT_s);

  static constexpr ADC GeOverflow = 5e5; // The overflow bin is almost the same for all Ge
  static constexpr ADC GeADCThreshold = 500;
}