#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include <unordered_map>
#include <typeindex>

/**
 * @brief Creates a unordered_map of all the object of a certain class (TH1F, TH2F...) inside a TFile, indexed by their name
 * @details
 * TFile* file(TFile::Open("file.root","read"));
 * auto list = file_get_map_of<TH1F>(file);
 * 
 * If no file is passed as parameter, reads the current file.
 * Internally perform a file->cd()
 */
template<class T>
std::unordered_map<std::string, T*> file_get_map_of(TFile* file = nullptr)
{
  // init
  std::unordered_map<std::string, T*> ret;
  T temp_obj; 

  // Check the files :
  if (file == nullptr) file = gFile;
  if (!file) { error("in file_get_map_of<", temp_obj.ClassName(), ">(TFile* file): file is nullptr"); return ret;}
  file->cd();

  // Get the class name :
  auto const & classname = temp_obj.ClassName();
  
  // Loop over the list of keys of every object in the TFile :
  auto list = file->GetListOfKeys();
  for (auto&& keyAsObj : *list)
  {
    auto key = dynamic_cast<TKey*>(keyAsObj);
    if(strcmp(key->GetClassName(), classname) == 0) 
    {
      T* obj = dynamic_cast<T*>(key->ReadObj());
      if (obj)
      {
        ret.emplace(obj->GetName(), obj);
      }
    }
  }

  return ret;
}

// static const std::unordered_map<std::type_index, std::string> typeRootMap = 
// {
//   // Bool :
//   {static_cast<std::type_index>(typeid(          true)), "O"},

//   // Integers :
//   {static_cast<std::type_index>(typeid(  static_cast<char>(1))), "B"}, {static_cast<std::type_index>(typeid( static_cast<uchar>(1))), "b"},
//   {static_cast<std::type_index>(typeid( static_cast<short>(1))), "S"}, {static_cast<std::type_index>(typeid(static_cast<ushort>(1))), "s"},
//   {static_cast<std::type_index>(typeid(   static_cast<int>(1))), "I"}, {static_cast<std::type_index>(typeid(  static_cast<uint>(1))), "i"},
//   {static_cast<std::type_index>(typeid(  static_cast<long>(1))), "G"}, {static_cast<std::type_index>(typeid( static_cast<ulong>(1))), "g"},

//   // Floating point :
//   {static_cast<std::type_index>(typeid(static_cast<double>(1))), "D"}, {static_cast<std::type_index>(typeid( static_cast<float>(1))), "F"},

//   // ROOT types :
//   {static_cast<std::type_index>(typeid(static_cast<Long64_t>(1))), "L"}, {static_cast<std::type_index>(typeid(static_cast<ULong64_t>(1))), "l"}
// };

// template<class T>
// std::string typeRoot(T const & t)
// {
//   std::type_index typeIndex{typeid(t)};
//   auto it = typeRootMap.find(typeIndex);
//   if (it != typeRootMap.end()) return it->second;
//   else                         return "Unknown";
// }

// template<class T>
// std::string typeRoot()
// {
//   T t;
//   std::type_index typeIndex{typeid(t)};
//   auto it = typeRootMap.find(typeIndex);
//   return (it != typeRootMap.end()) ? it->second : "Unknown";
// }

// /// @brief Create a branch for a given array and name
// /// @param name_size: The name of the leaf that holds the size of the array
// template<class T>
// TBranch* createBranchArray(TTree* tree, std::string const & name, T * array, std::string const & name_size, int buffsize = 32000)
// {
//   auto type_root_format = name+"["+name_size+"]/"+typeRoot(**array);
//   return (tree -> Branch(name.c_str(), array, type_root_format, buffsize));
// }
#include <string>
#include <typeinfo>

template<class T> std::string typeRoot()          {return "Unknown";}
template<> std::string typeRoot<bool>()           {return "O";}
template<> std::string typeRoot<char>()           {return "B";}
template<> std::string typeRoot<unsigned char>()  {return "b";}
template<> std::string typeRoot<short>()          {return "S";}
template<> std::string typeRoot<unsigned short>() {return "s";}
template<> std::string typeRoot<int>()            {return "I";}
template<> std::string typeRoot<unsigned int>()   {return "i";}
template<> std::string typeRoot<long>()           {return "G";}
template<> std::string typeRoot<unsigned long>()  {return "g";}
template<> std::string typeRoot<double>()         {return "D";}
template<> std::string typeRoot<float>()          {return "F";}
template<> std::string typeRoot<Long64_t>()       {return "L";}
template<> std::string typeRoot<ULong64_t>()      {return "l";}

/// @brief Create a branch for a given array and name
/// @param name_size: The name of the leaf that holds the size of the array
template<class T>
TBranch* createBranchArray(TTree* tree, std::string const & name, T * array, std::string const & name_size, int buffsize = 32000)
{
  std::string type_root_format = name + "[" + name_size + "]/" + typeRoot<std::remove_extent_t<T>>();
  return tree->Branch(name.c_str(), array, type_root_format.c_str(), buffsize);
}