#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include <unordered_map>

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

static const std::unordered_map<std::type_index, std::string> typeRootMap = 
{
  // Bool :
  {static_cast<std::type_index>(typeid(          true)), "O"},

  // Integers :
  {static_cast<std::type_index>(typeid(  char_cast(1))), "B"}, {static_cast<std::type_index>(typeid( uchar_cast(1))), "b"},
  {static_cast<std::type_index>(typeid( short_cast(1))), "S"}, {static_cast<std::type_index>(typeid(ushort_cast(1))), "s"},
  {static_cast<std::type_index>(typeid(   int_cast(1))), "I"}, {static_cast<std::type_index>(typeid(  uint_cast(1))), "i"},
  {static_cast<std::type_index>(typeid(  long_cast(1))), "G"}, {static_cast<std::type_index>(typeid( ulong_cast(1))), "g"},

  // Floating point :
  {static_cast<std::type_index>(typeid(double_cast(1))), "D"}, {static_cast<std::type_index>(typeid( float_cast(1))), "F"},

  // ROOT types :
  {static_cast<std::type_index>(typeid(static_cast<Long64_t>(1))), "L"}, {static_cast<std::type_index>(typeid(static_cast<ULong64_t>(1))), "l"}
};

template<class T>
std::string typeRoot(T const & t)
{
  auto const & typeIndex = static_cast<std::type_index>(typeid(t));
  auto it = typeRootMap.find(typeIndex);
  if (it != typeRootMap.end()) return it->second;
  else                         return "Unknown";
}

template<class T>
std::string typeRoot()
{
  T t;
  auto const & typeIndex = static_cast<std::type_index>(typeid(t));
  auto it = typeRootMap.find(typeIndex);
  if (it != typeRootMap.end()) return it->second;
  else                         return "Unknown";
}

/// @brief Create a branch for a given array and name
/// @param name_size: The name of the leaf that holds the size of the array
template<class T>
auto createBranchArray(TTree* tree, std::string const & name, T * array, std::string const & name_size, int buffsize = 64000)
{
  auto const & type_root_format = name+"["+name_size+"]/"+typeRoot(**array);
  return (tree -> Branch(name.c_str(), array, type_root_format.c_str(), buffsize));
}