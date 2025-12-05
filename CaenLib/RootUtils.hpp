#include "TFile.h"
#include "TKey.h"
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