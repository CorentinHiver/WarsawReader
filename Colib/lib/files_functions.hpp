#pragma once
#include "print.hpp"
#include "vector_functions.hpp"
#include "string_functions.hpp"

#include <algorithm>
#include <cmath>
#include <dirent.h>
#include <fstream>
#include <glob.h>
#include <iostream>
#include <map>
#include <unordered_map>
#include <tuple>

#if __GNUC__ >= 9 // g++ version
	#include <filesystem>
  #if __cplusplus >= 201703L
    #include <filesystem>
    namespace fs = std::filesystem;
  #else 
    #warning ("In version of c++ < 17, '...' fold expression not defined, and <filesystem> not loaded. Some parts of the code might not behave as expected...")
  #endif // C++ 17
#else // __GNUC__ < 9 
	#include <experimental/filesystem>
  #if __cplusplus >= 201703L
	  namespace fs = std::experimental::filesystem;
  #else 
    #warning ("In version of c++ < 17, '...' fold expression not defined, and <filesystem> not loaded. Some parts of the code might not behave as expected...")
  #endif // C++ 17
#endif
//----------------------------------------------------//
//       General files and folders manipulations      //
//----------------------------------------------------//
namespace Colib 
{
  /// @brief "/path/to/file.ext" -> "file.ext"
  std::string removePath(std::string const & file) 
  {
    size_t pos = file.find_last_of("/");
    return (pos == std::string::npos) ? file : file.substr(pos + 1);
  }

  /// @brief "/path/to/file.ext" -> "/path/to/"
  std::string getPath(std::string const & file) 
  { 
    size_t pos = file.find_last_of("/");
    return (pos == std::string::npos) ? "" : file.substr(0, pos + 1);
  }

  /// @brief "/path/to/file.ext" -> "/path/to/file"
  std::string removeExtension(std::string const & file) 
  {
    std::string name = removePath(file);
    size_t pos = name.find_last_of(".");
    // If no dot, or dot is the first character (hidden file)
    if (pos == std::string::npos || pos == 0) return file; 
    return getPath(file) + name.substr(0, pos);
  }

  /// @brief "/path/to/file.ext" -> ".ext"
  std::string extension(std::string const & file) 
  { 
    std::string name = removePath(file);
    size_t pos = name.find_last_of(".");
    if (pos == std::string::npos || pos == 0) return "";
    return name.substr(pos + 1);
  }


  /// @brief "/path/to/file.ext" -> "file"
  std::string getShortname(std::string const & file) 
  { 
    std::string name = removePath(file);
    size_t pos = name.find_last_of(".");
    if (pos == std::string::npos || pos == 0) return name;
    return name.substr(0, pos);
  }

  /// @brief Return /home/usr
  std::string getHome() {return std::getenv("HOME");}
  /// @brief Return /path/to/current/folder
  std::string getPwd () {return std::getenv("PWD" );}
  /// @brief Return /home/usr/
  std::string getHomePath() {return std::getenv("HOME")+std::string("/");}
  /// @brief Return /path/to/current/folder/
  std::string getPwdPath () {return std::getenv("PWD" )+std::string("/");}

  /// @brief "/path/to/file.ext" -> "/path/to/file.new_ext"
  std::string setExtension(std::string const & string, std::string const & extension)
  {
    return removeExtension(string)+"."+extension;
  }

  /// @brief appendFilename("/path/to/file.ext", "app") -> "/path/to/file_app.ext"
  std::string appendFilename(std::string const & string, std::string const & app) 
  {
    return removeExtension(string)+app+"."+extension(string);
  }  

  std::string getFullPath(std::string const & file)
  {
    if (file.empty()) return "";

    std::string path = getPath(file);

    // Path anchoring
    if (!path.empty()) 
    {
           if (path[0] == '~') path = getHome()    + "/" + path.substr(1);
      else if (path[0] != '/') path = getPwdPath() + "/" + path;
    } 
    else path = getPwdPath() + "/";

    // Path resolution
    // Pass removeVoids = true to clean double slashes ("//") and prevent empty elements
    auto folders = getList(path, "/", true); 
    std::vector<std::string> output;
    
    for (auto const & folder : folders)
    {
      if (folder == ".") continue;
      else if (folder == "..") {
        if (!output.empty()) output.pop_back();}
      else output.push_back(folder);
    }

    std::string ret = "/"; // Guaranteed absolute Unix path starts with root
    for (auto const & folder : output) ret += folder + "/";
    
    return ret;
  }

  std::string nicerPath(std::string const & path)
  {
    return getFullPath(path) + removePath(path);
  }
  
  bool fileIsEmpty(std::ifstream& file) { return file.peek() == std::ifstream::traits_type::eof();}
  
  void goToBeginning(std::ifstream& file) {file.seekg(0, std::ios::beg);}
  
  static const std::unordered_map<std::string, float> size_file_unit =
  {
    {"o",  1.},
    {"ko", 1024.},
    {"Mo", 1048576.},
    {"Go", 1073741824.},
    {"B",  1.},
    {"kB", 1024.},
    {"MB", 1048576.},
    {"GB", 1073741824.}
  };
  
  float sizeFileConversion(float const & size, std::string const & unit_i, std::string const & unit_o)
  {
    return size * (size_file_unit.at(unit_i)/size_file_unit.at(unit_o));
  }
  
  float sizeFile(std::ifstream& file, std::string const & unit = "o")
  {
    auto const init = file.tellg();
    file.seekg(0, std::ios::end);
    auto const ret = file.tellg();
    file.seekg(init);// Go back to inital place in the file
    return (static_cast<float>(ret)/size_file_unit.at(unit));
  }
  
  float sizeFile(std::string filename, std::string const & unit = "o")
  {
    std::ifstream f (filename, std::ios::binary);
    return sizeFile(f, unit);
  }
  
  std::tuple<double, std::string> sizeFileBestUnit(std::ifstream& file)
  {
    auto const & size = sizeFile(file);
    auto order_of_magnitude = int(std::log10(std::abs(size)));
    std::string unit;
         if (order_of_magnitude<3) unit = "o" ;
    else if (order_of_magnitude<6) unit = "ko";
    else if (order_of_magnitude<9) unit = "Mo";
    else                           unit = "Go";
    return std::tuple<double, std::string>(sizeFileConversion(size, "o", unit), unit); 
  }

  std::tuple<double, std::string> sizeFileBestUnit(std::string filename)
  {
    std::ifstream f (filename, std::ios::binary);
    return sizeFileBestUnit(f);
  }
  
  std::string sizeFileBestUnitString(std::ifstream& file)
  {
    auto const & size = sizeFile(file);
    auto order_of_magnitude = int(std::log10(std::abs(size)));
    std::string unit;
         if (order_of_magnitude<3) unit = "o" ;
    else if (order_of_magnitude<6) unit = "ko";
    else if (order_of_magnitude<9) unit = "Mo";
    else                           unit = "Go";
    return std::to_string(sizeFileConversion(size, "o", unit)) + " " + unit; 
  }

  std::string sizeFileBestUnitString(std::string filename)
  {
    std::ifstream f (filename, std::ios::binary);
    return sizeFileBestUnitString(f);
  }
    
  bool fileExists(std::string fileName)
  {
    std::string path = getFullPath(fileName);
    std::string name = removePath(fileName);
    print(path, name);
    struct dirent *file = nullptr;
    DIR *dp = nullptr;
    dp = opendir(path.c_str());
    if(dp == nullptr) return false;
    else
    {
      while ((file = readdir(dp)))
      {
        if (!strcmp(file -> d_name, name.c_str()))
        {
          closedir(dp);
          return true;
        }
      }
    }
    closedir(dp);
    return false;
  }
  
  std::string & appendIfMissing(std::string & string, char const & character)
  {
    if (string.back() != character) string.push_back(character);
    return string;
  }
  
  bool pathExists(std::string path)
  {
    appendIfMissing(path, '/');
    DIR *dp = nullptr;
    dp = opendir(path.c_str());
    bool ret = (dp != nullptr);
    std::string str = ((dp!=nullptr) ? "oui" : "non");
    if (dp) closedir(dp);
    return ret;
  }
  
  bool pathExists(std::string path, bool const & verbose)
  {
    appendIfMissing(path, '/');
    if (pathExists(path)) return true;
    if (verbose) print("Path", path, "not found...");
    return false;
  }
  
  void mkdir(std::string path, bool verbose = 0)
  {
    if (path=="")
    {
      print("Colib::mkdir(): No path !");
      return;
    }
#ifdef CoMT
    MT::lock_mutex lock(Colib::MT::mutex);
#endif //CoMT
    path = getFullPath(path);
    if (!pathExists(path))
    {
      if (verbose) print("Creating path", path);
      // mkdir -p to create the full path if needed (otherwise crashes if some directory of the path is missing)
      system(("mkdir -p "+path).c_str());
    }
  }

  std::string fileEnsurePath(std::string const & fileName)
  {
    std::string ret = getFullPath(fileName);
    mkdir(ret);
    return ret;
  }

  void rmFile(std::string const & filename)
  {
    std::remove(filename.c_str());
  }
  
  int nbFilesInFolder(std::string & folder)
  {
    int ret = -1;
    appendIfMissing(folder, '/');
    DIR *dp = nullptr;
    dp = opendir(folder.c_str());
    if(dp == nullptr) ret = -1;
    else
    {
      int i = 0;
      while ((readdir(dp))) i++;
      ret = i;
    }
    closedir(dp);
    return ret;
  }
  
  std::vector<std::string> listFilesInFolder(
      std::string const & folder = "./", 
      std::vector<std::string> const & extensions = {"*"},
      bool const & reorder = true
    )
  {
    std::vector<std::string> ret;
    struct dirent *file = nullptr;
    DIR *dp = nullptr;
    dp = opendir(folder.c_str());
    if (dp == nullptr) {print("Folder ", folder, " not found..."); return ret;}
    std::string name = "";
    while ((file = readdir(dp)))
    {
      name = file->d_name;
      std::string const & ext = extension(name);
      if (extensions[0] == "*" || find(extensions.begin(), extensions.end(), ext) != extensions.end()) ret.push_back(folder+name);
    }
    closedir(dp);
    delete file;
    if (reorder) std::sort(ret.begin(), ret.end());
    return ret;
  }
  
  std::vector<std::string> listFileNamesInFolder(
    std::string const & foldername = "./", 
    std::vector<std::string> const & extensions = {"*"},
    bool const & reorder = true
  )
  {
    std::vector<std::string> ret;
    struct dirent *file = nullptr;
    DIR *dp = nullptr;
    dp = opendir(foldername.c_str());
    if (dp == nullptr) {print("Folder ", foldername, " not found..."); return ret;}
    std::string name = "";
    while ( (file = readdir(dp)))
    {
      name = file->d_name;
      std::string const & ext = extension(name);
      if (extensions[0] == "*" || find(extensions.begin(), extensions.end(), ext) != extensions.end()) ret.push_back(name);
    }
    closedir(dp);
    delete file;
    if (reorder) std::sort(ret.begin(), ret.end());
    return ret;
  }
  
  int checkNewFile(std::string & folderName, std::string & lastFile)
  {
    int ret = -1;
    appendIfMissing(folderName, '/');
    DIR *dp = nullptr;
    dp = opendir(folderName.c_str());
    struct dirent *file = nullptr;
    if(dp == nullptr) ret = -1;
    else
    {
      int i = 0;
      while ((file = readdir(dp)))
      {
        lastFile = file -> d_name;
        i++;
      }
      ret = i;
    }
    closedir(dp);
    return ret;
  }
  
  std::vector<std::string> listFileReader(std::string const & filename)
  {
    std::vector<std::string> list;
  
    std::ifstream listFile(filename,std::ios::in);
    if(!listFile.good())
    {
      print("List file", filename, "not found !");
    }
    else
    {
      std::string line;
      while(getline(listFile,line))
      {
        list.push_back(line);
      }
    }
    return list;
  }
  
  bool hasGlobMeta(const std::string& expression) 
  {
    bool escape = false;
    for (char c : expression) 
    {
      if (escape) 
      {
        escape = false;
        continue;
      }
      if (c == '\\') 
      {
        escape = true;
        continue;
      }
      if (c == '*' || c == '?' || c == '[') return true;
    }
    return false;
  }

  std::vector<std::string> findFilesWildcard(const std::string& expression) 
  {
    if (!hasGlobMeta(expression)) return {expression};

    glob_t result;
    int retcode = glob(expression.c_str(), GLOB_TILDE, NULL, &result);
    if (retcode != 0) return {};

    std::vector<std::string> ret;
    ret.reserve(result.gl_pathc);
    for (size_t i = 0; i < result.gl_pathc; ++i) ret.push_back(result.gl_pathv[i]);
    globfree(&result);
    return ret;
  }
  
  void findFilesWildcard(std::string const & expression, std::vector<std::string> & vec)
  {
    auto const files = findFilesWildcard(expression);
    for (auto const & file : files) vec.push_back(file);
  }
}
