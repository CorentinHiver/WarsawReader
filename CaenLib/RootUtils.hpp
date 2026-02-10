#include "TFile.h"
#include "TKey.h"
#include "TTree.h"

#include <string>
#include <typeindex>
#include <typeinfo>
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
template <class T>
std::unordered_map<std::string, T*> file_get_map_of(TFile* file = nullptr)
{
    static_assert(std::is_base_of<TObject, T>::value,
                  "file_get_map_of<T>: T must inherit from TObject");

    std::unordered_map<std::string, T*> ret;

    // Resolve file
    if (!file) file = gFile;
    if (!file || file->IsZombie()) {
        error("file_get_map_of", "Invalid or null TFile");
        return ret;
    }

    file->cd();

    // Use ROOT RTTI instead of constructing a temporary object
    TClass* cls = TClass::GetClass(typeid(T));
    if (!cls) {
        error("file_get_map_of", "Could not resolve TClass for requested type");
        return ret;
    }

    const char* classname = cls->GetName();

    // Iterate over keys
    TIter next(file->GetListOfKeys());
    while (auto* key = static_cast<TKey*>(next())) {
        if (std::strcmp(key->GetClassName(), classname) != 0)
            continue;

        auto* obj = dynamic_cast<T*>(key->ReadObj());
        if (!obj)
            continue;

        ret.emplace(obj->GetName(), obj);
    }

    return ret;
}

template <typename T>
struct RootLeafType {
    static constexpr const char* value = nullptr;
};

template <> struct RootLeafType<bool>           { static constexpr const char* value = "O"; };
template <> struct RootLeafType<char>           { static constexpr const char* value = "B"; };
template <> struct RootLeafType<unsigned char>  { static constexpr const char* value = "b"; };
template <> struct RootLeafType<short>          { static constexpr const char* value = "S"; };
template <> struct RootLeafType<unsigned short> { static constexpr const char* value = "s"; };
template <> struct RootLeafType<int>            { static constexpr const char* value = "I"; };
template <> struct RootLeafType<unsigned int>   { static constexpr const char* value = "i"; };
template <> struct RootLeafType<long>           { static constexpr const char* value = "G"; };
template <> struct RootLeafType<unsigned long>  { static constexpr const char* value = "g"; };
template <> struct RootLeafType<float>          { static constexpr const char* value = "F"; };
template <> struct RootLeafType<double>         { static constexpr const char* value = "D"; };
template <> struct RootLeafType<Long64_t>       { static constexpr const char* value = "L"; };
template <> struct RootLeafType<ULong64_t>      { static constexpr const char* value = "l"; };

template <typename T>
constexpr const char* typeRoot()
{
    static_assert(RootLeafType<T>::value != nullptr, "typeRoot<T>: unsupported type for ROOT branch");
    return RootLeafType<T>::value;
}

//------------------------------------------------------------------------------
// Create array branch
//------------------------------------------------------------------------------

template <class T>
TBranch* createBranchArray(TTree* tree, const std::string& name, T* array, const std::string& name_size, int buffsize = 32000)
{
    if (!tree) {error("createBranchArray", "TTree is null"); return nullptr;}
    using element_t = std::remove_cv_t<std::remove_extent_t<T>>;
    const std::string leaflist = name + "[" + name_size + "]/" + typeRoot<element_t>();
    return tree->Branch(name.c_str(), array, leaflist.c_str(), buffsize);
}