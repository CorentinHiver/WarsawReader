#!/bin/bash

# cd LibCo

# base_url="https://raw.githubusercontent.com/CorentinHiver/Nuball2/tree/master/lib/"
# files=(errors.hpp files_functions.hpp libCo.hpp libRoot.hpp print.hpp randomCo.hpp string_functions.hpp vector_functions.hpp)

# for file in "${files[@]}"; do
#     wget -N $file "$base_url/$file"
# done

# cd ..


cd LibCo

base_url="https://raw.githubusercontent.com/CorentinHiver/Nuball2/master/lib"
files=(
    errors.hpp
    files_functions.hpp
    libCo.hpp
    libRootHeader.hpp
    libRoot.hpp
    print.hpp
    randomCo.hpp
    string_functions.hpp
    vector_functions.hpp
    Classes/FilesManager.hpp
    Classes/Timer.hpp
    MTObjects/MTObject.hpp
)

for file in "${files[@]}"; do
    wget -N "$base_url/$file"
done

cd ..