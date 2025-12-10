#!/bin/bash

# use it in root folder like $ ./LibCo/updateLibCo.sh

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
    Classes/Timeshifts.hpp
    MTObjects/MTObject.hpp
)

for file in "${files[@]}"; do
    mkdir -p "$(dirname "$file")"
    wget -N -O "$file" "$base_url/$file"
done

cd ..