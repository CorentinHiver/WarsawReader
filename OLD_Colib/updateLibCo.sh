#!/bin/bash

# use it in root folder like $ ./Colib/lib/updateColib.sh

cd Colib

# base_url="https://raw.githubusercontent.com/CorentinHiver/Colib/lib/master/lib"
base_url="https://gitlab.com/username/Colib/lib/-/raw/main/lib"
files=(
    errors.hpp
    files_functions.hpp
    Colib.hpp
    libRootHeader.hpp
    libRoot.hpp
    print.hpp
    randomCo.hpp
    string_functions.hpp
    vector_functions.hpp
    Classes/Timer.hpp
    Classes/Timeshifts.hpp
    MTObjects/MTObject.hpp
)

for file in "${files[@]}"; do
    mkdir -p "$(dirname "$file")"
    wget -N -O "$file" "$base_url/$file"
done

cd ..