STDLIB= -std=c++17
GCC_OPTIONS= -Wall -Wextra
OPTIONS=
OPT= -O3 -march=native
TRIGGER=

all: clean
all: writeTraces studyCFD caen2root rootReaderExample

writeTraces: writeTraces.cpp 
	g++ -o writeTraces writeTraces.cpp $(GCC_OPTIONS) `root-config --cflags` `root-config --glibs` $(OPT) $(STDLIB) $(TRIGGER) $(OPTIONS)

studyCFD: studyCFD.cpp 
	g++ -o studyCFD studyCFD.cpp $(GCC_OPTIONS) `root-config --cflags` `root-config --glibs` $(OPT) $(STDLIB) $(TRIGGER) $(OPTIONS)

caen2root: caen2root.cpp 
	g++ -o caen2root caen2root.cpp $(GCC_OPTIONS) `root-config --cflags` `root-config --glibs` $(OPT) $(STDLIB) $(TRIGGER) $(OPTIONS)

rootReaderExample: rootReaderExample.cpp 
	g++ -o rootReaderExample rootReaderExample.cpp $(GCC_OPTIONS) `root-config --cflags` `root-config --glibs` $(OPT) $(STDLIB) $(TRIGGER) $(OPTIONS)

debug: OPT= -g
debug: clean
debug: all

clean:
	rm -fr writeTraces
	rm -fr studyCFD
	rm -fr caen2root
	rm -fr rootReaderExample