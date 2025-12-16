STDLIB=-std=c++17
LIBS=
INCLUDES=
# LIBS=$(wildcard AnalysisLib/*.h* Caenlib/*.h* LibCo/*.h* Triggers/*.h*)
# INCLUDES=-IAnalysisLib -ICaenlib -ILibCo -ITriggers

all: writeTraces studyCFD caen2root rootReaderExample

writeTraces: writeTraces.cpp $(LIBS)
	g++ -o writeTraces writeTraces.cpp $(LIBS) -Wall -Wextra $(INCLUDES) `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

studyCFD: studyCFD.cpp $(LIBS)
	g++ -o studyCFD studyCFD.cpp $(LIBS) -Wall -Wextra $(INCLUDES) `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

caen2root: caen2root.cpp $(LIBS)
	g++ -o caen2root caen2root.cpp $(LIBS) -Wall -Wextra $(INCLUDES) `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

rootReaderExample: rootReaderExample.cpp $(LIBS)
	g++ -o rootReaderExample rootReaderExample.cpp $(LIBS) $(INCLUDES) `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

clean:
	rm -fr writeTraces
	rm -fr studyCFD
	rm -fr caen2root
	rm -fr rootReaderExample