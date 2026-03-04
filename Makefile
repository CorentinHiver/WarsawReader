STDLIB= #-std=c++17
LIBS=
INCLUDES= 
OPTIONS=-Wall -Wextra
# OPT= -O2
OPT= -O3
# OPT= -g
TRIGGER=NedaVeto

all: clean writeTraces studyCFD caen2root rootReaderExample

writeTraces: writeTraces.cpp $(LIBS)
	g++ -o writeTraces writeTraces.cpp $(LIBS) $(OPTIONS) `root-config --cflags` `root-config --glibs` $(INCLUDES) $(OPT) $(STDLIB)

studyCFD: studyCFD.cpp $(LIBS)
	g++ -o studyCFD studyCFD.cpp $(LIBS) $(OPTIONS) `root-config --cflags` `root-config --glibs` $(INCLUDES) $(OPT) $(STDLIB)

caen2root: caen2root.cpp $(LIBS)
	g++ -o caen2root caen2root.cpp $(LIBS) $(OPTIONS) `root-config --cflags` `root-config --glibs` $(INCLUDES) $(OPT) $(STDLIB)

rootReaderExample: rootReaderExample.cpp $(LIBS)
	g++ -o rootReaderExample rootReaderExample.cpp $(LIBS) $(OPTIONS) `root-config --cflags` `root-config --glibs` $(INCLUDES) $(OPT) $(STDLIB)

clean:
	rm -fr writeTraces
	rm -fr studyCFD
	rm -fr caen2root
	rm -fr rootReaderExample