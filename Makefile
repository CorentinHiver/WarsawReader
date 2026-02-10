STDLIB=# -std=c++17
LIBS=
INCLUDES=
# OPT= -O2
OPT= -O3
# OPT= -g

all: writeTraces studyCFD caen2root rootReaderExample

writeTraces: writeTraces.cpp $(LIBS)
	g++ -o writeTraces writeTraces.cpp $(LIBS) -Wall -Wextra $(INCLUDES) `root-config --cflags` `root-config --glibs` $(OPT) $(STDLIB)

studyCFD: studyCFD.cpp $(LIBS)
	g++ -o studyCFD studyCFD.cpp $(LIBS) -Wall -Wextra $(INCLUDES) `root-config --cflags` `root-config --glibs` $(OPT) $(STDLIB)

caen2root: caen2root.cpp $(LIBS)
	g++ -o caen2root caen2root.cpp $(LIBS) -Wall -Wextra $(INCLUDES) `root-config --cflags` `root-config --glibs` $(OPT) $(STDLIB)

rootReaderExample: rootReaderExample.cpp $(LIBS)
	g++ -o rootReaderExample rootReaderExample.cpp $(LIBS) $(INCLUDES) `root-config --cflags` `root-config --glibs` $(OPT) $(STDLIB)

clean:
	rm -fr writeTraces
	rm -fr studyCFD
	rm -fr caen2root
	rm -fr rootReaderExample