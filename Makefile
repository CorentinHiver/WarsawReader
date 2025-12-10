STDLIB=-std=c++17 -std=c++20
LIBS=AnalysisLib/* Caenlib/* LibCo/* Triggers/*

all: writeTraces studyCFD caen2root examples

writeTraces: writeTraces.cpp $(LIBS)
	g++ -o writeTraces writeTraces.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

studyCFD: studyCFD.cpp $(LIBS)
	g++ -o studyCFD studyCFD.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

caen2root: caen2root.cpp $(LIBS)
	g++ -o caen2root caen2root.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

examples: rootReaderExample.cpp $(LIBS)
	g++ -o rootReaderExample rootReaderExample.cpp `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

lib: 

# cpp17:
# 	@echo "Compiling with C++17 standard"
# 	$(MAKE) STDLIB="-std=c++17" all

# cpp20:
# 	@echo "Compiling with C++20 standard"
# 	$(MAKE) STDLIB="-std=c++20" all
