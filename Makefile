STDLIB=
LIBS=AnalysisLib/* Caenlib/* LibCo/* Triggers/*

all: writeTraces studyCFD caen2root

writeTraces: writeTraces.cpp $(LIBS)
	g++ -o writeTraces writeTraces.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

studyCFD: studyCFD.cpp $(LIBS)
	g++ -o studyCFD studyCFD.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

caen2root: caen2root.cpp $(LIBS)
	g++ -o caen2root caen2root.cpp -Wall -Wextra `root-config --cflags` `root-config --glibs` -O2 $(STDLIB)

cpp17:
	@echo "Compiling with C++17 standard"
	$(MAKE) STDLIB="-std=c++17" all
