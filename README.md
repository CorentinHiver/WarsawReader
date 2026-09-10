# WarsawReader
corentin.hiver@free.fr

NEW: need to install the Colib library
cd && git clone https://gitlab.com/CorentinHiver/Colib.git && cd Colib && sh Install.sh

Compilation:

To compile with make:
> make clean
> make -j4

Three executables are compiled by make:

caen2root   Main purpose of WarsawReader, converts data in .caendat to .root format, with possible event building, time synchronization, and various other options
writeTraces Reads the traces of the data in .caendat and writes in a .root file
studyCFD    After event-building, checks the time correlations between detectors

macros:
macros/Miscellaneous.C    Collection of functions to analysis the tree created by caen2root (examples, timeshifts, spectra, special requests ...)
macros/dTcalculate.C      Calculates the time resolution of hits written in a file with one detector used as time reference


make options:
TRIGGER="-DTRIGGER=[TriggerName]" (Read Triggers/)

Update using:
> sh UPDATE.sh

- General purpose tools : 

Tools/caendatInspect        Dumps raw data from .caendat datafile in the terminal
Tools/caendatIntegrityCheck Checks if the raw .caendat datafile is not corrupted

- Special tools :

writeTraces    Writes the traces in TGraphs with CFD
