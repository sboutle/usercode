LD            = g++


all: runAnalysis

runAnalysis: analysis.o ./eventVars_h.so
	$(LD) analysis.o ./eventVars_h.so -o runAnalysis `root-config --libs --cflags` 

analysis.o: analysis.cc analysis.hh eventVars.h 
	$(LD) -I `root-config --incdir` -c analysis.cc
#	$(LD) -I `root-config --incdir` -c analysis.cc  `root-config --libs --cflags`

#eventVars_h.so: analysis.cc analysis.hh eventVars.h 
#	$(LD) -I `root-config --incdir` -c analysis.cc

clean:
	rm -rf *.o runit *Dict.* *.def 
