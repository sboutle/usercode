LD            = g++


all: runit

runit: test.o 
	$(LD) -g plotMaker.o -o makePlots `root-config --libs --cflags`

test.o: plotMaker.C plotMaker.h
	$(LD) -I `root-config --incdir` -c plotMaker.C  `root-config --libs --cflags`

#MyDict.cxx: test.hh  Linkdef.h
#	rootcint -f $@ -c `root-config --cflags --libs` -p $^
 
#libMyLib.so: MyDict.cxx test.cc
#	g++ -shared -o$@ `root-config --cflags --libs --ldflags` 
#-I$(ROOTSYS)/include $^
#g++ -I `$ROOTSYS/bin/root-config --incdir` -o executable.bin  example.C `$ROOTSYS/bin/root-config --libs`

clean:
	rm -rf *o runit *Dict.* *.def *.so