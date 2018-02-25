#!/bin/bash
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ROOTSYS/lib
for filename in ../PCORAW/data_180215/*.pcoraw
do
	./PCOtoROOT "$filename" "../ROOT/$(basename "$filename" .pcoraw).root" 
done
