#!/bin/csh
#
set echo
#
cp subset.H /$HOME/include
#
g++ -c -g subset.C >& compiler.out
if ( $status != 0 ) then
  echo "Errors compiling subset.C."
  exit
endif
rm compiler.out
#
mv subset.o ~/libcpp/$ARCH/subset.o
#
echo "A new version of subset has been created."
