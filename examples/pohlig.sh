#!/bin/bash
HOME=Users/poojagarg/Downloads/c_library/givaro-3.8.0
g++ -g -DHAVE_CONFIG_H -I. -I/$HOME  -I/$HOME   -I/$HOME/src/kernel/system -I/$HOME/src/kernel/memory -I/$HOME/src/kernel/zpz -I/$HOME/src/kernel/integer -I/$HOME/src/kernel -I/$HOME/src/library/poly1 -I/$HOME/src/kernel/bstruct -I/$HOME/src/library/tools  -O2   -Wall -g -DNDEBUG -UGIVARO_DEBUG -UDEBUG -std=gnu++11 -c -o pohlig.o pohlig.cpp 
/bin/sh /$HOME/libtool  --tag=CXX   --mode=link g++ -O2   -Wall -g -DNDEBUG -UGIVARO_DEBUG -UDEBUG -std=gnu++11 -static  -o pohlig pohlig.o ../src/ellipticCurve.o ../src/attacksECC.o -L/$HOME/src -L/$HOME/src/.libs -lgivaro -lgmpxx -lgmp  
