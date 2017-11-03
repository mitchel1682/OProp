clear 
#CFLAGS='-lrt -Wall -pedantic -lpthread'
rm ./bin/* -f 
 
g++ -g -c nrlmsise/nrlmsise-00.c -o ./bin/nrlmsise-00.o #$CFLAGS
g++ -g -c nrlmsise/nrlmsise-00_data.c -o ./bin/nrlmsise-00_data.o #$CFLAGS
g++ -g propagator.cpp ./bin/nrlmsise-00.o ./bin/nrlmsise-00_data.o -o ./bin/main -larmadillo -llapack -lblas  #$CFLAGS