#GPP = mpicxx   #parallel compilation
GPP = g++       #single processor

INCLUDE = -I/home/montefra/include
#-I/usr/lib/openmpi/include

fftw = -lfftw3_mpi -lfftw3 -lm  #fftw options
gsl = -lgsl -lgslcblas -lm      #gsl options
options = -Wall -O2 #-static     #std options
debugmemory = -fmudflap -lmudflap #-L/usr/lib  #debugging options
debugoptions = -Wall -g #$(debugmemory) #-static #debug options

all: ft_1D powspec.o 1Dft_corrfunc_pk.o fftw_r2c_c2r.o

ft_1D: ft_1D.o 1Dft_corrfunc_pk.o
	$(GPP) -o $@ ft_1D.o 1Dft_corrfunc_pk.o $(debugoptions) $(gsl) -L$(HOME)/lib

ft_1D.o: ft_1D.cpp 1Dft_corrfunc_pk.h
	$(GPP) -o $@ -c $< $(debugoptions) -I$(HOME)/include

1Dft_corrfunc_pk.o: 1Dft_corrfunc_pk.cpp 1Dft_corrfunc_pk.h
	$(GPP) -c 1Dft_corrfunc_pk.cpp $(debugoptions)


clean:
	rm *.o

