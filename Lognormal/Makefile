GPP = mpicxx   #parallel compilation
#GPP = g++       #single processor

#FT directory
FTDIR = /home/montefra/data1/Cpp/FT
#contains functions often used functions
CFDIR =  /home/montefra/data1/Cpp/Commonfunctions
INCLUDE = -I$(FTDIR) -I$(CFDIR) -I/home/montefra/include 
#-I/usr/lib/openmpi/include
LINK = -L/home/montefra/lib 
#-I/usr/lib/openmpi/lib

fftw = -lfftw3_mpi -lfftw3 -lm #fftw options
gsl = -lgsl -lgslcblas -lm #gsl options
#mpi = -lmpi #mpi options
options = -Wall -O2 #-static #std options
debugmemory = -fmudflap -lmudflap #-L/usr/lib  #debugging options
debugoptions = -Wall -g #$(debugmemory) #-static #debug options

extraobjects = $(FTDIR)/1Dft_corrfunc_pk.o $(CFDIR)/io.o $(CFDIR)/mystrings.o $(CFDIR)/vector_array.o $(FTDIR)/fftw_r2c_c2r.o 
extrainclude = $(FTDIR)/1Dft_corrfunc_pk.h $(CFDIR)/io.h $(CFDIR)/mystrings.h $(CFDIR)/vector_array.h $(FTDIR)/fftw_r2c_c2r.h 

all: lognormal

lognormal: lognormal.o lognormal_functions.o $(extraobjects)
	$(GPP) -o lognormal_vpoisson lognormal.o lognormal_functions.o $(extraobjects) $(LINK) $(mpi) $(fftw) $(gsl) $(debugoptions)

lognormal.o: lognormal.cpp lognormal.h $(extrainclude)
	$(GPP) -c lognormal.cpp $(INCLUDE) $(debugoptions)

lognormal_functions.o: lognormal_functions.cpp lognormal.h $(CFDIR)/vector_array.h
	$(GPP) -c lognormal_functions.cpp $(INCLUDE) $(debugoptions)

clean:
	rm *.o 

cleanall:
	rm *.o lognormal
