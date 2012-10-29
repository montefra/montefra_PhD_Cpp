###############################################################
# Main makefile to compile the power spectrum codes           #
#                                                             # 
# Default compilation: make                                   #
# options:                                                    #                 
#   GPP=... to change the compiler                            #
#   nowarning=  to turn off warnings of the intel compiler    #
#     [useful -wd981 -wd869 -wd383 -wd593, default none]      #
#   optimize: optimization flag [default -g]                  #
#   includes: list of directories which contains includes for #
#     fftw3, gsl, tclap and cfitsio                           #
#   libs: list of directories containing the libraries for    #
#     fftw3, gsl and cfitsio                                  #
#                                                             # 
# options for "make install"                                  # 
#   prefix=... to change the installation prefix              #
#   [default $HOME/bin]                                       #
###############################################################

#parallel compilation
#GPP = mpic++
GPP = mpicxx
#single processor
#GPP = g++

prefix = $(HOME)/bin

#options
fftw = -lfftw3_mpi -lfftw3 
gsl = -lgsl -lgslcblas
cfitsio = -lcfitsio
nowarning = 
optimize = -g
compiler_options = -Wall $(optimize) $(nowarning) -DMPICH_IGNORE_CXX_SEEK   #options for the precompile and compiler
options = $(compiler_options) $(gsl) $(fftw) $(gsl) $(cfitsio) -lm 

#subdirectories
classes = $(CURDIR)/Power_spectrum_classes
SUBDIRS = $(CFdir) $(classes)

#includes
includes = $(HOME)/include 
#compiler includes
INCLUDE = $(patsubst %,-I%,$(SUBDIRS) $(includes))

#libraries
libs = $(HOME)/lib 
#link them
LINK =  $(patsubst %,-L%,$(libs))

#commond include file in this directory
commonh = parse_options.h common_functions.h 
#box power spectrum specific files
boxc = $(wildcard box*.cpp) 
boxo = $(patsubst %.cpp,%.o, $(boxc)) 
boxh = $(wildcard box*.h) $(commonh)
#fkp power spectrum specific files
fkpc = $(wildcard fkp*.cpp)
fkpo = $(patsubst %.cpp,%.o, $(fkpc)) 
fkph = $(wildcard fkp*.h) $(commonh)

#classes objects and includes
classc = $(wildcard $(classes)/*.cpp) 
classo = $(patsubst %.cpp,%.o, $(classc)) 
classh = $(wildcard $(classes)/*.h)

#names of the executables
exec = boxpower_spectrum fkppower_spectrum

#by default make all the exec
all: $(exec)

#dependances

boxpower_spectrum.o: boxpower_spectrum.cpp $(boxh) $(classo)
boxpower_spectrum_f.o: boxpower_spectrum_f.cpp $(boxh) $(read_snapshot.o) $(classo)

fkppower_spectrum.o: fkppower_spectrum.cpp $(fkph) $(classo)
fkppower_spectrum_f.o: fkppower_spectrum_f.cpp $(fkph) $(classo)

read_snapshot.o: read_snapshot.c read_snapshot.h

$(classo): $(classc) $(classh) common_functions.h

#create the exe
boxpower_spectrum: $(boxo) $(classo) read_snapshot.o
	$(GPP) -o $@ $(boxo) read_snapshot.o $(classo) $(options) $(LINK) 

fkppower_spectrum: $(fkpo) $(classo) 
	$(GPP) -o $@ $(fkpo) $(classo) $(options) $(LINK)
	
#-Wl,-rpath=$(LIB)

%.o: %.cpp
	$(GPP) -c $< -o $@ $(compiler_options) $(INCLUDE)
%.o: %.c
	$(GPP) -c $< -o $@ $(compiler_options) $(INCLUDE)

install:
	cp $(exec) $(prefix)

help:
	@head -n18 Makefile
.PHONY : clean cleanall
clean:
	rm *.o
cleanall:
	rm *.o $(classo) $(exec)

