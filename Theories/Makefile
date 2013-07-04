###############################################################
# Makefile to compile the MCMC code
#                                                              
# Default compilation: make                                   
# options:                                                     
#   GPP: to change the compiler [default: g++]                
#   optimize: optimization flag [default -g]                  
#   includes: directory containing includes for gsl, if not   
#     in standard path or in $HOME/include                                        
#   libs: directory containing library for gsl, if not in     
#     standard path or in $HOME/lib                                                                                 
###############################################################

#GPP = mpicxx
GPP = g++

# window matrix or not (if done with if might take longer)
#gsl options
gsl = -lgsl -lgslcblas -lm
optimize = -g
options = -Wall $(optimize)
debug = -fmudflap -lmudflap #-L/usr/lib  #debugging options

#contains functions often used functions
CFDIR = ../Commonfunctions
INCLUDE = -I$(CFDIR) -I$(HOME)/include
ifdef includes
  INCLUDE += -I$(includes)
endif
LINK = -L$(HOME)/lib
ifdef libs
  LINK += -L$(libs)
endif

CFo = $(CFDIR)/io.o 

exe = p22 #list of executables

all: $(exe)

#montecarlo markov chain
Objects_p22 = $(patsubst %.cpp,%.o, $(wildcard *.cpp)) 

$(exe): $(Objects_p22) $(wildcard *.h) $(CFo)
	$(GPP) $(options) -o $@ $(Objects_p22) $(CFo) $(gsl) $(LINK)

$(Objects_p22): %.o: %.cpp $(wildcard *.h) 
	$(GPP) -c $< -o $@ $(options) $(INCLUDE)

#create the object files for common functions
$(CFo): %.o: %.cpp %.h
	$(MAKE) GPP=$(GPP) MAKEFLAGS='$(compiler_options)' -C $(CFDIR)

#clean up
clean:
	rm *.o 

cleanall:
	rm *.o *~ $(exe)


