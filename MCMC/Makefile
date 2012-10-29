###############################################################
# Makefile to compile the MCMC code
#                                                              
# Default compilation: make                                   
# options:                                                     
#   GPP: to change the compiler [default: g++]                
#   optimize: optimization flag [default -g]                  
#   includes: directory containing includes for gsl, if not   
#     in standard path                                        
#   libs: directory containing library for gsl, if not in     
#     standard path                                           
#   PRECOPT=-DWINMAT: precompiler option to include the       
#     effect of the window matrix                             
###############################################################
#GPP = mpicxx
GPP = g++

# window matrix or not (if done with if might take longer)
#PRECOPT = -DWINMAT
#gsl options
gsl = -lgsl -lgslcblas -lm
#options = -Wall -O2 # -wd981      #std options
optimize = -g
options = -Wall $(optimize) $(PRECOPT)
debug = -fmudflap -lmudflap #-L/usr/lib  #debugging options

#contains functions often used functions
CFDIR =  /home/montefra/data1/Cpp/Commonfunctions
P22DIR = /data01/montefra/Cpp/Theories
INCLUDE = -I$(CFDIR) -I$(P22DIR) -I$(HOME)/include 
LINK = -L$(HOME)/lib 

#files in the current directory
cppfiles = $(wildcard *.cpp) 
ofiles = $(patsubst %.cpp,%.o, $(cppfiles)) 
hfiles = $(wildcard *.h)

oio = $(CFDIR)/io.o
cio = $(CFDIR)/io.cpp
hio = $(CFDIR)/io.cpp
op22 = $(P22DIR)/p22_class.o
cp22 = $(P22DIR)/p22_class.cpp
hp22 = $(P22DIR)/p22_class.h

ifeq ($(strip $(PRECOPT)),)
exe = mcmc #executables
else
exe = mcmc_w
endif


all: $(exe)

$(exe): $(ofiles) $(oio) $(hfiles) 
	$(GPP) $(options) -o $@ $(ofiles) $(oio) $(op22) $(gsl) $(LINK)

$(ofiles): %.o: %.cpp $(hfiles)
	$(GPP) -c $< -o $@ $(INCLUDE) $(options) 

#create the object files for common functions
$(oio): $(cio) $(hio)
	$(MAKE) GPP=$(GPP) MAKEFLAGS='$(compiler_options)' -C $(CFdir)

#create the object for the codes that computes p22
$(op22): $(cp22) $(hp22)
	$(MAKE) GPP=$(GPP) MAKEFLAGS='$(compiler_options)' -C $(CFdir)

#help: prints out the beginning of the Makefile
help:
	@head -n14 Makefile

#clean up
clean:
	rm *.o 

cleanall:
	rm *.o *~ $(exe)


