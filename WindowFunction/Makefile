GPP = g++

gsl = -lgsl -lgslcblas -lm      #gsl options
#options = -Wall -O2 # -wd981      #std options
options = -Wall -g             #std options

windowmat: windowmat.cpp
	$(GPP) -o windowmat windowmat.cpp $(gsl) $(options)
