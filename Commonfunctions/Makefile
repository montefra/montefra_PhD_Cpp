GPP = mpicxx   #parallel compilation
#GPP = mpig++       #single processor

objects = $(patsubst %.cpp,%.o, $(wildcard *.cpp))

all: $(objects)

#create objects for the files in this directory
$(objects): %.o : %.cpp %.h     
	$(GPP) -c $< -o $@ $(MAKEFLAGS)

