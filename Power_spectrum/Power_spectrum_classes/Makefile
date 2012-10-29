GPP = mpicxx   #parallel compilation
#GPP = mpig++       #single processor

objects = $(patsubst %.cpp,%.o, $(wildcard *.cpp))
classh = $(wildcard *.h)

all: $(objects)

#create objects for the files in this directory
$(objects): %.o : %.cpp $(classh)
	$(GPP) -c $< -o $@ $(MAKEFLAGS) $(INCLUDE)

