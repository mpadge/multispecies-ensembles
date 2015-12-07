CXX=clang++-3.5
CFLAGS=-c -std=c++11 
LIBS=-lboost_program_options
VPATH=./src
OBJECTS = utils.o pop-fns.o model.o
OBJECTS_050 = utils.o pop-fns.o model050.o
OBJECTS_THEORETICAL = model050-theoretical.o
OBJECTS_ENVCOR = utils.o pop-fns.o model-envcor.o
OBJECTS_TROPH = trophic-levels.o

all: model model050 theor envcor trophic

model: $(OBJECTS)
	$(CXX) $(OBJECTS) -o model $(LIBS) 

model050: $(OBJECTS_050)
	$(CXX) $(OBJECTS_050) -o model050 $(LIBS) 

theor: $(OBJECTS_THEORETICAL)
	$(CXX) $(OBJECTS_THEORETICAL) -o model-theor

envcor: $(OBJECTS_ENVCOR)
	$(CXX) $(OBJECTS_ENVCOR) $(LIBS) -o model-envcor

trophic: $(OBJECTS_TROPH)
	$(CXX) $(OBJECTS_TROPH) $(LIBS) -o trophic-levels

%.o: %.c++
	$(CXX) $(CFLAGS) $<

clean:
	rm -f *.o 
