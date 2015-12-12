CUBA_ROOT := $(CUBA_DIR)

ROOTCONFIG   := root-config

ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)

LHAPDFCONFIG   := lhapdf-config

LHAPDFCFLAGS   := $(shell $(LHAPDFCONFIG) --cflags)
LHAPDFLIBS     := $(shell $(LHAPDFCONFIG) --libs)


MYLIBS := -L./lib -L/usr/local -lRHEHpt
MYFLAGS := -I./src

#LOOPTOOLSLIBS := -L/usr/local/lib -looptools -lgfortran -lm

GSLLIBS := -lgsl -lgslcblas -lm


CUBAFLAGS := -I$(CUBA_ROOT)
CUBALIBS := -L$(CUBA_ROOT) -lcuba -lm

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp))

CXX := g++

OPTIM := -O3

DEBUG_FLAGS  := -O0 -g

PROFILE_FLAGS := -pg

C11 := -std=c++11

C14 := -std=c++1y

CXXFLAGS := $(CXXFLAGS) $(MYFLAGS) $(C11)  $(CUBAFLAGS) $(LHAPDFCFLAGS)  $(OPTIM) # $(DEBUG_FLAGS)
#$(ROOTCFLAGS)
#$(DEBUG_FLAGS) -Wall

DEBUG := -DDEBUG
GDB := -ggdb


lib/libRHEHpt.a: $(OBJECTS)
	@echo $(OBJECTS)
	ar -r $@ $^

	
%.o: %.cpp %.h
	${CXX} ${CXXFLAGS} -c $< -o $@

%.x: %.o lib/libRHEHpt.a
	${CXX} ${CXXFLAGS} -o $@ $< ${MYLIBS} ${GSLLIBS} ${CUBALIBS} ${LHAPDFLIBS}

clean:
	rm -f *.o *~ *.a *.x src/*.o src/*~ lib/libRHEHpt.a;
