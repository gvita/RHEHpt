CUBA_ROOT := $(CUBA_DIR)

ROOTCONFIG   := root-config

ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)

LHAPDFCONFIG   := lhapdf-config

LHAPDFCFLAGS   := $(shell $(LHAPDFCONFIG) --cflags)
LHAPDFLIBS     := $(shell $(LHAPDFCONFIG) --libs)


MYLIBS := -L./lib -L/usr/local -lRHEHpt
MYFLAGS := -I./src

GSLLIBS := -lgsl -lgslcblas -lm

CUBAFLAGS := -I$(CUBA_ROOT)
CUBALIBS := -L$(CUBA_ROOT) -lcuba -lm

HQTLIBS := -lgfortran -lhqt

OBJECTS := $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp))

CXX := g++

OPTIM := -O3

DEBUG_FLAGS  := -O0 -g

PROFILE_FLAGS := -pg

C11 := -std=c++11

C14 := -std=c++1y

CXXFLAGS := $(CXXFLAGS) $(MYFLAGS) $(C11)  $(CUBAFLAGS) $(LHAPDFCFLAGS) -fPIC $(DEBUG_FLAGS) # $(OPTIM)
#$(ROOTCFLAGS)
#$(DEBUG_FLAGS) -Wall

DEBUG := -DDEBUG
GDB := -ggdb
	
lib/libRHEHpt.so: $(OBJECTS)
	$(CXX) -shared -fPIC -o $@ $^
#	ar -r $@ $^
	
lib/libhqt.so:
	cd hqt; \
	make libhqt.so; \
	mv libhqt.so ../lib/.
	
%.o: %.cpp %.h
	${CXX} ${CXXFLAGS} -c $< -o $@

%.x: %.o lib/libRHEHpt.so lib/libhqt.so
	${CXX} ${CXXFLAGS} -o $@ $< ${MYLIBS} ${HQTLIBS} ${GSLLIBS} ${CUBALIBS} ${LHAPDFLIBS}

clean:
	rm -f *.o *~ *.a *.x src/*.o src/*~ lib/*;\
	cd hqt;\
	make clean;
