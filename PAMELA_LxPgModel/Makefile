ROOTINC := $(shell root-config --cflags)
ROOTLIB := $(shell root-config --libs)
LIBS         += $(ROOTLIB)
LIBS         += -lMinuit

CXX := g++ -g -fPIC -std=c++11 `root-config --cflags`
SOFLAGS := -shared -fPIC
ifeq ($(shell uname),Darwin)
SOFLAGS = -dynamiclib -single_module -undefined dynamic_lookup
endif

VPATH = num:src:main

# ---------- main targets

main: bin/Fit.exe
lib: lib/libNumerical.a
# ----------
NumericalSRC := DataPoints.cpp Spline3Interpolator.cpp FCMatrix.cpp EqSolver.cpp Vec.cpp FCMatrixBanded.cpp FCMatrixFull.cpp



NumericalOBJ := $(patsubst %.cpp,lib/%.o,$(notdir $(NumericalSRC)))

SRC := NTAMS02Data.cpp NTPAMELAData.cpp SolarModulation.cpp 
OBJ := $(patsubst %.cpp,lib/%.o,$(notdir $(SRC)))

# --------

lib/libNumerical.a: $(NumericalOBJ)
	ar rvcs $@ $^
	ranlib $@
	$(CXX) $(SOFLAGS) -o $(basename $@).so $^
	
bin/Fit.exe: NTFit.cpp $(OBJ) | lib/libNumerical.a
	$(CXX) $<  $(OBJ) $(ROOTINC) -I num -I src -L lib -l Numerical $(LIBS) -o $@
	
lib/%.o: %.cpp
	$(CXX) -I num $(ROOTINC) -c $< -o $@

clean:
	rm -f lib/*.o bin/Fit.exe
	rm -f lib/*.a lib/*.so






