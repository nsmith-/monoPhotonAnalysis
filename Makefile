CXX           = g++
CXXFLAGS      = $(OPT2) $(ROOTCFLAGS) -Wall -g -O2 -std=c++11
LD            = g++
LDFLAGS       = $(OPT2)

SRCS = monoPhotonAnalysis.cc monoPhotonAnalyzer.cc

ROOTINCLUDES = -I$(ROOTSYS)/include
ROOTLDLIBS = `$(ROOTSYS)/bin/root-config --cflags --glibs` -lEG

INCLUDES = -I. $(ROOTINCLUDES)
LDLIBS = -L. $(ROOTLDLIBS)

all : mpa
	@echo Built $@

.PHONY: clean

clean :
	rm -rf *.o

mpa: $(SRCS:.cpp=.o)
	$(LD) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) $^ $(LDLIBS) -o $@
	@echo Built $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -fPIC $^ 
	@echo Compiled $@

