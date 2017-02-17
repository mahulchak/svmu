#Default CXX is g++
CXXFLAGS+=-Wall -std=c++0x

PROGRAMS=fasplitter svmu scriptmaker checkCNV findDel findInvert

.PHONY: all
all: $(PROGRAMS)

fasplitter: fasplitter.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

svmu: mlib.cpp svmu.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

scriptmaker: script_maker.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

checkCNV: cnvlib.cpp ccnv.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

findDel: dellib.cpp findInDel.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

findInvert: findInvert.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@
.PHONY: clean
clean:
	rm -f $(PROGRAMS)
