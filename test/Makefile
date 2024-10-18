# -*- Makefile -*-

include ../Makefile.inc

TARGET = riemann1d vortex circular reconnection blastwave sheets tearing
OBJS   = $(addsuffix .o, $(TARGET))
LIBS   = $(TOPDIR)/libcommon.a

default: $(TARGET)

riemann1d: riemann1d.o $(LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

vortex: vortex.o $(LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

circular: circular.o $(LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

reconnection: reconnection.o $(LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

blastwave: blastwave.o $(LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

sheets: sheets.o $(LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

tearing: tearing.o $(LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

diffusion: diffusion.o $(LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f *.o *.out

cleanall:
	rm -f $(TARGET)
	rm -f *.o *.d *.out *.dat *.pyc

# dependency
$(OBJS): $(TOPDIR)/*.hpp $(TOPDIR)/*.cpp
