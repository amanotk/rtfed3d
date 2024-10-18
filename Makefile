# -*- Makefile -*-
.PHONY: default clean cleanall

include Makefile.inc

# target
LIBS	= libcommon.a
SRCS	= global.cpp bc_periodic.cpp solver.cpp mc2.cpp integrator.cpp \
	io_hdf5.cpp hdf5util.cpp
OBJS	= $(SRCS:%.cpp=%.o)
DEPS	= $(SRCS:%.cpp=%.d)

default: $(LIBS)

libcommon.a: $(OBJS)
	$(AR) r $@ $(OBJS)

test: $(LIBS)
	$(MAKE) -C ./test

clean:
	rm -f *.o *.out *.a
	$(MAKE) -C ./test clean

cleanall:
	rm -f *.o *.d *.a *.dat *.pyc
	$(MAKE) -C ./test cleanall

# dependency
-include $(DEPS)
