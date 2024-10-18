# -*- Makefile -*-

# include path
TOPDIR = $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BASE   = $(TOPDIR)/../base
BLITZ  = $(HOME)/usr/include

# compilers
CXX      = mpicxx
CXXFLAGS = -O3 -MMD -I$(BLITZ) -I$(BASE) -I$(TOPDIR)
LDFLAGS  = -L$(HOME)/usr/lib -lblitz

%.o	: %.cpp
	$(CXX) -c $(CXXFLAGS) $<
