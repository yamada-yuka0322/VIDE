PROGS= showZobov

CXX=g++
CC=gcc

include config.mk

SOURCES= showZobov.cpp loadZobov.cpp zobovConf.c 

LIBS= $(LDFLAGS) 

all: $(PROGS)

showZobov: showZobov.o loadZobov.o zobovConf.o

testZobov: testZobov.o loadZobov.o

depend: $(SOURCES)
	@echo "[DEPENDS] $^"
	@$(CC) $(CPPFLAGS) -M -MM $^ > .mydepends

distclean: clean
	@rm -f .mydepends

clean: 
	@rm -f *.o
	@rm -f $(PROGS)

.mydepends: $(SOURCES) 
	@touch .mydepends
	@make depend

%.prog:
	@echo "[L] $@"
	@$(CXX) -o $@ $^ $(LIBS)

zobovConf.c zobovConf.h: showZobov.ggo Makefile
	@echo "[OPT] $@"
	gengetopt  -i $<  -f zobovConf  -a zobovConf_info  -F zobovConf -C

%.o: %.c
	@echo "[C] $< ..."
	@$(CC) -c -o  $@ $< $(CFLAGS)

%.o: %.cpp
	@echo "[C++] $< ..."
	@$(CXX) -c -o $@ $< $(CXXFLAGS)

include .mydepends
