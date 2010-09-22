CC= gcc
CXX= g++
CPPFLAGS=  
LDFLAGS= -L/usr/local/lib -lCosmoTool -lgsl -lgslcblas
CXXFLAGS= $(CPPFLAGS) -ggdb -O3 -ffast-math
CFLAGS= $(CPPFLAGS) -ggdb -O3 -ffast-math
