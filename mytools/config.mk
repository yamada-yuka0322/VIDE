CC= gcc
CXX= g++
CPPFLAGS=  
LDFLAGS= -L/usr/local/lib -lCosmoTool -lgsl -lgslcblas
CXXFLAGS= $(CPPFLAGS) -ggdb -O0 -ftree-vectorize -ftree-vectorizer-verbose=2 -ffast-math
CFLAGS= $(CPPFLAGS) -ggdb -O0 -ftree-vectorize -ftree-vectorizer-verbose=2 -ffast-math
