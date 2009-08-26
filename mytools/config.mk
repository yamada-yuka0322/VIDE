CC= gcc
CXX= g++
CPPFLAGS=  -I$(HOME)/Science/Software/CosmoToolbox/install/include
LDFLAGS= -L$(HOME)/Science/Software/CosmoToolbox/install/lib -lCosmoTool
CXXFLAGS= $(CPPFLAGS) -ggdb -O0 -ffast-math
CFLAGS= $(CPPFLAGS) -ggdb -O0 -ffast-math
