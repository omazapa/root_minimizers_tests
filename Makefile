CXX = g++
CCFLAGS = $(shell root-config --cflags) -ggdb


LD = g++
LDFLAGS = 

LIBS = $(shell root-config --libs) -lMLP -lMinuit -lTreePlayer -lTMVA -lTMVAGui -lXMLIO -lScipy  -lMLP -lm $(python3-config --ldflags --embed)

scipy:
	$(CXX) scipy.C $(CCFLAGS) -o test $(LIBS)

