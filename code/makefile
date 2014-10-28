CXX=g++

CXXFLAGS+=-c -Wall -Werror -fPIC -Wno-unused-local-typedefs -O2 -std=gnu++11 -pthread -fopenmp -I $(MONGO_ROOT)/src/mongo/client -I $(FLANN_ROOT) -I $(EIGEN_ROOT) -frounding-math
LDFLAGS+=-L $(MONGO_ROOT) -lmongoclient -lboost_thread -lboost_filesystem -lboost_program_options -lboost_system -lgomp -lm -lpthread -lCGAL -lCGAL_Core -lgmp

SOURCES=main.cpp metric.cpp dbscan.cpp mongo.cpp hist.cpp eval.cpp graph.cpp clustering.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=dbscan

FLANN_ROOT=/home/seydanator/Downloads/flann-1.8.4-src/src/cpp/
MONGO_ROOT=/home/seydanator/Downloads/mongo-cxx-driver-v2.4
EIGEN_ROOT=/home/seydanator/Downloads/eigen-eigen-3.2/
all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) $(MONGO_ROOT)/libmongoclient.a -o $@
.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	@- $(RM) $(EXECUTABLE)
	@- $(RM) $(OBJECTS)

distclean: clean
