CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)
Repo = $(shell git rev-parse --show-toplevel)

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include 
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR) -lRooFit -lRooFitCore -L$(INC)/include 

CPPFLAGS += -g -std=c++11

TARGET = testapp 
TARGET2 = HLTeffApp 

SRC = app/testapp.cc src/Pi0Events.cc src/PhoEvents.cc src/HLTeff.cc
SRC2 = app/HLTeffApp.cc src/Pi0Events.cc 

OBJ = $(SRC:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)

all : $(TARGET) $(TARGET2)

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET2) : $(OBJ2)
	$(LD) $(CPPFLAGS) -o $(TARGET2) $(OBJ2) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^


%.o : %.cc
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@
	@echo $<


clean :
	rm -f *.o app/*.o src/*.o $(TARGET) $(TARGET2) *~
