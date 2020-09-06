CXX ?= g++
AR  ?= ar
MAKE ?= make

LIBS := -lm

CFLAGS := -I./ \
		  -I../ \
		  -I../* \
		  -Wl,--no-as-needed

LDFLAGS += $(shell pkg-config opencv --cflags --libs)

OBJS := *.o */*.o

.PHONY: clean build start all

all : clean build

build:
	$(CXX) $(CFLAGS) $(LIBS) $(LDFLAGS) -c *.cpp 

clean:
	rm -f *.o

