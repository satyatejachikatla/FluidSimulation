CXX ?= g++
AR  ?= ar
MAKE ?= make

LIBS := -lm

CFLAGS := -I./ \
		  -I./Fluid \
		  -Wl,--no-as-needed

LDFLAGS += $(shell pkg-config opencv --cflags --libs)

OBJS := *.o */*.o

.PHONY: clean build start all

all : clean build start

build:
	$(MAKE) -C Fluid
	$(CXX) $(CFLAGS) $(LIBS) -c main.cpp $(LDFLAGS)
	$(CXX) $(CFLAGS) $(LIBS) $(OBJS) -o run $(LDFLAGS)

start:
	./run

clean:
	$(MAKE) -C Fluid clean
	rm -f *.o
	rm -f run
	rm *.jpg

