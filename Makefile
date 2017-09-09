CCPP = g++

SSE_FLAG = -msse2 -DWITH_SSE

CFLAGS = -O3 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -fstack-protector-strong -Wformat -Werror=format-security -fPIC -march=native $(SSE_FLAG) -Iinclude
#CFLAGS = -w -ggdb3  -Iinclude
LDFLAGS = -lopencv_core -lopencv_highgui


TARGET_LIB = libctypesCPM.so # target lib

SHARED_LDFLAGS = -shared -Wl,-soname,$(TARGET_LIB)# linking flags


SOURCES_CPP := $(shell find . -name '*.cpp')
OBJ := $(SOURCES_CPP:%.cpp=%.o)
HEADERS := $(shell find . -name '*.h')

all: CPM $(TARGET_LIB)

.cpp.o:  %.cpp %.h
	$(CCPP) -c -o $@ $(CFLAGS) $+

%.pylib.o:  %.cpp
	$(CCPP) -c -o $@ $(CFLAGS) -DPYLIB $+

	
CPM: $(HEADERS) CPM.o main.o
	$(CCPP) -o $@ CPM.o main.o $(LDFLAGS)


$(TARGET_LIB): CPM.pylib.o ctypesCPM.pylib.o
	$(CCPP) $(CFLAGS) ${SHARED_LDFLAGS} -o $@ $^

clean:
	rm -f $(OBJ) *.o CPM $(TARGET_LIB)

