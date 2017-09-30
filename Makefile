CCPP = g++

SSE_FLAG = -msse2 -DWITH_SSE

CFLAGS = -O2 -march=native -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -fstack-protector-strong -Wformat -Werror=format-security -fPIC -march=native $(SSE_FLAG) -Iinclude
#CFLAGS = -w -ggdb3  -Iinclude
LDFLAGS = -lopencv_core -lopencv_highgui


TARGET_LIB = libctypesCPM.so # target lib

# linking flags
#SHARED_LDFLAGS = -shared -Wl,-soname,$(TARGET_LIB) -Wl,--version-script=$(TARGET_LIB:%.so=%.version) 
SHARED_LDFLAGS = -shared -Wl,-soname,$(TARGET_LIB)


SOURCES_CPP := $(shell find . -name '*.cpp')
OBJ := $(SOURCES_CPP:%.cpp=%.o)
HEADERS := $(shell find . -name '*.h')

all: CPM $(TARGET_LIB)

.cpp.o:  %.cpp %.h
	$(CCPP) -c -o $@ $(CFLAGS) $+

%.pylib.o:  %.cpp
	$(CCPP) -c -o $@ $(CFLAGS) -fvisibility=hidden -DPYLIB $+

	
CPM: $(HEADERS) CPM.o main.o
	$(CCPP) -o $@ CPM.o main.o $(LDFLAGS)


$(TARGET_LIB): ctypesCPM.pylib.o
	$(CCPP) $(CFLAGS) ${SHARED_LDFLAGS} -o $@ $^

clean:
	rm -f $(OBJ) *.o CPM $(TARGET_LIB)

