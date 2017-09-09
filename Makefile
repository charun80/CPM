CCPP = g++

SSE_FLAG = -msse2 -DWITH_SSE

CFLAGS = -O3 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -fstack-protector-strong -Wformat -Werror=format-security -fPIC -march=native $(SSE_FLAG) -Iinclude
#CFLAGS = -w -ggdb3  -Iinclude
LDFLAGS = -lopencv_core -lopencv_highgui

SHARED_LDFLAGS = -shared # linking flags

TARGET_LIB = libctypesCPM.so # target lib

SOURCES_CPP := $(shell find . -name '*.cpp')
OBJ := $(SOURCES_CPP:%.cpp=%.o)
HEADERS := $(shell find . -name '*.h')

all: CPM $(TARGET_LIB)

.cpp.o:  %.cpp %.h
	$(CCPP) -c -o $@ $(CFLAGS) $+

CPM: $(HEADERS) CPM.o main.o
	$(CCPP) -o $@ CPM.o main.o $(LDFLAGS)


$(TARGET_LIB): CPM.o ctypesCPM.o
	$(CCPP) $(CFLAGS) ${LDFLAGS} ${SHARED_LDFLAGS} -o $@ $^

clean:
	rm -f $(OBJ) CPM $(TARGET_LIB)

