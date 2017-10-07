CCPP = g++

OPTMIZATION_FLAGS = -O2 -march=native -DWITH_SSE
WARNING_FLAGS = -Wall -Wdate-time -Wformat -Werror=format-security  # -Wextra
DEBUGGING_FLAGS = #-ggdb3

CFLAGS = -fsigned-char -fno-strict-aliasing -D_FORTIFY_SOURCE=2 -fstack-protector-strong -fPIC $(OPTMIZATION_FLAGS) $(WARNING_FLAGS) $(DEBUGGING_FALGS) -Iinclude

TARGET_LIB = libctypesCPM.so # target lib

# linking flags
LDFLAGS = -lopencv_core -lopencv_highgui
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

