CCPP = g++


OPTMIZATION_FLAGS = -O2 -march=native
WARNING_FLAGS = -Wall -Wdate-time -Wformat -Werror=format-security  # -Wextra
DEBUGGING_FLAGS = #-ggdb3

CFLAGS = -fsigned-char -fno-strict-aliasing -D_FORTIFY_SOURCE=2 -fstack-protector-strong -fPIC $(OPTMIZATION_FLAGS) $(WARNING_FLAGS) $(DEBUGGING_FALGS) -Iinclude
LDFLAGS = -lopencv_core -lopencv_highgui

SOURCES_CPP := $(shell find . -name '*.cpp')
OBJ := $(SOURCES_CPP:%.cpp=%.o)
HEADERS := $(shell find . -name '*.h')

all: CPM

.cpp.o:  %.cpp %.h
	$(CCPP) -c -o $@ $(CFLAGS) $+

CPM: $(HEADERS) $(OBJ)
	$(CCPP) -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm -f $(OBJ) CPM

