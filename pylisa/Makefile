LDFLAGS := -lm
CFLAGS := -std=c99 -Wall -Wextra -pedantic -fPIC #-g

ifeq ($(OS_OPTION), windows)
	GXX := x86_64-w64-mingw32-gcc
	EXT := dll
else ifeq ($(OS_OPTION), darwin)
	GXX := clang -arch x86_64 -arch arm64
	EXT := dylib
else
	GXX := gcc
	EXT := so
endif

all: test_mie mie_so

test_mie: test_mie.o mie.o
	$(GXX) $^ $(LDFLAGS) -o $@

test_mie.o: test_mie.c mie.h
	$(GXX) $(CFLAGS) -c $< -o $@

mie.o: mie.c mie.h
	$(GXX) $(CFLAGS) -c $< -o $@

mie_so: mie.c mie.h
	$(GXX) $(CFLAGS) $(LDFLAGS) -shared $< -o mie.$(EXT)

clean:
	rm -rf *.o *.exe *.so *.dll *.dylib test_mie

.PHONY: clean all