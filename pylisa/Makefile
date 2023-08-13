LDFLAGS := -lm
CFLAGS := -std=c99 -Wall -Wextra -pedantic -fPIC -O2

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

mie: mie.c
	$(GXX) $(CFLAGS) $(LDFLAGS) -shared $< -o mie.$(EXT)

clean:
	rm -rf *.o *.exe *.so *.dll

.PHONY: clean all