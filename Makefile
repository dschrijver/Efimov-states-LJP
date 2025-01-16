COMPILER = gcc
OPT_FLAGS = -std=c17 -O3
DEBUG_FLAGS = -Wall -Wextra -Wno-unused-result
LIBRARIES = -lm
C_SOURCES = $(wildcard main.c src/*.c)

all: clean main.out
	./main.out

main.out: ${C_SOURCES}
	$(COMPILER) $^ -o $@ $(OPT_FLAGS) $(DEBUG_FLAGS) $(LIBRARIES) 

clean:
	rm -f main.out
