C_SOURCES = $(wildcard main.c lib/*.c)
# -lm flag links the math library
# -f flag deletes if present

run: clean main.out
	@./main.out

main.out: ${C_SOURCES}
	@gcc $^ -o main.out -lm 

opt: clean clang
	@./main.out

clang: ${C_SOURCES}
	@clang $^ -o main.out -lm 

clean:
	@rm -f main.out
