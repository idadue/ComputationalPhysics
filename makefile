all: compile execute

compile:
	c++ -o main.out main.cpp

execute:
	./main.out -1 2 -1
