all: compile execute

compile:
	c++ -o main.out main.cpp

execute:
	./main.out Hello 10 -1 2 -1
