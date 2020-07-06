all:	main

main:	main.o
	g++ -Wall -Wextra -pedantic -fopenmp -g part1-redux.o -o d-sphere

main.o:	part1-redux.cpp
	g++ -Wall -Wextra -pedantic -fopenmp -g -c part1-redux.cpp

part2:	part2.o
	g++ -Wall -Wextra -pedantic -fopenmp -g part2.o -o simd

part2.o:	part2.cpp
	g++ -Wall -Wextra -pedantic -mavx -g -c part2.cpp
clean:
	rm -f *.o simd-vec simd part1 part2 d-sphere
