all:	main
main:	main.o
	g++ -Wall -Wextra -pedantic -pthread -g main.o -o k-nn
main.o:	main.cpp KDTree.hpp
	g++ -Wall -Wextra -pedantic -pthread -g -c main.cpp
clean:
	rm -f *.o k-nn
