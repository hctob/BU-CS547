all:	main
main:	main.o
	g++ -Wall -Wextra -pedantic -pthread -g main.o -o integrate

main.o:	main.cpp
	g++ -Wall -Wextra -pedantic -pthread -g -c main.cpp

clean:
	rm -f *.o integrate
