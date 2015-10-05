

all: oddtool
	make clean_main

oddtool: sources/ODD_output.o sources/ODD_input.o sources/ODD_stepper.o sources/ODD_solution.o sources/ODD_utils.o sources/ODD_interpolator.o  main.o 
	g++ sources/ODD_output.o sources/ODD_input.o sources/ODD_stepper.o sources/ODD_solution.o sources/ODD_interpolator.o sources/ODD_utils.o main.o -o oddtool -L /usr/lib/x86_64-linux-gnu/ -l boost_regex 
	
sources/ODD_output.o: sources/ODD_output.cpp
	g++ -c -O2 -std=c++11 sources/ODD_output.cpp -o sources/ODD_output.o

sources/ODD_input.o: sources/ODD_input.cpp
	g++ -c -O2 -std=c++11 sources/ODD_input.cpp -o sources/ODD_input.o

sources/ODD_stepper.o: sources/ODD_stepper.cpp
	g++ -c -O2 -std=c++11 sources/ODD_stepper.cpp -o sources/ODD_stepper.o

sources/ODD_solution.o: sources/ODD_solution.cpp
	g++ -c -O2 -std=c++11 sources/ODD_solution.cpp -o sources/ODD_solution.o

sources/ODD_utils.o: sources/ODD_utils.cpp
	g++ -c -O2 -std=c++11 sources/ODD_utils.cpp -o sources/ODD_utils.o

sources/ODD_interpolator.o: sources/ODD_interpolator.cpp
	g++ -c -O2 -std=c++11 sources/ODD_interpolator.cpp -o sources/ODD_interpolator.o
	
	
main.o: main.cpp
	g++ -c -O2 -std=c++11 main.cpp -o main.o
	
	
clean:
	rm sources/*.o
	rm *.o

clean_main:
	rm *.o




