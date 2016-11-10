# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ 

CPP_EX = prog 

all: display clean prog

display:
	clear
    
# Deletes all binary files
clean:
	rm -rf prog main.o
        
# Compilation
prog: main.o
	$(CCC) main.o -o prog 

main.o: main.cpp
	$(CCC) -c -Wall main.cpp -o main.o
   

