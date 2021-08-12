Ohmic circuit solver is a simple program that computes the current (I) value across each branch of a given (linear) circuit. The program utilises Graph data structure and various Graph algorithms for its efficient execution. 
The program takes the following circuit data as input:
•	number of nodes
•	branch connections
•	resistor and EMF positions, and their respective values. 

Proper exception handling has been implemented to catch faulty input data.

The input data is stored in the form of adjacency matrix.  The input is then processed using a combination of network theory techniques to obtain a system of linear equations. 

Briefly, the following steps are taken in order to achieve the same:
•	The graph network is converted into a Spanning Tree
•	Links are connected one at a time, to obtain a fundamental set of cycles 
•	Further pruning is done in order to isolate each cycle. Refer method isolateLoop()
•	Each cycle is then traversed, and the direction of current is set up simultaneously, i.e. the cycles are converted from undirected graph to a directed graph. 

The system of equation is then solved using Gauss elimination algorithm (refer method gauss() ) to obtain the current values. 

NOTE: The circuit element should only comprise of resistor and battery/emf. The program is not designed to handle capacitor, inductor, voltage source, current source or any non-linear element (like diode, transistor, etc.)

REFRENCES: https://core.ac.uk/download/pdf/53745212.pdf
