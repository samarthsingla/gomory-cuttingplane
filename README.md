# Gomory Cutting Plane Method
This project implements the *Gomory Cutting plane* method for solving *Integer Linear Programming problems* in canonical form 

The code assumes the following format:
maximise (c'x) 
Subject to 1) Ax <= b, 2) x >= 0, 3) x is an integer  
(x is the variable vector, c' indicates the transpose of c vector) 

Input format:
- The first line of input contains two space-separated integers n and m - The dimensions of c
and b, respectively.
- The second line has m space-separated integers - b1,b2,...,bm
- The third line has n space-separated integers - c1, c2,..., cn : The cost vector
- The ith line of the next m lines contains the ith row of the coefficient matrix A.

Requirements: Python3 and numpy  

Team Members:  
Samarth Singla  
Nishant Yadav
