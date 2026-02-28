1. Description
	MyTm is designed for fast melting point calculation. MyTm can achieve the commonly used MD direct simulation melting point method。


2. Requirements
	System:		Linux
	compiler:		intel icc/icx
	library:		GSL (v. 2.8 or later) and MKL (intel 2019 or later)


3. Project Structure
	src/			source files
	examples/		test examples
	bin/			 pre-compiled program
	Makefile		Makefile file
	README		Installation Guide

	
3. Installation Process
	(a) If the GSL and MKL libraries are not indexed by the system, please set GSL_DIR, GSL_LIB, MKL_DIR, and MKL_LIB in the Makefile
	(b) make clean
	(c) make 

	The bin directory provides a pre compiled program that can be used directly.



4. Command
	$ MyTm -r 
	Perform a melting point calculation.
	$ MyTm -h 
	Get command instructions.
	$ MyTm -m method
	Output the scheduled input.dat file.
	method （integer）:		0 Direct-heating method
						1 Void method
						2 Modified void method
						3 Two-phase method
						4 Sandwich method
						5 Z method
						6 Modified Z method

