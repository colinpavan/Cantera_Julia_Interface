# Cantera_Julia_Interface
# cpavan 2023-05-19

This project is a julia interface for Cantera
I have been building it and adding functionality as I need it, but currently in only contains some of the functionality
Examples of how to use can be found in the demo folder.
It is currently only tested on Linux.

To use, the library libcantera_shared.so must exist and its parents folder added to the dll search path (see demo file).
This library can be constructed by directly compiling the cantera c++ source code
follow instructions on: https://cantera.org/

