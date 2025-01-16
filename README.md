# Efimov state LJP
This project aims to find the Efimov ground state of a three-particle system with an LJP interaction. Machine learning is applied to achieve this, by using a small artificial neural network to efficiently store the three-particle wave function. 

## Running the program
- The code is explicitely written for Linux.
- The only required programs are gcc and make. On Ubuntu, these can be installed using:
```console
sudo apt install build-essential
```
- Download the repository and move into the root directory containing the [Makefile](Makefile).
- To run the program, simply run the command ```make``` in your terminal.

## Changing settings
All important settings can be changed in [constants.h](constants.h) in the root directory. 
