# Introduction
Q# Implementation of Shor's Quantum Factoring Algorithm for my PHYS 575 Final Project

Shor's algorithm is based on finding the period of the unitary function b^x mod N using the Quantum Fourier Transform. Once the a frequency is found, continued fraction expansion is used to find the period, and the Euclid's greatest common denomintor is used to find the actual factors. 

This code implements the QFT using controlled Rd gates, and implements the unitary function based on the binary decomposition of x. It uses a library function in the Microsoft Quantum Software Development Kit to perform the modular multiplcations. The continued fraction and GCD code is custom. 

There are quick printed out validation tests of the most functions. They can be run by changing the `@EntryPoint()` of the code. 

# Running the code

To run the code :
1. Install Microsoft [Visual Studio Code](https://code.visualstudio.com/download)
2. Install [.NET Code SDK  ](https://www.microsoft.com/net/download)
3. Install the [Microsoft Quantum Development Kit](https://marketplace.visualstudio.com/items?itemName=quantum.quantum-devkit-vscode) - should link to inside of VS Code
4. Clone this repo using git 
5. From the command line, cd into the repository folder and type `dotnet run`  

Updated instructions might be found [here](https://docs.microsoft.com/en-us/quantum/quickstarts/install-command-line?tabs=tabid-vscode)

## Example Output

Factoring 15:
```
PS C:\Users\DasCorCor\Dropbox\PHYS 575\code\ShorsFactoringAlgorithm> dotnet run
Running Period Finding Attempt #1
State is prepared!
Finished ModExp of bitOffset 0
Finished ModExp of bitOffset 1
Finished ModExp of bitOffset 2
Finished ModExp of bitOffset 3
Finished ModExp of bitOffset 4
Finished ModExp of bitOffset 5
Modular Exponentiation Complete!
QFT Complete!
Measured Frequency: |0>
Period Measured to be: 0, index: 0
Measured Period is inconsistent, retrying ...
Running Period Finding Attempt #2
State is prepared!
Finished ModExp of bitOffset 0
Finished ModExp of bitOffset 1
Finished ModExp of bitOffset 2
Finished ModExp of bitOffset 3
Finished ModExp of bitOffset 4
Finished ModExp of bitOffset 5
Modular Exponentiation Complete!
QFT Complete!
Measured Frequency: |16>
Period Measured to be: 4, index: 1
Computing GCD(15, 48)
Computing GCD(15, 50)
Factors are: 3 and 5
```

Factoring 21 (WARNING - Takes Awhile!): 