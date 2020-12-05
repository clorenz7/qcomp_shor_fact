# Introduction
Q# Implementation of Shor's Quantum Factoring Algorithm for my PHYS 575 Final Project

Shor's algorithm is based on finding the period of the unitary function b^x mod N using the Quantum Fourier Transform. 

This code implements the QFT using controlled Rd gates, and implements the unitary function based on the binary decomposition of x. It uses a library function in the Microsoft Quantum Software Development Kit to perform the modular multiplcations. 

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

