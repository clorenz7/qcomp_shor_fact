namespace ShorsFactoringAlgorithm {
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;

    operation QFT(inQubits: Qubit[], n: Int) : Unit  {
        // Computes the quantum fourier transform
        for (dIdx in 1..n) {
            let qIdx = n-dIdx;
            H(inQubits[qIdx]);
            for (d in 1..(qIdx)) {
                // Message("d={d}, qIdx={qIdx}" );
                Controlled R1Frac([inQubits[qIdx-d]], (d, 1, inQubits[qIdx]));
            }
        }
    }

    operation modExpU(inQubits: Qubit[], ancilla: LittleEndian,  modulus: Int, baseInt: Int, nQubits: Int): Unit {
        // computes b^x mod N
        // Via bitwise modular exponentiation
        
        mutable factor = baseInt; // b^(2^0) = b
        mutable factor = baseInt; // e.g. b^(2^0) = b

        for (bitOffset in 0..nQubits-1) {
            if (bitOffset > 0) {
                set factor = ModI(factor*factor, modulus); // (b^x)^2 = b^(2x)
                // Square the factor as we go up the qubits
                set factor = ModI(factor*factor, modulus); // e.g. (b^x)^2 = b^(2x)
            }
            // Only apply the unitary multiplication if the control qubit is |1>
            Controlled MultiplyByModularInteger([inQubits[bitOffset]], (factor, modulus, ancilla));
            Message($"Finished ModExp of bitOffset {bitOffset}");

        }
    }

    operation periodFinding(inQubits: Qubit[], ancilla: LittleEndian,  modulus: Int, baseInt: Int, nQubits: Int): Int {

        // Perform the modular exponentiation
        modExpU(inQubits, ancilla,  modulus, baseInt, nQubits);
        Message("Modular Exponentiation Complete!");

        // Take the QFT of the input Qubits
        QFT(inQubits, nQubits);
        Message("QFT Complete!");

        // Measure
        // But need to account for the swap that the QFT Does
        let measInt = MeasureInteger(
            BigEndianAsLittleEndian(BigEndian(inQubits))  // Was LittleEndian(inQubits)
        );  

        // Return the result. 
        return measInt;
    }   

        }
    }


    // --------- UNIT TESTS ---------------------

    // @EntryPoint()
    operation TestQFT() : Unit {
        
        let n = 3;
        using (qubits = Qubit[n]) {

            // H(qubits[1]);
            H(qubits[1]);

            Message("Pre-QFT State: ");
            DumpMachine();
            QFT(qubits, n);
            Message("Post-QFT State: ");
            DumpMachine();

            mutable qBE = BigEndian(qubits);
            let qLE = BigEndianAsLittleEndian(qBE);

            ApplyToEach(Reset, qubits);

            H(qubits[1]);
        
            Message("Pre-QFT State: ");
            DumpMachine();
            QFTLE(qLE);
            Message("Post-QFTLE State: ");          

            DumpMachine();
            ApplyToEach(Reset, qubits);

        }
    }

    @EntryPoint()
    // @EntryPoint()
    operation TestModMult() : Unit {
        
        let n = 4;
        using (qubits = Qubit[2*n]) {

            X(qubits[0]); // Set to |1> which we will multiply on
            X(qubits[5]); // Set control qubit to 1. Will multiply by 7^2
            X(qubits[4]); // Set control qubit. Will multiply by 7

            Message("Pre-Multiply State: ");
            DumpMachine();  // Expect 33 = 2^5 + 2^0
                        
            let yLE = LittleEndian(Subarray([0,1,2,3], qubits));
            let x = Subarray([4,5,6,7], qubits);

            let modulus = 15;
            // let baseInt = IntAsBigInt(7);
            let baseInt = 7;

            // modMultU(x, yLE, modulus, baseInt, 1);
            modExpU(x, yLE, modulus, baseInt, n);
            X(qubits[5]);  // # Flip control qubits back. 
            X(qubits[4]);
            Message("Post-Multiply State: ");
            DumpMachine();  // Expect 4 since 7^2 mod 15 = 4, or 13 since 7^3 mod 15 = 13

            let result = MeasureInteger(yLE);
            Message($"Post-Multiply and Measure State: |{result}>");
            // DumpMachine();  // Expect 4 since 7^2 mod 15 = 4, or 13 since 7^3 mod 15 = 13           

            ApplyToEach(Reset, qubits);
        }
    }

    // @EntryPoint()
    operation TestPeriodFinding() : Unit {

        let n = 6;
        let modulus = 15;
        let baseInt = 7;
        using (qubits = Qubit[2*n]) {

            let chunks = Chunks(n, qubits);
            // Set the ancilla qubit to $|1> for multiplication 
            X(chunks[1][0]);
            let x = chunks[0];
            let LX = Length(x);
            Message($"Length of x is: {LX}");
            // Prepare maximum superposition
            ApplyToEach(H, x);
            let y = LittleEndian(chunks[1]); 
            Message("State is prepared!");

            let period = periodFinding(x, y,  modulus, baseInt, n);

            Message($"Period Measured to be: |{period}>");
            ApplyToEach(Reset, qubits);
        }
    }

}
