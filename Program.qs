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
        
        for (dIdx in 1..n) {
            let qIdx = n-dIdx;
            H(inQubits[qIdx]);
            for (d in 1..(qIdx)) {
                Message("d={d}, qIdx={qIdx}" );
                Controlled R1Frac([inQubits[qIdx-d]], (d, 1, inQubits[qIdx]));
            }
        }
    }

    operation modExpU(inQubits: Qubit[], ancilla: LittleEndian,  modulus: Int, baseInt: Int, nQubits: Int): Unit {
        // computes b^x mod N
        
        mutable factor = baseInt; // b^(2^0) = b

        for (bitOffset in 0..nQubits-1) {
            if (bitOffset > 0) {
                set factor = ModI(factor*factor, modulus); // (b^x)^2 = b^(2x)
            }
            Controlled MultiplyByModularInteger([inQubits[bitOffset]], (factor, modulus, ancilla));
        }
    }



    // @EntryPoint()
    operation Main() : Unit {
        
        let n = 3;
        using (qubits = Qubit[n]) {

            // H(qubits[1]);

            Message("Pre-QFT State: ");
            DumpMachine();
            QFT(qubits, n);
            Message("Post-QFT State: ");
            DumpMachine();

            ApplyToEach(Reset, qubits);
        }
    }

    @EntryPoint()
    operation TestModMult() : Unit {
        
        let n = 4;
        using (qubits = Qubit[2*n]) {

            X(qubits[0]); // Set to |1> which we will multiply on
            X(qubits[5]); // Set control qubit to 1. 
            X(qubits[4]);

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

            ApplyToEach(Reset, qubits);
        }
    }

}
