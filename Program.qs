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
        
        mutable factor = baseInt; // e.g. b^(2^0) = b

        for (bitOffset in 0..nQubits-1) {
            if (bitOffset > 0) {
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
    operation gcd(a: Int, b: Int): Int {

        mutable factor = 0;
        mutable holder = 0;
        mutable small = a;
        mutable large = b;

        if (large < small) {
            set holder = small;
            set small = large;
            set large = holder;
        }

        // The basic strategy is to repeatedly subtract the smaller number from
        // the larger since cx-cy = c(x-y), eg GCD remains after subtraction. 
        repeat  {
            set factor = DividedByI(large, small); // Do division to speed things up
            set holder = small;
            set small = large - factor*small;
            set large = holder;
        } until (large == small);

        return small;
    }

    operation continuedFracAsRatio(contFrac: Int[]): (Int, Int) {
        mutable num = 1;
        mutable nFrac = Length(contFrac);
        mutable denom = contFrac[nFrac-1];
        mutable holder = 0;

        for ( fromEnd in 2.. nFrac) {
            set holder = denom;
            set denom = num + denom*contFrac[nFrac-fromEnd];
            set num = holder;
        } 

        return (num, denom);
    }


    operation estimatePeriodWithContinuedFrac(y: Int, nQubits: Int): (Int, Int) {
        mutable contFracRep = new Int[0];
        mutable num = PowI(2, nQubits);
        mutable denom = y;
        let actual = IntAsDouble(y)/IntAsDouble(num);
        mutable holder = 0;
        mutable factor = 0;
        mutable delta = 1.0;
        mutable j = 0;
        mutable period = 0;
        let stopThresh = 1.0/(2.0*IntAsDouble(num));

        repeat {
            set factor = DividedByI(num, denom);
            set contFracRep += [factor];
            set holder = denom;
            set denom = num - denom*factor;
            set num = holder;

            set (j, period) = continuedFracAsRatio(contFracRep);

            if (denom == 0) {
                set delta = 0.0;
            } else {
                set delta = AbsD(IntAsDouble(j)/IntAsDouble(period) - actual );
            }

        } until ( delta < stopThresh );

        return (period, j);
    }


    // --------- UNIT TESTS ---------------------

    // @EntryPoint()
    operation TestQFT() : Unit {
        
        let n = 3;
        using (qubits = Qubit[n]) {

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

    // @EntryPoint()
    operation TestContFrac() : Unit {
        // mutable contFrac = [1,3]; // 1/(1+ 1/3) = 3/4
        mutable contFrac = [4,12,4]; // 1/(1+ 1/3) = 49/200
        mutable num = 0; 
        mutable denom=0;
        set (num, denom) = continuedFracAsRatio(contFrac);
        Message($"Frac is: {num} / {denom}");

        // let (period, j) = estimatePeriodWithContinuedFrac(3, 2); // = 3/4 expect 4, 3
        let (period, j) = estimatePeriodWithContinuedFrac(48, 6); // = 48/64 expect 4,3
        Message($"Period, index is: {period} , {j}");
    }
}