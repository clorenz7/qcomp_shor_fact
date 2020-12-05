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
        // In-place computation of the quantum fourier transform 
        //  Transform taken via sequential controlled of Rd gates       
        // Note that for littleEndian input, the result will be big endian
        for (dIdx in 1..n) {
            let qIdx = n-dIdx;
            H(inQubits[qIdx]);
            for (d in 1..(qIdx)) {
                Controlled R1Frac([inQubits[qIdx-d]], (d, 1, inQubits[qIdx]));
            }
        }
    }

    operation modExpU(inQubits: Qubit[], ancilla: LittleEndian,  modulus: Int, baseInt: Int, nQubits: Int): Unit {
        // computes b^x mod N
        // via bitwise modular exponentiation
        
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

    operation measureModularFrequency(inQubits: Qubit[], ancilla: LittleEndian,  modulus: Int, baseInt: Int, nQubits: Int): Int {
        // Performs modular exponentiation, takes the QFT and measures the qubits

        // Perform the modular exponentiation
        modExpU(inQubits, ancilla,  modulus, baseInt, nQubits);
        Message("Modular Exponentiation Complete!");

        // Take the QFT of the input Qubits
        QFT(inQubits, nQubits);
        Message("QFT Complete!");

        // Measure, but need to account for the swap that the QFT Does
        let measInt = MeasureInteger(
            BigEndianAsLittleEndian(BigEndian(inQubits))  
        );  

        return measInt;
    }   

    operation ShorsFactoringAlgorithm(N: Int, nQubits: Int, baseInt: Int): (Int, Int) {
        // Factors N based on Shor's quantum factoring algorithm
        //  nQubits should be such that 2^nQubits approx N^2. 
        //  baseInt: a randomly selected number coprime to N
        let modulus = N;
        mutable measInt = 0;
        mutable foundPeriod = false;
        mutable period = 0;
        mutable j = 0;
        mutable count = 1;

        repeat {
            Message($"Running Period Finding Attempt #{count}");
            using (qubits = Qubit[2*nQubits]) {

                let chunks = Chunks(nQubits, qubits);
                // Set the ancilla qubit to $|1> for multiplication 
                X(chunks[1][0]);
                let x = chunks[0];
                // Prepare maximum superposition
                ApplyToEach(H, x);
                let y = LittleEndian(chunks[1]); 
                Message("State is prepared!");

                set measInt = measureModularFrequency(x, y,  modulus, baseInt, nQubits);
                Message($"Measured Frequency: |{measInt}>");

                ApplyToEach(Reset, qubits);  // Uncompute
            }

            set (period, j) = estimatePeriodWithContinuedFrac(measInt, nQubits, N);
            Message($"Period Measured to be: {period}, index: {j}");

            // check that period is not odd, 0 and that modular exp is consistent
            let checkVal = ExpModI(baseInt, period, N);
            set foundPeriod = (
                ModI(period, 2) == 0 and
                period != 0 and 
                checkVal == 1
            );
            set count += 1;
            if (not foundPeriod) {
                Message("Measured Period is inconsistent, retrying ...");
            }

        } until (foundPeriod); 

        let baseToHalfPeriod = PowI(baseInt, DividedByI(period,2));

        let factor1 = gcd(N, baseToHalfPeriod-1); // a^(r/2)-1
        let factor2 = gcd(N, baseToHalfPeriod+1);

        return (factor1, factor2);

    }

    operation gcd(a: Int, b: Int): Int {
        // Computes the greatest common divisor of two numbers
        Message($"Computing GCD({a}, {b})");

        // Initialize and find the larger number
        mutable factor = 0;
        mutable holder = 0;
        mutable small = b;
        mutable large = a;

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
        } until (large == small or small == 0);

        return large;
    }

    operation continuedFracAsRatio(contFrac: Int[]): (Int, Int) {
        // Converts a continued fraction represented as [a,b,c,d]
        // = 1/(a+1/(b+1/(c+1/d)))) as a simple fraction
        mutable num = 1;
        mutable nFrac = Length(contFrac);
        mutable denom = contFrac[nFrac-1];
        mutable holder = 0;

        for ( fromEnd in 2.. nFrac) {
            // Put to a common denominator and invert
            // e.g. 1/(a + n/d) = d/(a*d+n)
            set holder = denom;
            set denom = num + denom*contFrac[nFrac-fromEnd];
            set num = holder;
        } 

        return (num, denom);
    }


    operation estimatePeriodWithContinuedFrac(y: Int, nQubits: Int, N: Int): (Int, Int) {
        // Uses the measured quantum state to estimate the period
        // by estimating y/2^nQubits as j/period using continued fractions
        //   e.g. y = x_0 + j*period
        // Returns: (period, j)

        // Make sure to not divide by zero! 
        if ( y == 0 ) {
            return (0, 0);  // period of 0 will get picked up and retried later
        }
        
        mutable num = PowI(2, nQubits);
        mutable denom = y;
        // Compute floating point of fraction and stopping criteria
        let actual = IntAsDouble(y)/IntAsDouble(num);
        let stopThresh = 1.0/(2.0*IntAsDouble(num));

        // Initialize variables
        mutable contFracRep = new Int[0];
        mutable holder = 0;
        mutable factor = 0;
        mutable delta = 1.0;
        mutable j = 0;
        mutable period = 0;
        
        // Iteratively create a continued fraction representation until:
        //   It is a perfect representation OR we are very close to actual floating point
        repeat {
            // n/d = 1/(f+ (d-n)/n) where f = floor(d/n) -> assumes n < d
            set factor = DividedByI(num, denom);
            set contFracRep += [factor];
            set holder = denom;
            set denom = num - denom*factor;
            set num = holder;

            // Calculate the estimated fraction so far
            set (j, period) = continuedFracAsRatio(contFracRep);

            if (denom == 0) {  // perfect representation
                set delta = 0.0;
            } else {
                // See how far the current floating point fraction is from the actual
                set delta = AbsD(IntAsDouble(j)/IntAsDouble(period) - actual );
            }

        } until ( delta < stopThresh and period < N);

        return (period, j);
    }


    // --------- UNIT TESTS ---------------------
    //  (Verify that each function works correctly)

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
            let baseInt = 7;

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

            let period = measureModularFrequency(x, y,  modulus, baseInt, n);

            Message($"Period Measured to be: |{period}>");
            ApplyToEach(Reset, qubits);
        }
    }

    // @EntryPoint()
    operation TestContFrac() : Unit {
        // mutable contFrac = [1,3]; // 1/(1+ 1/3) = 3/4
        mutable contFrac = [4,12,4]; // 1/(4+ 1/(12+ 1/4))) = 49/200
        mutable num = 0; 
        mutable denom=0;
        set (num, denom) = continuedFracAsRatio(contFrac);
        Message($"Frac is: {num} / {denom}");

        // let (period, j) = estimatePeriodWithContinuedFrac(3, 2, 15); // = 3/4 expect 4, 3
        let (period, j) = estimatePeriodWithContinuedFrac(48, 6, 15); // = 48/64 expect 4,3
        Message($"Period, index is: {period} , {j}");
    }

    // @EntryPoint()
    operation TestGCD(): Unit {
        mutable val=1;
        set val = gcd(808, 4343); // Should be 101
        Message($"GCD is:{val}");
    }

    @EntryPoint()
    operation TestShors15() : Unit {
        // Test that 15 = 3x5 using 6 qubits by computing (7^x mod 15)
        let nQubits = 6;
        let baseInt = 7;
        let (factor1, factor2) = ShorsFactoringAlgorithm(15, nQubits, baseInt);
        Message($"Factors are: {factor1} and {factor2}");
    }

    // @EntryPoint()
    operation TestShors21() : Unit {
        // Test that 21 = 7x3 using 7 qubits by computing (5^x mod 21)
        let nQubits = 7;  // Using less qubits than suggested, more likely to retry, but faster attemps
        let baseInt = 5;
        let (factor1, factor2) = ShorsFactoringAlgorithm(21, nQubits, baseInt);
        Message($"Factors are: {factor1} and {factor2}");
    }

}