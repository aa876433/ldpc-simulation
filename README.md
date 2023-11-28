# LDPC Simulation

## Overview
This is a C language implementation of an LDPC (Low-Density Parity-Check) decoder. LDPC codes are a type of error-correcting codes widely used in data communication and storage systems to improve the reliability of data transmission.

## Features
- Supports various LDPC decoding algorithms, including SPA (Sum-Product Algorithm), Layered SPA, MS (Min-Sum Algorithm), and Layered MS Algorithm.
- Provides functions for matrix simplification, encoding, decoding, and error rate calculation.
- Supports different Signal-to-Noise Ratio (SNR) settings and configurable number of decoding iterations.

## Usage
1. **Initialize LDPC Decoder**: Referencing the H.txt format, use the ldpc_init function to initialize the decoder.
2. **Set Configuration**: Set parameters such as SNR, number of iterations and decode algorithm using the `ldpc_set_config` function.
4. **Calculate Error Rates**: Use the `ldpc_simulation` function to calculate the bit error rate and block error rate.

## Compilation and Execution
- **Compile the Code**: Compile the source code using a C compiler (such as gcc).
- **Execute**: Run the compiled executable and perform tests.

## Notes
- Ensure understanding of each parameter and function in the code before use.
- Modifications and adjustments may be necessary depending on specific application requirements.
