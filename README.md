# PAM-4 Link Simulation with Least-Squares Equalization

A MATLAB simulation of a PAM-4 baseband link through a bandwidth-limited single-pole channel. It compares the received signal before and after an 8-tap least-squares (LS) equalizer, then measures symbol errors over an AWGN sweep.

## Contents

- `PAM4.m` - simulation source. The complex frequency-domain plot uses magnitude so MATLAB does not discard its imaginary component or print the complex-plot warning.
- This README - system description, receiver metrics, and results from the run supplied with this project.

The source PDF used to develop the technical explanation is not included.

## Signal chain

1. Generate 20,000 random bits and group them into 10,000 two-bit symbols.
2. Map the symbols to the four voltage levels `{-2, -2/3, +2/3, +2}`. The bit labels use natural binary order: `00`, `01`, `10`, `11`.
3. Hold each symbol for 20 samples. The sample rate is 40 kHz, so the symbol rate is 2 kbaud and the symbol period is 0.5 ms.
4. Apply a first-order low-pass channel with cutoff at one quarter of the rectangular pulse's main-lobe bandwidth. This bandwidth restriction rounds the transitions and spreads energy into adjacent symbol periods, creating intersymbol interference (ISI).
5. Apply an ideal receiver low-pass filter with bandwidth equal to the pulse main-lobe bandwidth.
6. Train an 8-tap, symbol-spaced FIR equalizer on the noiseless received waveform. Search all 20 sample phases and 8 target alignments; select the phase and alignment that minimize mean-squared error (MSE).
7. Add white Gaussian noise, filter and equalize it through the same receiver, sample at the selected phase, and slice to the nearest PAM-4 level.

## Equalizer and receiver metrics

The LS equalizer finds tap weights that make its output approximate the known transmitted symbol sequence. Searching the sampling phase selects where within each symbol interval the receiver should sample. Searching the target index shifts the desired sequence against the equalizer output; this lets the finite tap span balance precursor and postcursor ISI.

The selected sampling phase is **17** (MATLAB's 1-based sample index), and the target index is **5** of 8. The minimum training MSE is **3e-6**. The measured effective signal gain is **0.999997**, and the residual ISI power is **3e-6**. The near-unity gain means the desired symbol amplitude is essentially preserved; residual ISI is the remaining error after subtracting that gain-scaled ideal symbol sequence.

The reported **noise enhancement factor (NEF) is 0.195483**. In this script it is computed from the mean squared magnitude of the combined receiver low-pass filter and upsampled equalizer response. It scales the input noise power to estimate noise power at the decision point. An NEF below one means this filtered receiver chain reduces average white-noise power under the script's frequency-grid normalization; it is not a claim that an equalizer alone always reduces noise. Equalizers can boost frequencies that the channel attenuates, which can increase noise, so the preceding receiver filter and equalizer must be considered together.

The run reports pre-equalization SIR of **2.96 dB** and post-equalization SIR of **58.19 dB**, an improvement of about **55.23 dB**. SIR describes signal power relative to interference power; it excludes the added AWGN. The post-equalizer SIR uses the measured residual ISI power.

For the theoretical curve, the script combines residual ISI and NEF-scaled noise into an effective SINR, then uses the PAM-4 Q-function approximation `BER ≈ (1.5/2) Q(sqrt(0.2 × SINR))`. The reported simulated BER is computed as symbol-error rate divided by two. Since the mapping is natural binary rather than Gray coded, this conversion is an approximation: a decision between the two middle levels can change both bits. The theoretical curve is likewise an approximate comparison for this bit mapping, rather than an exact bit-error calculation.

## Results

Results below reproduce the supplied run. The error count is the number of symbol decisions that differ from the transmitted symbols, out of 10,000 symbols. The SINR values in the supplied console output are linear values despite its `dB` label; the table converts them using `10 log10(SINR)`.

| Eb/N0 (dB) | Symbol errors | Reported BER (SER/2) | Theoretical BER | SINR (dB, converted) |
|---:|---:|---:|---:|---:|
| 0  | 5660 | 2.830000e-1 | 2.809032e-1 | -2.91 |
| 2  | 5151 | 2.575500e-1 | 2.576943e-1 | -0.91 |
| 4  | 4568 | 2.284000e-1 | 2.295733e-1 | 1.09 |
| 6  | 3897 | 1.948500e-1 | 1.962533e-1 | 3.09 |
| 8  | 3143 | 1.571500e-1 | 1.581434e-1 | 5.09 |
| 10 | 2320 | 1.160000e-1 | 1.169201e-1 | 7.09 |
| 12 | 1497 | 7.485000e-2 | 7.608194e-2 | 9.09 |
| 14 | 796  | 3.980000e-2 | 4.084406e-2 | 11.09 |
| 16 | 318  | 1.590000e-2 | 1.634111e-2 | 13.09 |
| 18 | 88   | 4.400000e-3 | 4.149016e-3 | 15.09 |
| 20 | 17   | 8.500000e-4 | 5.181000e-4 | 17.09 |
| 22 | 1    | 5.000000e-5 | 2.122628e-5 | 19.09 |
| 24 | 0    | 0           | 1.500202e-7 | 21.09 |

The measured signal power is **2.246756**. At 24 dB, this finite 10,000-symbol run observed no symbol errors; that does not imply the true BER is zero. The theoretical estimate at that point is approximately **1.50e-7**, and a longer simulation would be needed to estimate such a low error rate reliably.

The simulated and theoretical curves follow the same decreasing trend across the sweep. Differences at high Eb/N0 are expected when only a small number of errors occur, and the natural-binary mapping and SER-to-BER approximation also affect their agreement.

## Run it

Open `PAM4.m` in MATLAB and run the script. It uses MATLAB functions including `upsample`, `qfunc`, and `eyediagram`; the Communications Toolbox may be needed for the latter functions, depending on the MATLAB installation. The bit stream and noise are randomly generated, so later runs will not reproduce the supplied error counts exactly.

The script displays waveform, BER, frequency-response, and eye-diagram figures and prints receiver and error-rate metrics to the Command Window.
