# Eigenmannia: Notched Constant-Envelope Waveforms

An open-source, independent reconstruction of spectrally notched, constant-envelope
radar waveform design, with the goal of running it in real time on Open Research
Institute (ORI) field-programmable gate array (FPGA) development boards.

## Where this came from

On November 21, 2026, at the Radio Club of America (RCA) Technical Symposium in
New York, Spencer Markowitz presented "Optimal Waveform Design in the Presence of
Interference" (slide footer: "Spectral notching for chirp waveforms"). The talk
showed how to design a wideband radar waveform that:

1. keeps its transmitted spectrum about 30 dB down in protected bands where other
   users operate,
2. keeps range sidelobes low, so small targets are not hidden next to large ones,
3. keeps a narrow main lobe (about 1 m range resolution from 150 MHz of bandwidth), and
4. keeps a constant envelope, for maximum energy on target from a saturated amplifier.

This repository is **not** the speaker's code and does not claim to represent his
results. It is our interpretation of the methodology as we heard it, written so that
it can be checked, improved and, we hope, useful to him and to others. Parameters he
did not state are marked as assumptions in the code.

## Our interpretation of the method

**The problem with the obvious approaches** (as presented):

- A linear FM chirp has a sinc-shaped range response. A window lowers the sidelobes
  but widens the main lobe.
- A strong narrowband user inside the band raises the sidelobe floor and hides weaker
  targets.
- Jumping the chirp over a protected band creates sharp transitions that spread energy
  everywhere, giving only about 8 dB of attenuation.

**Method 1: alternating projections.** Alternate between two constraints until the
notch target is met:

1. FFT, set the spectral magnitude to the desired notched shape, keep the phase.
2. Inverse FFT, set the time-domain magnitude to a constant, keep the phase.

Stop on step 2, so the transmitted waveform is exactly constant envelope. This reaches
the notch, but says nothing about the autocorrelation, so the range sidelobes suffer.

**Method 2: design the spectrum for the autocorrelation you want.** The power
spectral density g and the autocorrelation r are a Fourier pair, so constraints on
both are linear in g and the design problem is convex:

- minimize the peak sidelobe level s
- subject to r = IDFT(g), r₀ = 1, |r_m| ≤ s outside a main-lobe zone of half-width K,
  r_m = 0 beyond the pulse length, g ≥ 0, and g at or below the notch level in
  protected bands.

Then realize g with a constant envelope using the alternating projections of method 1.

**"Spoiling"** (widening K) lets the optimizer trade resolution for lower sidelobes in
the best possible way: the optimal version of choosing a window.

**Reported result** (speaker's claim, not reproduced yet): with a matched receiver,
the notched link survives 19 dB more interference.

**What we still have to guess:** frequency-grid oversampling, main-lobe zone width,
initialization, stopping rule, and how the 19 dB figure was measured.

## Status

Phase 0 (Python reconstruction): first run complete. See
[docs/statement-of-work.md](docs/statement-of-work.md) for the plan.

First results (N = 256 samples, 300 MHz complex sampling, 150 MHz sweep,
25 MHz protected band at +30 MHz, K = 4; our parameters, not the speaker's):

| Design | Mean notch (dB) | Peak sidelobe (dB) | 3 dB width (m) | Interference rejection vs. linear FM (dB, our metric) |
| --- | --- | --- | --- | --- |
| Linear FM | 0.1 | −13.5 | 1.0 | 0 |
| Jump over band | −8.1 | −10.8 | 1.0 | 7.6 |
| Alternating projections, flat target | −16.5 | −9.8 | 1.0 | 15.8 |
| Convex target + alternating projections | −15.8 | −23.7 | 2.0 | 15.2 |

All four are exactly constant envelope. Early observations:

- The jump-over approach gives about 8 dB, matching the figure in the talk.
- At this short time-bandwidth product (128), neither notched design reaches −30 dB
  on average, and the notch edges are much shallower than the center.
- The convex design lowers the near-in peak sidelobe by about 10 dB, at the cost of
  2 m resolution (K = 4) and a higher far-sidelobe floor. The convex optimum was
  −28.1 dB, so realization loses about 4 dB.
- The convex target has sharp spikes at band and notch edges that a constant-envelope
  waveform cannot follow. A smoothness or integrated-sidelobe term is the next step.
- Doubling the pulse to N = 512 (time-bandwidth product 256) did **not** deepen the
  notch: alternating projections stayed near −15 to −16 dB mean. The shortfall
  against the speaker's −30 dB is therefore not just time-bandwidth product. Our
  projection loop, masks or notch metric likely differ from his. This is the first
  thing to resolve.
- The convex solve took 7 s at N = 256 and 69 s at N = 512 with a general-purpose
  solver. The realization step took under 0.2 s. This supports the real-time plan:
  precompute the convex designs and run only the projection step live.

## Running it

```bash
pip install numpy scipy matplotlib cvxpy
python prototype/notched_waveform.py --N 256 --K 4 --outdir results
```

Writes metrics as JSON to standard output and a spectrum and range-response plot to
`results/`.

## License and attribution

Open source (license to be chosen by ORI). Credit for the methodology as presented
belongs to Spencer Markowitz; any errors of interpretation are ours.
