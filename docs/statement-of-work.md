# Statement of Work: Open-Source Notched Constant-Envelope Waveform Design

**Organization:** Open Research Institute (ORI)
**Project code name:** Eigenmannia
**Date:** November 21, 2026
**Status:** Draft

## 1. Background

At the Radio Club of America (RCA) Technical Symposium in New York on November 21,
2026, Spencer Markowitz presented "Optimal Waveform Design in the Presence of
Interference." He showed a method for designing wideband radar waveforms that avoid
transmitting in protected bands (about 30 dB down), keep range sidelobes low, keep fine
range resolution, and stay constant envelope. The method designs the power spectrum
for a desired autocorrelation as a convex problem, then realizes it with a
constant-envelope waveform by alternating projections.

In Q&A, he said the design currently takes a couple of seconds per waveform. Running it
pulse to pulse would need better numerical methods, and prior knowledge of the
interference would help.

That gap is a good fit for ORI. A real-time, open-source hardware implementation would
benefit radar researchers, spectrum-sharing policy work, and amateur radio
experimenters. The same mathematics applies to communications waveforms that must
avoid occupied spectrum.

## 2. Objective

Build an independent, open-source implementation of spectrally notched,
constant-envelope waveform design, demonstrate it in real time on ORI FPGA hardware,
and offer the results to the speaker and the community.

## 3. Principles

- **Independence.** We do not use or claim the speaker's code or unpublished results.
  The repository documents where we heard the talk and our own interpretation.
- **Credit.** The methodology as presented is credited to the speaker. Errors of
  interpretation are ours.
- **Value to the speaker.** Share the working design with him early and invite his review
  and comments. The goal is something useful to him, not something competing with him.
- **Labeled assumptions.** Every parameter he did not state is marked as an assumption.

## 4. Scope and phases

### Phase 0: Python reconstruction (first run done November 21)

- Convex power-spectrum design (cvxpy) and alternating-projection realization (NumPy).
- Compare against baselines: plain linear FM, jump over the band, alternating
  projections to a flat notched target.
- Metrics: notch depth (mean and worst case), peak and integrated sidelobe level,
  main-lobe width, envelope ripple, interference rejection.

### Phase 1: Refine and characterize

- Add an integrated-sidelobe or smoothness term so the convex target is realizable.
  The first run produced spiky targets at band and notch edges.
- Sweep time-bandwidth product, notch width and position, and main-lobe zone width
  (spoiling) to map notch depth against sidelobe level.
- Define and document an interference-tolerance metric. Attempt to reproduce the
  speaker's "19 dB more interference" figure, and state any difference in method.
- Add a mismatched receive filter (least-squares) as an option, and compare.
- Gather the references cited on the speaker's slides.

### Phase 2: Real-time algorithm

- Measure convergence: iterations needed against notch depth.
- Warm-start each pulse from the previous solution.
- Precompute convex spectrum designs for a library of likely notch patterns, and run
  only the fast alternating-projection step live.
- Use a fixed iteration count with a quality threshold for deterministic timing.
- Target: a new waveform within one pulse repetition interval.

### Phase 3: FPGA implementation

- Implement the alternating-projection loop (FFT, magnitude replacement, inverse FFT,
  constant-envelope projection) in FPGA fabric on the ZCU102.
- Feed the notch mask from a spectrum sensor. An ORI channelizer (as in Haifuraiya) is
  a candidate for supplying prior knowledge of the interference.
- Transmit and receive through the ADRV9002. **Constraint to check:** as far as we
  know, the ADRV9002 supports up to about 40 MHz of bandwidth, not 150 MHz. At 40 MHz,
  range resolution is about 3.75 m. Wider bandwidth needs a different radio front end.

### Phase 4: Demonstration

- Loopback and cabled tests with an injected narrowband interferer.
- Over-the-air demonstration in an amateur band, if Part 97 permits the emission.
  Pulse and radar-type emissions need to be checked against the band rules before
  any transmission.
- Write up the work for Inner Circle, and possibly as a QEX or conference paper.

## 5. Deliverables

| # | Deliverable | Phase |
| --- | --- | --- |
| D1 | Public repository with README, provenance and interpretation of the method | 0 |
| D2 | Python reference implementation with baselines and metrics | 0 to 1 |
| D3 | Characterization report: trade-off curves, interference tolerance | 1 |
| D4 | Real-time algorithm specification and timing measurements | 2 |
| D5 | FPGA design (VHDL), testbench, and build instructions | 3 |
| D6 | Hardware demonstration video and data | 4 |
| D7 | Article or paper | 4 |

## 6. Open questions

- Code name: chosen (Eigenmannia, after the glass knifefish whose jamming avoidance response shifts its electric signal away from a neighbor's frequency so both can keep sensing).
- License (ORI's usual choice for software and hardware designs).
- Whether to contact the speaker now, after Phase 1, or after a hardware demonstration.
- Choice of radio front end if 150 MHz of bandwidth is wanted.
- Volunteers: DSP and Python (Phases 0 to 2), VHDL (Phase 3), RF test (Phase 4).

## 7. Risks

- The constant-envelope realization may fall well short of the convex optimum at
  practical time-bandwidth products. Phase 0 already showed about 4 dB of loss and
  shallow notch edges.
- Real-time convergence may need more iterations than one pulse interval allows,
  making precomputed libraries necessary.
- Regulatory limits on radar emissions in amateur bands may restrict over-the-air
  testing.
