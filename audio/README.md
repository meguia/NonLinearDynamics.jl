# Generated instrument sounds

Four-second, 48 kHz, mono PCM16 WAVs from the default parameters of the
[NLD13–NLD16 ODE notebooks](../src/Pluto/03_ODE).

| File | Excitation | Frequency scale |
| --- | --- | --- |
| [reed_resonator.wav](reed_resonator.wav) | Rayleigh self-oscillation, one bore mode | 440 Hz |
| [bowed_resonator.wav](bowed_resonator.wav) | Smooth bow friction, one body mode | 196 Hz |
| [struck_resonator.wav](struck_resonator.wav) | Initial strike, cubic spring, second mode | 220 Hz |
| [two_resonators.wav](two_resonators.wav) | Nonlinear reed feedback, two acoustic modes, quasiperiodic preset | 185 Hz |

NLD13–NLD15 record the resonator displacement proxy `q`; NLD16 records the summed
modal pressure `p₁ + p₂`. Coupling and nonlinearity shift actual frequencies from
the reference scale. The first three files include their attack and decay; the
two-resonator example discards the first 5,000 units of dimensionless time to
expose the sustained regime. Each is low-pass filtered from a 4× sampled solution, stripped
of DC, peak-normalized to 0.7, and given 10 ms fades. Their normalized loudness
should not be used to compare the model's physical amplitudes.

Regenerate from the repository root:

```sh
julia --project=. scripts/render_audio.jl
```

The same models are available with sliders, symbolic MTK equations, direct ODE
functions, and live RealTimeAudioDiffEq controls in the four notebook tracks.
