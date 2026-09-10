# Generated instrument sounds

Four-second, 48 kHz, mono PCM16 WAVs from the default parameters of the
[NLD13–NLD15 ODE notebooks](../src/Pluto/03_ODE).

| File | Excitation | Frequency scale |
| --- | --- | --- |
| [reed_resonator.wav](reed_resonator.wav) | Rayleigh self-oscillation, one bore mode | 220 Hz |
| [bowed_resonator.wav](bowed_resonator.wav) | Smooth bow friction, one body mode | 196 Hz |
| [struck_resonator.wav](struck_resonator.wav) | Initial strike, cubic spring, second mode | 220 Hz |

The recorded coordinate is the resonator displacement proxy `q`. Coupling and
nonlinearity shift actual frequencies from the reference scale. Files include
attack and decay; each is low-pass filtered from a 4× sampled solution, stripped
of DC, peak-normalized to 0.7, and given 10 ms fades. Their normalized loudness
should not be used to compare the model's physical amplitudes.

Regenerate from the repository root:

```sh
julia --project=. scripts/render_audio.jl
```

The same models are available with sliders, symbolic MTK equations, direct ODE
functions, and live RealTimeAudioDiffEq controls in the four notebook tracks.
