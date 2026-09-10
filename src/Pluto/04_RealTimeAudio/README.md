# RealTimeAudio

Define the ODE and a DESource, choose its output states, and change parameters, time scale, and gain while it runs.

Open these files with Pluto; their first cell activates and installs the shared
repository environment automatically. See the [course index](../../../README.md)
for the matching topics in the other tracks.

Playback is off on opening. Check **Play** to start and uncheck it to stop.
Start with a low speaker volume. Sliders change the existing source; they do not
recreate it. Stop and play again to repeat a strike. Closing a browser tab alone
may leave the notebook running: stop playback or shut down the notebook in Pluto.

These examples use the registered RealTimeAudioDiffEq 0.1 API in the shared
manifest. A local audio output device is needed for live playback; the other
tracks provide WAV playback and export. The audio callback solves a short interval
for each output buffer, using the state reached by the preceding buffer.

`channel_map=[3,3]` sends the resonator to both channels. For the introductory
oscillator and Duffing example, `[1,1]` sends the displacement. `set_ts!(source,
2pi*f0)` makes one unit-frequency cycle correspond to `f0` Hz; coupling and
nonlinearity can shift the resulting pitch. The introductory NLD03 uses zero
damping to sustain a tone. The live callback uses its selected gain without the
offline WAV filter or peak normalization.

NLD16 uses `channel_map=[[4,6],[4,6]]` to send the summed modal pressures to both
channels. Keep Ω₂ = 2.5 to compare its periodic and quasiperiodic blowing presets;
allow the transient to settle after changing a parameter. Its source problem sets
tighter solver tolerances before playback, using the pinned version 0.1 API.

The Play cell returns a PlutoHooks cleanup function that stops its source when
the cell is rerun or removed or the notebook is shut down. `list_devices()` and
`get_device_index(...)` from RealTimeAudioDiffEq can be used to choose an output
other than the default in that cell. The stream uses the device’s default sample
rate, while the WAV examples use 48 kHz.

See the [package API](https://github.com/antonioortegabrook/RealTimeAudioDiffEq.jl)
and [PlutoHooks documentation](https://juliapluto.github.io/PlutoHooks.jl/src/notebook.html).

- [NLD03_Flows2D](NLD03_Flows2D.jl)
- [NLD10_Flows2D_Forced](NLD10_Flows2D_Forced.jl)
- [NLD13_Reed_Resonator](NLD13_Reed_Resonator.jl)
- [NLD14_Bowed_Resonator](NLD14_Bowed_Resonator.jl)
- [NLD15_Struck_Resonator](NLD15_Struck_Resonator.jl)
- [NLD16_Two_Resonators](NLD16_Two_Resonators.jl)
