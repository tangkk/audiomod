# audiomod

A compact C++ toolkit for **audio effects** and **audio analysis**.

`audiomod` combines classic DSP building blocks—time stretching, pitch shifting, modulation, dynamics, reverb, EQ—with analysis utilities such as **pYIN pitch tracking**, **global tuning estimation**, **key detection**, **chord estimation**, and **f0-to-note conversion**.

It can be used as:
- a command-line tool for offline WAV processing
- a static/shared library embedded in your own app
- a reference project for learning audio DSP

## Features

### Audio effects
- Time stretch
- Pitch shift / formant shift / gender change
- Vocoder / robotic / whisper voice effects
- Vibrato / tremolo / chorus / flanger / delay / phaser / autowah
- Compressor / limiter / autogain
- Reverb / EQ / gain

### Audio analysis
- Loudness measurement
- Envelope extraction
- pYIN note / pitch tracking
- Global tuning estimation
- Key detection
- Chord estimation
- f0-to-note conversion
- Humming estimation helper

## References

This project borrows ideas and code from:
- [Rubber Band](https://github.com/breakfastquay/rubberband)
- [Audio Effect Book Code](https://code.soundsoftware.ac.uk/projects/audio_effects_textbook_code/repository)
- [DAFX Book Code](https://www.dafx.de/DAFX_Book_Page_2nd_edition/matlab.html)

## License

Core parts are derived from and modified from Rubber Band, so this project is released under the **GNU GPL**. See [`COPYING`](COPYING).

---

## Build

### macOS / Linux

```bash
mkdir -p build
cd build
cmake ..
make
```

This builds:
- `audiomod-exe` — main CLI
- `audiomix-exe` — simple audio mixing CLI
- `libaudiomod.a` — static library

### Library-only build

```bash
mkdir -p build
cd build
cmake .. -DBUILD_LIBS=ON
cmake --build . --config Release
```

On Linux / Android, this can also produce a shared library.

### Docker

A CentOS-based Docker build environment is included:

```bash
docker build -t audiomod-build .
docker run --rm -it -v "$PWD":/app audiomod-build bash
```

Then build normally inside the container.

### Mobile builds

- Android: see `build-android.sh`
- iOS: see `build-ios.sh`

---

## CLI usage

```bash
./audiomod-exe <module> <input> <output> <args...>
```

Notes:
- Most effect modules read **WAV** and write **WAV**.
- Several analysis modules write **text output** instead of audio.
- `f0tonote` reads a **text input file** instead of WAV.
- If input is `-`, the program reads WAV data from **stdin**.

Run without enough arguments to see the current module list compiled into `main/main.cc`.

---

## Examples

### Pitch shift

```bash
./audiomod-exe normal_pitchshift stereo.wav out.wav 4 1 2048
```

Shift pitch up by 4 semitones.

### Time stretch

```bash
./audiomod-exe time_stretch in.wav out.wav 1.25 1 2048
```

Stretch audio to 1.25× duration.

### Reverb

```bash
./audiomod-exe reverb in.wav out.wav 0.8 0.9 0.5 0.9 0.1
```

### Loudness meter

```bash
./audiomod-exe loudnessmeter in.wav loudness.txt
```

Writes a single LUFS value.

### pYIN pitch tracking

```bash
./audiomod-exe pyin vocal.wav pyin.txt 2048 256 0.15
```

Output format:

```text
start,duration,f0
```

### Key detection

```bash
./audiomod-exe keydetection in.wav key.txt 8192 8192 -1
```

### Chord estimation

```bash
./audiomod-exe chordestimate in.wav chords.txt 8192 8192
```

### f0 to note conversion

```bash
./audiomod-exe f0tonote f0.csv notes.txt 44100 256 0.15
```

Input format:

```text
rms,f0
```

Output format:

```text
start,duration,pitch
```

---

## Library usage

Include `audiomod.h`, create the module you want, then process audio block by block.

```cpp
modbase *m_modbase = nullptr;
std::unique_ptr<phasevocoder> m_phasevocoder;

m_phasevocoder = std::make_unique<phasevocoder>(
    sample_rate,
    num_channels,
    1,
    shift_pitch,
    NORMAL_SHIFT,
    coremode,
    fftsize
);

m_modbase = m_phasevocoder.get();
```

Processing loop:

```cpp
for (int i = 0; i < file_length; i += block_size) {
    int num_samples = input->read(buff, block_size);
    m_modbase->processBlock(buff, num_samples);
    if (m_modbase->outputReady()) {
        output->write(buff, num_samples);
    }
}
```

For more complete examples, see:
- `main/main.cc`
- `include/audiomod.h`

---

## Project layout

- `main/` — CLI entrypoints and WAV helpers
- `include/` — public headers
- `src/` — DSP implementations
- `cmake/` — build helpers
- `Dockerfile` — reproducible Linux build environment

---

## Notes

The source of truth for supported modules and CLI arguments is currently `main/main.cc`. If you add or change modules, update this README accordingly.
