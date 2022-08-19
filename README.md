# audiomod quickstart

## Introduction

**audiomod** is a project for audio modifications, including audio manipulators such as time-stretching, pitch-shifing, formant-changing, and audio filters such as vibrato, tremolo, ring-modulation, compression, reverb, equalizer, etc. 
It is a good starting point for audio signal processing developpers to build more advanced applications, or for students to learn and practice their audio coding techniques.

## References
This project refers deeply to the algorithms and code of the following opensourced projects:
- [Rubberband](https://github.com/breakfastquay/rubberband)
- [Audio Effect Book Code](https://code.soundsoftware.ac.uk/projects/audio_effects_textbook_code/repository)
- [Dafx Code](https://www.dafx.de/DAFX_Book_Page_2nd_edition/matlab.html)
 
This project tries to unified these three sources into one framework,
and provide "build-abilities" for desktop and mobile platforms as well.

## License
Since the core part of the code are modified based on [Rubberband](https://github.com/breakfastquay/rubberband), 
the license of this project follows GNU General Public License (GPL).

## Build standalone executables for Linux and Mac

under the source folder do:
```
mkdir build
cd build
cmake ..
make
```
This will process an executable called `audiomod-exe`, as well as 
a static library `libaudiomod.a`.

usage:
```
usage: ./audiomod-exe dafx_name infile outfile <args> (dafx: constant, time_stretch, normal_pitchshift, formant_pitchshift, gender_change, vocoder, vocoder_chord, robotic, whisper, vibrato, delay, flanger, chorus, ringmod, tremolo, compressor, limiter, autogain, reverb, autowah, phaser, loudnessmeter, equalizer, gain, envelope, )
```
The `<args>` here can be looked up in `main.cc`. For example:
`./audiomod-exe normal_pitchshift stereo.wav out.wav 4 1 2048`
runs a 4-semitone up pitch-shift to stereo.wav, using a phasevocoder with a 2048-point FFT. Note that currently only `.wav` format is supported.


## Build for Android

For an android build, you need to use cmake's cross-platform build process, 
which depends on the [Android NDK](https://developer.android.com/ndk). Here is an example for arm64:
```
cmake .. -DCMAKE_TOOLCHAIN_FILE=/path/to/ndk/18.1.5063045/build/cmake/android.toolchain.cmake \
    -DANDROID_NDK=/path/to/ndk/18.1.5063045/ \
    -DANDROID_ABI="arm64-v8a" \
    -DANDROID_STL="c++_shared" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_LIBS=ON
	
cmake --build . --config Release
```

Please refer to `build-android.sh` for more details on this process.

## Build for iOS

For an ios build, you need to have the Xcode installed, and then make use of cmake's
cross-platform build process, which depends on the [ios-cmake](https://github.com/leetal/ios-cmake).
Here is an example:
```
cmake .. -G Xcode -DCMAKE_TOOLCHAIN_FILE=../ios-cmake/ios.toolchain.cmake -DPLATFORM=OS -DBUILD_LIBS=ON
cmake --build . --config Release
```
Please refer to `build-ios.sh` for more details on this process.


## SDK usage

To incorporate audiomod into your own project, first you need to build it for your target platform. After this step you get a `libaudiomod.a` or `libaudiomod.so` (for Android)

Include the `audiomod.h` in your source code. Then you create an object of the sound filter, and initialize it like this:
```
modbase * m_modbase = nullptr;

std::unique_ptr<phasevocoder> m_phasevocoder;

m_phasevocoder = std::make_unique<phasevocoder>(sample_rate, num_channels, 1, shift_pitch, NORMAL_SHIFT, coremode, fftsize);
        m_modbase = m_phasevocoder.get();
```
you should provide the parameters needed accordingly. Then you could process the sound with the following code:
```
for (int i = 0; i < file_length; i+=block_size) {
	int num_samples = input->read(buff, block_size);
	m_modbase->processBlock(buff, num_samples);
	if (m_modbase->outputReady()) {
		output->write(buff, num_samples);
	}
}
```
Refer to `main.cc` for more details on the exact procedure.

Finally, make sure to link the `libaudiomod` library in your building process.



