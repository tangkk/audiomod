# run this from the build/ dir
cmake .. -G "Unix Makefiles"
cmake --build . --config Release
cp libaudiomod.a ../lib/mac/
# cp libaudiomod.dylib ../lib/mac/
