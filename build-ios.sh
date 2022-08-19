# run this from the build/ dir
cmake .. -G Xcode -DCMAKE_TOOLCHAIN_FILE=/path/to/ios-cmake/ios.toolchain.cmake -DPLATFORM=OS -DBUILD_LIBS=ON
cmake --build . --config Release
cp Release-iphoneos/libaudiomod.a ../tmpbuild/
if [[ "$PWD" =~ "audiomod/build" ]]
then
    echo "remove all in build folder for another build"
    rm -rf *
fi

cmake .. -G Xcode -DCMAKE_TOOLCHAIN_FILE=/path/to/ios-cmake/ios.toolchain.cmake -DPLATFORM=SIMULATOR64 -DBUILD_LIBS=ON
cmake --build . --config Release

lipo -create ../tmpbuild/libaudiomod.a Release-iphonesimulator/libaudiomod.a -output libaudiomod.a

cp libaudiomod.a ../lib/ios/
# cp Release-iphoneos/libaudiomod.dylib ../lib/ios/
