# run this from the build/ dir
cmake .. -DCMAKE_TOOLCHAIN_FILE=/path/to/Android/sdk/ndk/versionxx.x.xxxx/build/cmake/android.toolchain.cmake \
    -DANDROID_NDK=/path/to/Android/sdk/ndk/versionxx.x.xxxx/ \
    -DANDROID_ABI="armeabi-v7a" \
    -DANDROID_STL="c++_shared" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_LIBS=ON
cmake --build . --config Release
cp libaudiomod.a ../lib/android-armv7/
cp libaudiomod.so ../lib/android-armv7/

if [[ "$PWD" =~ "audiomod/build" ]]
then
    echo "remove all in build folder for another build"
    rm -rf *
fi

# run this from the build/ dir
cmake .. -DCMAKE_TOOLCHAIN_FILE=/path/to/Android/sdk/ndk/versionxx.x.xxxx/build/cmake/android.toolchain.cmake \
    -DANDROID_NDK=/path/to/Android/sdk/ndk/versionxx.x.xxxx/ \
    -DANDROID_ABI="arm64-v8a" \
    -DANDROID_STL="c++_shared" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_LIBS=ON
cmake --build . --config Release
cp libaudiomod.a ../lib/android-arm64/
cp libaudiomod.so ../lib/android-arm64/