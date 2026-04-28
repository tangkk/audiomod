cmake_minimum_required(VERSION 3.5.1)
project(audiomod)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake; ${CMAKE_MODULE_PATH};")


if (ANDROID) 

elseif(UNIX AND NOT APPLE) # LINUX
  set(LINUX TRUE)  
elseif(APPLE) # iOS or MAC
  
else()

endif()

file(GLOB audiomix_TEST
    "main/audiomix.cc"
    "main/wavfile.cc"
    "main/wavfile.h"
)

file(GLOB audiomod_TEST
    "main/main.cc"
    "main/wavfile.cc"
    "main/wavfile.h"
)

file(GLOB AE_SOURCES
    # audio module sources
    "src/phasevocoder/*.cc"
    "src/delayline/*.cc"
    "src/distortion/*.cc"
    "src/dynamics/*.cc"
    "src/filtering/*.cc"
    "src/modulation/*.cc"
    "src/sinusoidal/*.cc"
    "src/sourcefilter/*.cc"
    "src/spatial/*.cc"
    "src/reverb/*.cc"
    "src/equalizer/*.cc"
    "src/gain/*.cc"

    # meter sources
    "src/meter/*.cc"

    # common sources
    "src/common/gen/*.cc"
    "src/common/system/*.cc"
    "src/common/speex/*.c"
    "src/common/pommier/*.c"
    "src/common/kissfft/*.c"
    "src/common/dsp/*.cc"
    "src/common/base/*.cc"
    "src/common/filters/*.cc"
    "src/common/utils/*.c"
    "src/common/realtime/*.cc"
    "src/common/*.cc"
    "src/common/*.c"


    # analyzer sources
    "src/envelope/*.cc"
    "src/pyin/*.cc"
    "src/chromathing/*.cc"
    "src/chromathing/*.c"
    )

file(GLOB AE_HEADERS
    "include/*.h"
    "include/dafx/*.h"
    "include/analyzer/*.h"
    "include/common/*.h"
)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_library(audiomod STATIC ${AE_SOURCES})
SET_TARGET_PROPERTIES(audiomod PROPERTIES
   OUTPUT_NAME audiomod
   CLEAN_DIRECT_OUTPUT 1
   XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY ""
   )

if (ANDROID OR LINUX)
add_library(audiomodshared SHARED ${AE_SOURCES})
SET_TARGET_PROPERTIES(audiomodshared PROPERTIES
   OUTPUT_NAME audiomod
   CLEAN_DIRECT_OUTPUT 1
   XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY ""
   )
endif()

set (AE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src")
foreach (_headerFile ${AE_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND AE_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES AE_INCLUDE_DIRS)

include_directories(${AE_INCLUDE_DIRS})

option(BUILD_LIBS "Build the static library" OFF)

if(BUILD_LIBS)
 
if (ANDROID OR LINUX)
    set_property(TARGET audiomod audiomodshared PROPERTY CXX_STANDARD 14)
else()
    set_property(TARGET audiomod PROPERTY CXX_STANDARD 14)
endif()

else(BUILD_LIBS)

    add_executable(audiomod-exe ${audiomod_TEST})
    target_link_libraries(audiomod-exe audiomod)

    add_executable(audiomix-exe ${audiomix_TEST})
    target_link_libraries(audiomix-exe audiomod)

    if(ANDROID OR LINUX)
        set_property(TARGET audiomod-exe audiomix-exe audiomod audiomodshared PROPERTY CXX_STANDARD 14)
    else()
        set_property(TARGET audiomod-exe audiomix-exe audiomod PROPERTY CXX_STANDARD 14)
    endif()
    
    if(APPLE)
        # set(CMAKE_EXE_LINKER_FLAGS "-libc++")
        install (TARGETS audiomod-exe audiomix-exe DESTINATION bin)

    else(APPLE)
        # for static linkage you need to have static gcc installed, do
        # yum install glibc-static libstdc++-static
        set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
        install (TARGETS audiomod-exe audiomix-exe DESTINATION bin)
    endif(APPLE)
endif(BUILD_LIBS)


