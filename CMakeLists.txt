cmake_minimum_required(VERSION 3.5)

project(ChordDetectorAndChromagram)

list(APPEND CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -Wno-sign-compare -Woverloaded-virtual -Wwrite-strings -Wno-unused")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE -O3)
set(CMAKE_CXX_STANDARD 17)

set(KISSFFT
	libs/kiss_fft130/kiss_fft.c
)

set(SOURCES
	src/ChordDetector.cpp
	src/ChordDetector.h
	src/Chromagram.cpp
	src/Chromagram.h
)

add_library(chromagram ${SOURCES} ${KISSFFT})
target_include_directories(chromagram PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src libs/kiss_fft130/)
target_compile_definitions(chromagram PUBLIC -DUSE_KISS_FFT)
