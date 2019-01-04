cmake_minimum_required(VERSION 3.13)
project(DeBrujinBloom)

set(CMAKE_CXX_STANDARD 14)

#add_executable(DeBrujinBloom main.cpp BloomFilter.cpp BloomFilter.h hashFn.cpp hashFn.h base.h)
add_executable(DeBrujinBloom main.cpp StringBloomFilter.cpp StringBloomFilter.h MurmurHash3.cpp MurmurHash3.h)
