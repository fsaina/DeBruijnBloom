cmake_minimum_required(VERSION 3.13)
project(DeBruijnBloom)

set(CMAKE_CXX_STANDARD 14)

add_executable(DeBruijnBloom main.cpp BloomFilter.cpp BloomFilter.h MurmurHash3.cpp MurmurHash3.h ExactDeBruijnGraph.cpp ExactDeBruijnGraph.h KmerUtil.cpp KmerUtil.h)
