#!/bin/bash
g++ -std=c++11 src/* -pthread -o ising; ./ising; cat results
