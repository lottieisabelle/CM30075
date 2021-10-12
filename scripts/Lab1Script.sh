#!bin/bash

# Compile Lab 1
g++ -o lab1executable main_lab1.cpp framebuffer.cpp linedrawer.cpp -lm
# Execute stage
./lab1executable

echo "Image saved."