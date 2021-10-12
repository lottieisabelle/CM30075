#!bin/bash

# Compile Lab 2
g++ -o lab2executable main_lab2.cpp framebuffer.cpp linedrawer.cpp polymesh.cpp -lm
# Execute stage
./lab2executable

echo "Image saved."