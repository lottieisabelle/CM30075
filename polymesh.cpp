/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "polymesh.h"

using namespace std;

PolyMesh::PolyMesh(char *file)
{
  Transform *transform = new Transform();

  this->do_construct(file, transform);
}

PolyMesh::PolyMesh(char *file, Transform *transform)
{
  this->do_construct(file, transform);
}

void PolyMesh::do_construct(char *file, Transform *transform)
{
  // open the file
  std::ifstream f_reader(file);
  std::string line;

  // read the lines of the file
  if(f_reader.is_open()){
    // parse the header
    std::getline(f_reader,line);
    // printf("Header: %s\n",line.c_str()); TODO REMOVE DEBUG LINE

    // parse the number of vertices line
    std::getline(f_reader,line);
    // printf("Line 2: %s\n",line.c_str()); TODO REMOVE DEBUG LINE

    std::string vertices;
    for (int i = 15; i < 20; i++) {
      vertices = vertices + line[i];
    }
    vertex_count = std::stoi(vertices);
    cout << "vertex_count : " << vertex_count << "\n";

    // parse the number of faces line
    std::getline(f_reader,line);
    // printf("Line 3: %s\n",line.c_str()); TODO REMOVE DEBUG LINE

    std::string triangles;
    for (int i = 13; i < 18; i++) {
      triangles = triangles + line[i];
    }
    triangle_count = std::stoi(triangles);
    cout << "triangle_count : " << triangle_count << "\n";

    // loop through and process all lines that correspond to vertex coordinates
    for (int i = 0; i < vertex_count; i++) {
      
      triangles = triangles + line[i];
    }

    // QUESTION : do the numbers in the triangle lines correspond to lines in the file or index of coordinate in list of that data alone?

  }

  f_reader.close();
}
