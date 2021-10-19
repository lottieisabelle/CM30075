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

int get_vertex_count(string line)
{
  stringstream ss(line);
  string word1;
  string word2;
  int vertex_count;
  ss >> word1 >> word2 >> vertex_count;
  return vertex_count;
}

int get_face_count(string line)
{
  stringstream ss(line);
  string word1;
  string word2;
  int face_count;
  ss >> word1 >> word2 >> face_count;
  return face_count;
}

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

    // parse the number of vertices line
    std::getline(f_reader,line);
    vertex_count = get_vertex_count(line);
    cout << "vertex_count : " << vertex_count << "\n";

    // parse the number of faces line
    std::getline(f_reader,line);
    triangle_count = get_face_count(line);
    cout << "triangle_count : " << triangle_count << "\n";

    // create arrays to store vertices and triangles
    PolyMesh::vertex = new Vertex[PolyMesh::vertex_count];
    PolyMesh::triangle = new TriangleIndex[PolyMesh::triangle_count];

    // loop through and process all lines that correspond to vertex coordinates
    for (int i = 0; i < vertex_count; i++) {
      std::getline(f_reader,line);
      int axis_count = 0;
      istringstream ss(line.c_str());
      float x;
      float y;
      float z;

      string value;
      while (ss >> value){
        if (axis_count == 0){
          x = std::stof(value);
        } else if (axis_count == 1){
          y = std::stof(value);
        } else if (axis_count == 2){
          z = std::stof(value);
        }
        axis_count++;
      }
      // store vertex data in array
      PolyMesh::vertex[i] = Vertex(x, y, z);

      // apply transform
      transform->apply(PolyMesh::vertex[i]);

    }
    cout << "All vertices processed. \n";

    // loop through and process all lines that correspond to triangle coordinate pointers
    for (int i = 0; i < triangle_count; i++) {
      std::getline(f_reader,line);
      int axis_count = 0;
      istringstream ss(line.c_str());

      string index_value;
      while (ss >> index_value){
        // store triangle face vertex data in array
        if(axis_count != 0){

          int number = std::stoi(index_value);
          number -= 1;
          PolyMesh::triangle[i][axis_count-1] = number;
        }
        axis_count++;
      }
    }
    cout << "All faces processed. \n";
  }

  f_reader.close();
}
