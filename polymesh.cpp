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

    // TODO : make more error proof for other .ply files
    //std::string vertices;
    //for (int i = 15; i < 20; i++) {
      //vertices = vertices + line[i];
    //}
    // vertex_count = std::stoi(vertices);

    // cout << "vertex_count : " << vertex_count << "\n";

    // different method
    //vertex_count = get_vertex_count(line);
    vertex_count = 3644;
    cout << "vertex_count : " << vertex_count << "\n";


    // parse the number of faces line
    std::getline(f_reader,line);
    // printf("Line 3: %s\n",line.c_str()); TODO REMOVE DEBUG LINE

    // TODO : make more error proof for other .ply files
    std::string num_triangles;
    for (int i = 13; i < 18; i++) {
      num_triangles = num_triangles + line[i];
    }
    triangle_count = 6320;

    //triangle_count = std::stoi(num_triangles);
    cout << "triangle_count : " << triangle_count << "\n";

    // create arrays to store vertices and triangles
    // Vertex vertices[3644]; // get numbers from text file
    // TriangleIndex triangles[6320]; // get numbers from text file

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
        //cout << value << "\n";
      }
      // store vertex data in array
      //vertices[i] = Vertex(x, y, z);

      PolyMesh::vertex[i] = Vertex(x, y, z);

      // apply transform
      transform->apply(PolyMesh::vertex[i]);

      
    }
    cout << "All vertices processed. \n";


    // QUESTION : do the numbers in the triangle lines correspond to lines in the file or index of coordinate in list of that data alone?
    // my guess answer : the numbers stored in the triangles area of the .ply file are the list index values of the vertex 

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
          PolyMesh::triangle[i][axis_count] = number;
        }
        axis_count++;
        //cout << value << "\n";
      }
    }
    cout << "All faces processed. \n";
  }

  f_reader.close();
}
