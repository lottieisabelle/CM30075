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
    printf("Header: %s\n",line.c_str());

    // parse the number of vertices line
    std::getline(f_reader,line);
    printf("Line 2: %s\n",line.c_str());

    for (int i = 15; i < 19; i++) {
      

      cout << i << "\n";
    }
    
  }

  f_reader.close();
}
