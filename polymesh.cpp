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
#include "object.h"

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

PolyMesh::PolyMesh(char *file, int x)
{
  Transform *transform = new Transform();

  this->do_construct(file, transform, x);
}

PolyMesh::PolyMesh(char *file, Transform *transform, int x)
{
  this->do_construct(file, transform, x);
}

void PolyMesh::do_construct(char *file, Transform *transform, int flag)
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
          number -= flag;
          
          PolyMesh::triangle[i][axis_count-1] = number;
        }
        axis_count++;
      }
    }
    cout << "All faces processed. \n";
  }

  f_reader.close();
}

void PolyMesh::intersection(Ray ray, Hit &hit)
{
  // for each triangle
  for(int i = 0; i < this->triangle_count; i+=1){
    // get the vertex data
    Vertex a = this->vertex[this->triangle[i][0]];
    Vertex b = this->vertex[this->triangle[i][1]];
    Vertex c = this->vertex[this->triangle[i][2]];

    // create vectors

    Vector ab = a.getDirection(b);
    Vector ab = a.getDirection(c);
    Vector bc = b.getDirection(c);
    Vector ca = c.getDirection(a);

    // find the normal to the plane abc
    Vector N;
    // get cross product of ab and ac and store in n
    ab.cross(ac,N);
    N.normalise();

    if(ray.direction.dot(N) == 0){
      continue;
    }
    
    // get direction vector between the camera (0,0,0) and point on plane e.g. a
    Vector dir = ray.position.getDirection(a);

    float d = dir.dot(N)/ray.direction.dot(N);
    
    if (d < 0){
      // not intersecting
      continue;
    } 
  
    // if you multiply ray.direction by d, then you get the point on the plane
    // point on plane where shooting ray intersects
    Vertex P (ray.position.x + ray.direction.x*d, ray.position.y + ray.direction.y*d, ray.position.z + ray.direction.z*d);

    // now need to know if the point P is inside the triangle on the plane
    // get vectors PA, PB and PC
    Vector PA = P.getDirection(a);
    Vector PB = P.getDirection(b);
    Vector PC = P.getDirection(c);

    // calculate cross products to get normal at each vertex of triangle
    Vector a_normal;
    PA.cross(ab,a_normal);
    Vector b_normal;
    PB.cross(bc,b_normal);
    Vector c_normal;
    PC.cross(ca,c_normal);

    // if point is inside shape, all normals will be pointing in the same direction
    // if dot product between 2 vectors is greater than zero, they are pointing in the same direction 
    if ((a_normal.dot(b_normal) > 0) && (b_normal.dot(c_normal) > 0)){
      if(d < hit.t){
        hit.t = d;
        hit.flag = true;
        hit.normal = N;
        hit.position = P;
        hit.what = this;
      }
    }  

  }

}

float* PolyMesh::colour_hit(Hit &hit)
{
  
  float red = 0.5 * (hit.normal.x+1);
  float green = 0.5 * (hit.normal.y+1);
  float blue = 0.5 * (-hit.normal.z+1);

  //printf("red: %f", red);
  //printf(" green: %f", green);
  //printf(" blue: %f", blue);

  float* colour {new float[3] {red,green,blue}};

  return colour;
}

float* PolyMesh::colour_no_hit(Ray ray)
{
  // based on y direction of ray
  // function takes ray not hit
  // create value called t which is the 0.5 * (ray.direction.y + 1)
  // for each colour
  // colour is then (1- t) * (colour values of white) + t * (colour values of the other colour in the gradient)
  float t = 0.5 * (ray.direction.y-0.00222 +1);

  float red = (1-t) * 1 + t * 0.5;
  float green = (1-t) * 1 + t * 0.7;
  float blue = (1-t) * 1 + t * 1.0;

  //printf("red: %f", red);
  //printf(" green: %f", green);
  //printf(" blue: %f", blue);

  float* colour {new float[3] {red,green,blue}};

  return colour;
}

float* PolyMesh::calculate_lighting(Hit &hit, float Ia, Lighting light, int flag)
{
  float red;
  float green;
  float blue;
  if (flag == 1){
    // only ambient light due to shadows
    red = Ia*ambient[0];
    green = Ia*ambient[1];
    blue = Ia*ambient[2];

  } else if (flag == 2){
    // ambient, diffuse and specular
    
    Vector L = hit.position.getDirection(light.position);
    L.normalise();
    Vector I = light.position.getDirection(hit.position);
    I.normalise();
    Vector R;
    hit.normal.reflection(I,R);
    R.normalise();
    Vector V = hit.position.getDirection(Vertex (0,0,0));
    V.normalise();
    int n = 40;

    red = Ia*ambient[0] + light.diffuse_intensity *  ( diffuse[0]*(hit.normal.dot(L))  +  specular[0]*  pow(R.dot(V),n));
    green = Ia*ambient[1] + light.diffuse_intensity *  ( diffuse[1]*(hit.normal.dot(L))  +  specular[1]*  pow(R.dot(V),n));
    blue = Ia*ambient[2] + light.diffuse_intensity *  ( diffuse[2]*(hit.normal.dot(L))  +  specular[2]*  pow(R.dot(V),n));

  } else if (flag == 3){
    // background (currently black 11/11 10:54)
    red = 0;
    green = 0;
    blue = 0;
  }
  
  /*
  printf("red: %f", red);
  printf(" green: %f", green);
  printf(" blue: %f", blue);
  printf("\n\n");*/

  float* colour {new float[3] {red,green,blue}};

  return colour;
}

void PolyMesh::set_coeffs(float ar, float ag, float ab, float dr, float dg, float db, float sr, float sg, float sb)
{
  PolyMesh::ambient = new float[3];
  PolyMesh::diffuse = new float[3];
  PolyMesh::specular = new float[3];

  ambient[0] = ar;
  ambient[1] = ag;
  ambient[2] = ab; 

  diffuse[0] = dr;
  diffuse[1] = dg;
  diffuse[2] = db;

  specular[0] = sr;
  specular[1] = sg;
  specular[2] = sb; 

}