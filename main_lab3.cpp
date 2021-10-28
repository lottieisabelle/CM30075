/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the entry point function for the program you need to create for lab two.
 * You should not need to modify this code.
 * It creates a framebuffer, loads an triangle mesh object, calls the drawing function to render the object and then outputs the framebuffer as a ppm file.
 *
 * On linux.bath.ac.uk:
 *
 * Compile the code using g++ -o lab3executable main_lab3.cpp framebuffer.cpp linedrawer.cpp polymesh.cpp sphere.cpp -lm
 *
 * Execute the code using ./lab3executable
 *
 * This will produce an image file called test.ppm. You can convert this a png file for viewing using
 *
 * pbmropng test.ppm > test.png
 *
 * You are expected to fill in the missing code in polymesh.cpp.
 */

#include "framebuffer.h"
#include "linedrawer.h"
#include "polymesh.h"
#include "ray.h"
#include "sphere.h"

#include <iostream>
#include <math.h>
#include <float.h>

#define screen_width 64
#define screen_height 64

Vector getDirection(Vertex a, Vertex b){
  return Vector ((b.x-a.x),(b.y-a.y),(b.z-a.z));
}

int main(int argc, char *argv[])
{

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(screen_width,screen_height);

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  Transform *transform = new Transform(1.0f, 0.0f, 0.0f, 0.0f,0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 1.0f, 7.0f,0.0f,0.0f,0.0f,1.0f);

  // Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot.ply", transform);


  Vertex camera (0,0,0);

  // map each pixel to a value between -1 and 1
  float xInt = 2.0f/(float) screen_width;
  float yInt = 2.0f/(float) screen_height;

  // for each section in the mapped size image
  for (float ray_x = -1.0; ray_x < 1.0; ray_x+=xInt){
    for (float ray_y = -1.0; ray_x < 1.0; ray_y+=yInt){
      float ray_z = 1;
      float closest_plot = FLT_MAX;

      // for each triangle
      for(int i = 0; i < pm->triangle_count; i+=1){
        // get the vertex data
        Vertex a = pm->vertex[pm->triangle[i][0]];
        Vertex b = pm->vertex[pm->triangle[i][1]];
        Vertex c = pm->vertex[pm->triangle[i][2]];

        // create vectors
        Vector ab ((b.x-a.x),(b.y-a.y),(b.z-a.z));
        Vector ac ((c.x-a.x),(c.y-a.y),(c.z-a.z));
        Vector bc ((c.x-b.x),(c.y-b.y),(c.z-b.z));

        // find the normal to the plane abc
        Vector N;
        // get cross product of ab and ac and store in n
        ab.cross(ac,N);
        N.normalise();

        // generate the shooting ray
        Vector D (ray_x,ray_y,ray_z);
        Ray ray (camera, D);
        ray.direction.normalise();

        // get direction vector between the camera (0,0,0) and point on plane e.g. a
        Vector dir (a.x-camera.x,a.y-camera.y,a.z-camera.z);

        // if you multiply ray.direction by d, then you get the point on the plane
        float d = dir.dot(N)/ray.direction.dot(N);

        Vertex P (camera.x + ray.direction.x*d, camera.y + ray.direction.y*d,camera.z+ray.direction.z*d);

        // now need to know if the point P is inside the triangle on the plane
        // get vectors PA, PB and PC
        Vector PA = getDirection(P,a);
        Vector PB = getDirection(P,b);
        Vector PC = getDirection(P,c);

        // calculate cross products to get normal at each vertex of triangle
        Vector a_normal;
        ab.cross(PA,a_normal);
        Vector b_normal;
        bc.cross(PB,b_normal);
        Vector c_normal;
        ac.cross(PC,c_normal);

        // if point is inside shape, all normals will be pointing in the same direction
        // if dot product between 2 vectors is greater than zero, they are pointing in the same direction
        if (a_normal.dot(b_normal) > 0 && b_normal.dot(c_normal) > 0){
          closest_plot = d;
        }   
      }

      fb->plotDepth(xInt,yInt,closest_plot);
      printf("triangle done");

    }
  }

  printf("writing to file");

  /* // each  verticy is a b c, make a plae, then find intersection from each plane
  // for each triangle in traingale count



  Sphere ball = Sphere(Vertex (5,5,5), 1);

  // for each pixel in image
  for (int wx = 0; wx < screen_width; wx+=1){
    for (int wy = 0; wy < screen_height; wy+=1){
      // generate a ray
      Vertex p = Vertex(0,0,0);
      int x = -1 + wx/1024;
      int y = -1 + wy/1024;
      int z = 1;
      Vector d = Vector(x,y,z);
      d.normalise();

      Ray ray = Ray(p,d);
      float t = 999999999.0;
      Object closest = Object();

      Hit touch = Hit();

      ball.intersection(ray, touch);

      // put loop for checking each triangle plane
      float r = 0;
      float g = 0;
      float b = 0;

      if (touch.t < t){
        t = touch.t;
        closest = ball;
        r = 1;
      }

      //pos = ray.position(t);
      //fb->plotPixel(wx,wy,r,g,b);
      fb->plotDepth(wx,wy,t);

    }
  } */

  // Output the framebuffer.
  //fb->writeRGBFile((char *)"test.ppm");
  fb->writeDepthFile((char *)"test.ppm");


  return 0;
  
}
