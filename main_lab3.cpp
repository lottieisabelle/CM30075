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
 * Compile the code using g++ -o lab3executable main_lab3.cpp framebuffer.cpp polymesh.cpp -lm
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

#define screen_width 2048
#define screen_height 2048

Vector getDirection(Vertex a, Vertex b){
  return Vector ((b.x-a.x),(b.y-a.y),(b.z-a.z));
}

int main(int argc, char *argv[])
{

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(screen_width,screen_height);

  /* original transformation matrix
  Transform *transform = new Transform(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, -1.0f, 
    0.0f, 0.0f, 1.0f, 7.0f,
    0.0f, 0.0f, 0.0f, 1.0f
    );*/

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  Transform *transform = new Transform(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, -1.0f, 3.0f, 
    0.0f, 1.0f, 0.0f, 7.0f,
    0.0f, 0.0f, 0.0f, 1.0f
    );

  // Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot.ply", transform);

  Vertex camera (0,0,0);

  // map each pixel to a value between -1 and 1
  float xInt = 2.0f/(float) screen_width;
  float yInt = 2.0f/(float) screen_height;

  // for each section in the mapped size image
  for (float ray_x = -1.0; ray_x < 1.0; ray_x+=xInt){
    for (float ray_y = -1.0; ray_y < 1.0; ray_y+=yInt){
      //printf("ray x: %f\n",ray_x);
      //printf("ray_y: %f\n\n",ray_y);

      bool hit = false;
      float ray_z = 1;
      float closest_plot = 99999999;
      //printf("closest plot: %f\n",closest_plot);

      // for each triangle
      for(int i = 0; i < pm->triangle_count; i+=1){
        // get the vertex data
        Vertex a = pm->vertex[pm->triangle[i][0]];
        Vertex b = pm->vertex[pm->triangle[i][1]];
        Vertex c = pm->vertex[pm->triangle[i][2]];

        // create vectors
        Vector ab = getDirection(a,b);
        Vector ac = getDirection(a,c);
        Vector bc = getDirection(b,c);
        Vector ca = getDirection(c,a);

        // find the normal to the plane abc
        Vector N;
        // get cross product of ab and ac and store in n
        ab.cross(ac,N);
        N.normalise();

        // generate the shooting ray
        Vector D (ray_x+0.00222,ray_y+0.00222,ray_z+0.00222);
        Ray ray (camera, D);
        ray.direction.normalise();

        // get direction vector between the camera (0,0,0) and point on plane e.g. a
        Vector dir = getDirection(camera,a);
        //dir.normalise();

        // if you multiply ray.direction by d, then you get the point on the plane
        float d = dir.dot(N)/ray.direction.dot(N);
        if (d < 0){
          // not intersecting
          continue;
        } else {
          hit = true;
        }

        // point on plane where shooting ray intersects
        Vertex P (camera.x + ray.direction.x*d, camera.y + ray.direction.y*d, camera.z + ray.direction.z*d);

        // now need to know if the point P is inside the triangle on the plane
        // get vectors PA, PB and PC
        Vector PA = getDirection(P,a);
        Vector PB = getDirection(P,b);
        Vector PC = getDirection(P,c);

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
          if(d < closest_plot){
            closest_plot = d;
          }
        }  
        /*
        int w = (ray_x+1)*(screen_width/2);
        int h = (ray_y+1)*(screen_height/2);
        fb->plotDepth(w,h,closest_plot);
        */

      }

      int w = (ray_x+1)*(screen_width/2);
      int h = (ray_y+1)*(screen_height/2);

      if (closest_plot < 99999999){
        fb->plotDepth(w,h,closest_plot);
      } else {
        fb->plotDepth(w,h,0);
      }
      
    }
  }

  //printf("writing to file");

  // Output the framebuffer.
  //fb->writeRGBFile((char *)"test.ppm");
  fb->writeDepthFile((char *)"test.ppm");


  return 0;
  
}
