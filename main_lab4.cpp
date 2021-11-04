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
 * Compile the code using g++ -o lab4executable main_lab4.cpp framebuffer.cpp polymesh.cpp -lm
 *
 * Execute the code using ./lab4executable
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

#define screen_width 212
#define screen_height 212

int main(int argc, char *argv[])
{
  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(screen_width,screen_height);

  // original transformation matrix
  /*Transform *transform = new Transform(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f, 
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
    );*/

  // original transformation matrix
  Transform *transform = new Transform(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, -1.0f, 
    0.0f, 0.0f, 1.0f, 7.0f,
    0.0f, 0.0f, 0.0f, 1.0f
    );

    // new matrix from soren in discord

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  /*Transform *transform = new Transform(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, -1.0f, 3.0f, 
    0.0f, 1.0f, 0.0f, 7.0f,
    0.0f, 0.0f, 0.0f, 1.0f
    );*/

  // Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot.ply", transform);

  // map each pixel to a value between -1 and 1
  float xInt = 2.0f/(float) screen_width;
  float yInt = 2.0f/(float) screen_height;

  // for each section in the mapped size image
  for (float ray_x = -1.0; ray_x < 1.0; ray_x+=xInt){
    for (float ray_y = -1.0; ray_y < 1.0; ray_y+=yInt){
      
      Hit hit;
      hit.flag = false;
      hit.t = 99999999;
      float ray_z = 1;
      
      // generate the shooting ray
      Ray ray (Vertex (0,0,0), Vector (ray_x+0.00222,ray_y+0.00222,ray_z+0.00222));
      ray.direction.normalise();

      pm->intersection(ray, hit);

      int w = (ray_x+1)*(screen_width/2);
      int h = screen_height-1-(ray_y+1)*(screen_height/2);

      if (hit.flag==true){
        fb->plotDepth(w,h,hit.t);
        // calculate colour on object based on normal
        float* colour = pm->colour_hit(hit);
        //printf("red: %f", colour[0]);
        //printf(" green: %f", colour[1]);
        //printf(" blue: %f", colour[2]);
        //printf("\n\n");
        fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
      } else {
        fb->plotDepth(w,h,0);
        float* colour = pm->colour_no_hit(ray);
        fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
      }

    }
  }

  // Output the framebuffer.
  fb->writeDepthFile((char *)"test depth.ppm");
  fb->writeRGBFile((char *)"test colour.ppm");

  return 0;
  
}
