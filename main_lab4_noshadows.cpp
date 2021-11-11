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
#include <sstream>
#include <math.h>
#include <float.h>

#define screen_width 150
#define screen_height 150

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
  /*Transform *transform = new Transform(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, -1.0f, 
    0.0f, 0.0f, 1.0f, 7.0f,
    0.0f, 0.0f, 0.0f, 1.0f
    );*/

  // new matrix from soren in discord
  // [[1, 0, 0, 0], [0, -0.00000004371139, -1, 0], [0, 1, -0.00000004371139, 0], [0, -1.5, 7, 1]]
  /*Transform *transform = new Transform(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, -0.00000004371139f, -1.0f, 0.0f, 
    0.0f, 1.0f, -0.00000004371139f, 0.0f,
    0.0f, -1.5f, 7.0f, 1.0f
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

  // set surface coefficients for lighting for each colour
  pm->set_coeffs(0, 0.8, 0.8, 0, 0.8, 0.8, 0.4, 0.4, 0.4);

  // create light - set ambient light intensity and diffuse intensity of light
  Lighting light (0.5,0.8, Vertex (-2,-1,1));

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

      pm->intersection(ray, hit, 1);
      
      int w = (ray_x+1)*(screen_width/2);
      //int h = screen_height-1-(ray_y+1)*(screen_height/2);
      int h = (ray_y+1)*(screen_height/2);
      if (hit.flag==true){

        // TODO calculate lighting here ?
        float* colour = pm->calculate_lighting(hit, light);
        //printf("DEBUG1");

        fb->plotDepth(w,h,hit.t);
        // calculate colour on object based on normal
        //float* colour = pm->colour_hit(hit);

        // TODO : remove print statements
        /*
        printf("red: %f", colour[0]);
        printf(" green: %f", colour[1]);
        printf(" blue: %f", colour[2]);
        printf("\n\n");*/


        fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
      } else {
        fb->plotDepth(w,h,0);
        float* colour = pm->colour_no_hit(ray);
        //printf("DEBUG2");

        //float* colour = pm->calculate_lighting(hit, light);

        fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
      }

    }
  }

  // Output the framebuffer.
  fb->writeDepthFile((char *)"test depth.ppm");
  fb->writeRGBFile((char *)"test colour.ppm");

  return 0;
  
}