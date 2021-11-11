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
  printf("hello\n");
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

  // set surface coefficients for lighting for each component
  pm->set_coeffs(0, 0.8, 0.8, 0, 0.8, 0.8, 0.4, 0.4, 0.4);

  // create light - set ambient light intensity and diffuse intensity of light
  //Lighting light (1,1, Vertex (2,1,-1));
  Lighting light (0.5,0.8, Vertex (-2,-1,1));

  // map each pixel to a value between -1 and 1
  float xInt = 2.0f/(float) screen_width;
  float yInt = 2.0f/(float) screen_height;

  // for each section in the mapped size image
  for (float ray_x = -1.0; ray_x < 1.0; ray_x+=xInt){
    for (float ray_y = -1.0; ray_y < 1.0; ray_y+=yInt){
      
      Hit shooting_hit;
      shooting_hit.flag = false;
      shooting_hit.t = 99999999;
      float ray_z = 1;
      
      // generate the shooting ray
      Ray shooting_ray (Vertex (0,0,0), Vector (ray_x+0.00222,ray_y+0.00222,ray_z+0.00222));
      shooting_ray.direction.normalise();
      pm->intersection(shooting_ray, shooting_hit, 1);

      // calculate shadows here
      // create shadow data storage
      Hit shadow_hit;
      shadow_hit.flag = false;
      shadow_hit.t = 99999999;
      Vertex shadow_point = Vertex (ray_x+0.00222,ray_y+0.00222,ray_z+0.00222);
      Vector shadow_dir = pm->getDirection(shadow_point, light.position);
      Ray shadow_ray (shadow_point, shadow_dir);

      pm->intersection(shadow_ray, shadow_hit, 2);
      if(shadow_hit.flag ==true){
        //printf("true");
      }

      int w = (ray_x+1)*(screen_width/2);
      int h = (ray_y+1)*(screen_height/2);
      // TODO : reverse y direction?
      //int h = screen_height-1-(ray_y+1)*(screen_height/2);

      float* colour;
      // calculate lighting
      if (shooting_hit.flag==true){
        fb->plotDepth(w,h,shooting_hit.t);

        // calculate lighting here
        float* colour = pm->calculate_lighting(shooting_hit, light);
        
        // calculate colour on object based on normal
        //float* colour = pm->colour_hit(hit);

        fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
      } else {
        fb->plotDepth(w,h,0);
        float* colour = pm->colour_no_hit(shooting_ray);
        

        //float* colour = pm->calculate_lighting(hit, light);

        fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
      }
      
      // calculate shadows
      if (shadow_hit.flag==true){
        //printf("true\n");
        float* colour = pm->calculate_lighting(shooting_hit, light);
        //printf("colour: %f, %f, %f\n", colour[0],colour[1],colour[2]);
        // reduce all colour values by 0.6, set to 0 if below 0
        float temp;
        for (int i = 1; i < 3; i++){
          temp = colour[i];
          temp = temp - 0.5;
          if(temp<0){
            temp = 0;
          }
          colour[i] = temp;
        }
        //printf("shadow colour: %.4f, %.4f, %.4f\n", colour[0],colour[1],colour[2]);
        fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
      }
      
      


    }
  }

  // Output the framebuffer.
  fb->writeDepthFile((char *)"test depth.ppm");
  fb->writeRGBFile((char *)"test colour.ppm");

  return 0;
  
}
