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
 * Compile the code using g++ -o lab5executable main_lab5.cpp framebuffer.cpp polymesh.cpp -lm
 *
 * Execute the code using ./lab5executable
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
#include "scene.h"
#include "colour.h"

#include <iostream>
#include <sstream>
#include <math.h>
#include <float.h>

#define screen_width 500
#define screen_height 500

int main(int argc, char *argv[])
{
  // create a scene that contains the ambient light
  // Scene *picture = new Scene(0.5);
  Scene picture (0.5);

  // create lights
  Lighting light (0.8, Vertex (-4,-1,1));

  

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(screen_width,screen_height);

  // identity transformation matrix
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

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  Transform *transform = new Transform(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, -1.0f, 2.0f, 
    0.0f, 1.0f, 0.0f, 5.5f,
    0.0f, 0.0f, 0.0f, 1.0f
    );

  // Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot_small.ply", transform, 0);
  //PolyMesh *pm = new PolyMesh((char *)"teapot_big.ply", transform, 1); 

  // create sphere
  Sphere ball (Vertex(2,2,5), 2.5);

  // create list of objects
  //Object list = [pm, ball];

  Object ball = ball;
  Object teapot = *pm;
  ball.next = &teapot;

  // set surface coefficients for lighting for each component and colour
  pm->set_coeffs(0, 0.8, 0.8, 0, 0.8, 0.8, 0.4, 0.4, 0.4);  

  // map each pixel to a value between -1 and 1
  float xInt = 2.0f/(float) screen_width;
  float yInt = 2.0f/(float) screen_height;

  // for each section in the mapped size image
  for (float ray_x = -1.0; ray_x < 1.0; ray_x+=xInt){
    for (float ray_y = -1.0; ray_y < 1.0; ray_y+=yInt){

      int w = (ray_x+1)*(screen_width/2);
      int h = (ray_y+1)*(screen_height/2);
      
      Hit shooting_hit;
      shooting_hit.flag = false;
      shooting_hit.t = 99999999;
      float ray_z = 1;
      
      // generate the shooting ray
      Ray shooting_ray (Vertex (0,0,0), Vector (ray_x+0.00222,ray_y+0.00222,ray_z+0.00222));
      shooting_ray.direction.normalise();

      Colour colour = picture.raytracer(shooting_ray, shooting_hit);







      pm->intersection(shooting_ray, shooting_hit);
      ball.intersection(shooting_ray, shooting_hit);
      // TODO : reverse y direction?
      //int h = screen_height-1-(ray_y+1)*(screen_height/2);

      // determine which lighting calculation is needed
      if (shooting_hit.flag==true){
        // calculate if shadows here
        Hit shadow_hit;
        shadow_hit.flag = false;
        shadow_hit.t = 99999999;
        
        Vertex shadow_point = Vertex (shooting_hit.position.x-0.00222, shooting_hit.position.y-0.00222, shooting_hit.position.z-0.00222);
        Vector shadow_dir = pm->getDirection(shadow_point, light.position);
        Ray shadow_ray (shadow_point, shadow_dir);

        pm->intersection(shadow_ray, shadow_hit);
        ball.intersection(shadow_ray, shadow_hit);

        if (shadow_hit.flag==true){
          // only ambient lighting
          float* colour = pm->calculate_lighting(shooting_hit, picture.ambient_intensity, light, 1);
          fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
        } else {
          // ambient, diffuse and specular lighting
          float* colour = pm->calculate_lighting(shooting_hit, picture.ambient_intensity, light, 2);
          fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
        }
      } else {
        // background - currently just black
        float* colour = pm->calculate_lighting(shooting_hit, picture.ambient_intensity, light, 3);
        fb->plotPixel(w,h,colour[0],colour[1],colour[2]);
      }

      // plot pixel here

      // plot depth
      if (shooting_hit.flag==true){
        fb->plotDepth(w,h,shooting_hit.t);
      } else {
        fb->plotDepth(w,h,0);
      }

    }
  }

  // Output the framebuffer.
  fb->writeDepthFile((char *)"test depth.ppm");
  fb->writeRGBFile((char *)"test colour.ppm");

  return 0;
  
}
