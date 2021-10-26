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
 * Compile the code using g++ -o lab2executable main_lab3.cpp framebuffer.cpp linedrawer.cpp polymesh.cpp -lm
 *
 * Execute the code using ./lab2executable
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

#define screen_width 2048
#define screen_height 2048


int main(int argc, char *argv[])
{

  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(screen_width,screen_height);

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  Transform *transform = new Transform(1.0f, 0.0f, 0.0f, 0.0f,0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 1.0f, 7.0f,0.0f,0.0f,0.0f,1.0f);

  // Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot.ply", transform);

  // For each triangle in the model,
  for (int i = 0; i< pm->triangle_count; i += 1)
  {
    // The following lines project the point onto the 2D image from 3D space.
    float x0 = (pm->vertex[pm->triangle[i][0]].x/pm->vertex[pm->triangle[i][0]].z)*2000.0 + 1024.0;
    float y0 = (pm->vertex[pm->triangle[i][0]].y/pm->vertex[pm->triangle[i][0]].z)*-2000.0 + 1024.0;
    float x1 = (pm->vertex[pm->triangle[i][1]].x/pm->vertex[pm->triangle[i][1]].z)*2000.0 + 1024.0;
    float y1 = (pm->vertex[pm->triangle[i][1]].y/pm->vertex[pm->triangle[i][1]].z)*-2000.0 + 1024.0;
    float x2 = (pm->vertex[pm->triangle[i][2]].x/pm->vertex[pm->triangle[i][2]].z)*2000.0 + 1024.0;
    float y2 = (pm->vertex[pm->triangle[i][2]].y/pm->vertex[pm->triangle[i][2]].z)*-2000.0 + 1024.0;
    
    // then draw the three edges.
    draw_line(fb, (int)x0, (int)y0, (int)x1, (int)y1);
    draw_line(fb, (int)x1, (int)y1, (int)x2, (int)y2);
    draw_line(fb, (int)x2, (int)y2, (int)x0, (int)y0);

  }

  Sphere ball = Sphere(Vertex (55,55,55), 25);

  
  
  Ray rays = Ray[][];
  rays = new Ray[screen_width][screen_height];

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
      int t = INFINITY;

      Hit new_t = Hit();

      ball->intersection(ray, new_t);


    }
  }


  // code from slides
  for each (x, y) in (screen_width, screen_height){
	  ray = camera_model(x, y); t = infinity; closest = null;

    for each obj in scene {
      if (obj.intersect(ray, &new_t) == true){
        if (new_t < t) { t = new_t; closest = obj;}
      }
    }
	  pos = ray.position(t); draw(x, y, closest.colour(pos));
  }
  // code from slides end


  // Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");

  

  return 0;
  
}
