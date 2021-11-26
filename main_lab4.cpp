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
 * Compile the code using g++ -o lab4executable main_lab4.cpp framebuffer.cpp polymesh.cpp sphere.cpp phong.cpp directional_light.cpp -lm
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
#include "ray.h"
#include "hit.h"
#include "polymesh.h"
#include "sphere.h"
#include "light.h"
#include "directional_light.h"
#include "material.h"
#include "phong.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

void object_test(Ray ray, Object *objects, Hit &best_hit)
{
  Object *obj = objects;

  best_hit.flag = false;


  while(obj != 0)
  {
    Hit obj_hit;
    obj_hit.flag=false;
	  
    obj->intersection(ray, obj_hit);

    
    if (obj_hit.flag)
    {
      if (obj_hit.t > 0.0f)
      {
        if (best_hit.flag == false)
	      {
	        best_hit = obj_hit;
	      } else if (obj_hit.t < best_hit.t)
	      {
	        best_hit = obj_hit;
	      }
      }
    }
    
    obj = obj->next;
  }


  return;
}

void raytrace(Ray ray, Object *objects, Light *lights, Colour &colour, float &depth, int d)
{
  if(d <= 0){
    return;
  }
  // first step, find the closest primitive

  Hit shadow_hit;
  Hit best_hit;
  object_test(ray, objects, best_hit);

  
  // if we found a primitive then compute the colour we should see
  if(best_hit.flag)
  {
    best_hit.what->material->compute_base_colour(colour);
    depth = best_hit.t;
    Light *light = lights;

    while (light != (Light *)0)
    {
      Vector viewer;
      Vector ldir;

      viewer.x = -best_hit.position.x;
      viewer.y = -best_hit.position.y;
      viewer.z = -best_hit.position.z;
      viewer.normalise();

      bool lit;
      lit = light->get_direction(best_hit.position, ldir);

      if(ldir.dot(best_hit.normal)>0)
      {
	      lit=false;//light is facing wrong way.
      }

      if(lit)
      {
      
        Ray shadow_ray;

        shadow_ray.direction.x = -ldir.x;
        shadow_ray.direction.y = -ldir.y;
        shadow_ray.direction.z = -ldir.z;
        shadow_ray.position.x = best_hit.position.x + (0.0001f * shadow_ray.direction.x);
        shadow_ray.position.y = best_hit.position.y + (0.0001f * shadow_ray.direction.y);
        shadow_ray.position.z = best_hit.position.z + (0.0001f * shadow_ray.direction.z);

        object_test(shadow_ray, objects, shadow_hit);

		    if(shadow_hit.flag==true)
		    {
			    if (shadow_hit.t < 1000000000.0f)
			    {
				    lit = false; //there's a shadow so no lighting, if realistically close
			    }
        }
      }

      if (lit)
      {
        Colour intensity;
		    Colour scaling;

		    light->get_intensity(best_hit.position, scaling);

		    best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);

        intensity.scale(scaling);

		    colour.add(intensity);
      }

      light = light->next;
    }

    // compute reflection ray if material supports it.
    if(best_hit.what->material->k_reflection > 0.0f)
    {
      Vector r_dir;
      best_hit.normal.reflection(ray.direction, r_dir);
      r_dir.normalise();

      Vertex r_pos;
      r_pos.x = best_hit.position.x + 0.00222 * r_dir.x;
      r_pos.y = best_hit.position.y + 0.00222 * r_dir.y;
      r_pos.z = best_hit.position.z + 0.00222 * r_dir.z;

      Ray r_ray (r_pos, r_dir);

      float kr = best_hit.what->material->k_reflection;

      Colour kr_col;
      kr_col.r = kr;
      kr_col.g = kr;
      kr_col.b = kr;

      Colour col;
      raytrace(r_ray, objects, lights, col, depth, d-1);

      col.scale(kr_col);
      colour.add(col);
      
    }

    // compute refraction ray if material supports it.
    if(best_hit.what->material->index_refraction > 0.0f)
    {
      float ior = best_hit.what->material->index_refraction;
      // tray.dir = refraction(ray.dir, hit.normal, hit.ior);
      
      // cos θi = N.I
      // I = incident ray direction vector
      // N = normal to surface direction vector
      float cos_i = best_hit.normal.dot(ray.direction);

      // cos θt = sqrt(1 – (1/η2) * (1 - cos2 θi) )
      float n2 = ior*ior;
      float cos_i2 = cos_i*cos_i;

      float cos_t = sqrt(1 - (1/n2) * (1 - cos_i2) );

      // T = 1/η * I – (cos θt – (1/η)* cos θi ) * N
      Vector T_dir = 1/ior * ray.direction - (cos_t - (1/ior) * cos_i) * best_hit.normal;
      
      // tray.pos = hit.pos + small_e * tray.dir;
      Vertex T_pos;
      T_pos.x = best_hit.position.x + 0.00222 * T_dir.x;
      T_pos.y = best_hit.position.y + 0.00222 * T_dir.y;
      T_pos.z = best_hit.position.z + 0.00222 * T_dir.z;
      
      Ray T_ray (T_pos, T_dir);

      // col += hit.object.kt * raytrace(tray, depth-1);
      float kt = best_hit.what->material->k_refraction;

      Colour kt_col;
      kt_col.r = kt;
      kt_col.g = kt;
      kt_col.b = kt;

      Colour col;
      raytrace(T_ray, objects, lights, col, depth, d-1);

      col.scale(kt_col);
      colour.add(col);

    }

  } else
  {
    depth = 7.0f;
    colour.r = 0.0f;
    colour.g = 0.0f;
    colour.b = 0.0f;
  }	
}

int main(int argc, char *argv[])
{
  int width = 256;
  int height = 256;
  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
				       0.0f, 0.0f, 1.0f, -2.7f,
				       0.0f, 1.0f, 0.0f, 5.0f,
				       0.0f, 0.0f, 0.0f, 1.0f);

  //  Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot_smaller.ply", transform);

  // create material properties for teapot
  Phong bp1;

	bp1.ambient.r = 0.0f;
	bp1.ambient.g = 0.4f;
	bp1.ambient.b = 0.4f;
	bp1.diffuse.r = 0.0f;
	bp1.diffuse.g = 0.4f;
	bp1.diffuse.b = 0.4f;
	bp1.specular.r = 0.4f;
	bp1.specular.g = 0.4f;
	bp1.specular.b = 0.4f;
	bp1.power = 40.0f;

	pm->material = &bp1;

  pm->material->k_reflection = 0.4f;
  // no refraction
  pm->material->index_refraction = 0.0f;
  pm->material->k_refraction = 0.0f;

  // create spheres and material properties of spheres
  Vertex v;
  v.x = 2.0f;
  v.y = 2.0f;
  v.z = 5.0f;
  
  Sphere *sphere = new Sphere(v, 0.75f);
  Phong bp2;

  bp2.ambient.r = 0.0f;
	bp2.ambient.g = 0.2f;
	bp2.ambient.b = 0.0f;
	bp2.diffuse.r = 0.0f;
	bp2.diffuse.g = 0.4f;
	bp2.diffuse.b = 0.0f;
	bp2.specular.r = 0.4f;
	bp2.specular.g = 0.4f;
	bp2.specular.b = 0.4f;
	bp2.power = 40.0f;

	sphere->material = &bp2;

  sphere->material->k_reflection = 0.8f;
  sphere->material->k_refraction = 0.8f;
  sphere->material->index_refraction = 1.33f; // water

  Vertex v2;
  v2.x = -1.0f;
  v2.y = 1.0f;
  v2.z = 3.0f;
  
  Sphere *sphere2 = new Sphere(v2,0.5f);

  Phong bp3;
  // create material properties for spheres

  bp3.ambient.r = 0.0f;
	bp3.ambient.g = 0.2f;
	bp3.ambient.b = 0.0f;
	bp3.diffuse.r = 0.0f;
	bp3.diffuse.g = 0.4f;
	bp3.diffuse.b = 0.0f;
	bp3.specular.r = 0.4f;
	bp3.specular.g = 0.4f;
	bp3.specular.b = 0.4f;
	bp3.power = 40.0f;

 	sphere2->material = &bp3;

  sphere2->material->k_reflection = 0.8f;
  sphere2->material->k_refraction = 0.8f;
  // glass
  sphere2->material->index_refraction = 1.52f;
    
  // link objects
  pm->next = sphere;
  sphere->next = sphere2;
  
  // generate shooting ray from camera point
  Ray ray;

  ray.position.x = 0.0001f;
  ray.position.y = 0.0f;
  ray.position.z = 0.0f;

  // create directional light
  DirectionalLight *dl = new DirectionalLight(Vector(1.01f, -1.0f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f));

  for (int y = 0; y < height; y += 1)
  {
    for (int x = 0; x < width; x += 1)
    {
      float fx = (float)x/(float)width;
      float fy = (float)y/(float)height;

      Vector direction;

      ray.direction.x = (fx-0.5f);
      ray.direction.y = (0.5f-fy);
      ray.direction.z = 0.5f;
      ray.direction.normalise();

      Colour colour;
      float depth;

      // reflection recursion depth
      int d = 3;

      raytrace(ray, pm, dl, colour, depth, d);

      fb->plotPixel(x, y, colour.r, colour.g, colour.b);
      //fb->plotDepth(x,y, depth);
    }

    cerr << "*" << flush;
  }
  
  // Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");
  //  fb->writeDepthFile((char *)"depth.ppm");
  return 0;
  
}
