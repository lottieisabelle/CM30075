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
 * Compile the code using g++ -o lab4executable main_lab4.cpp framebuffer.cpp polymesh.cpp sphere.cpp phong.cpp directional_light.cpp plane.cpp -lm
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
#include "plane.h"

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

void raytrace(Ray ray, Object *objects, Light *lights, Colour &colour, float &depth, int r_d, int t_d)
{
  if(r_d <= 0){
    return;
  }
  if(t_d <=0){
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
    if(best_hit.what->material->bool_reflection)
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
      raytrace(r_ray, objects, lights, col, depth, r_d-1, t_d);

      col.scale(kr_col);
      colour.add(col);
      
    }

    // compute refraction ray if material supports it.
    if(best_hit.what->material->bool_refraction)
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
      //float x = best_hit.what->material->k_reflection;
      //float kt = 1.0f - x;

      //cout << kt << "\n";
      float kt = best_hit.what->material->k_refraction;

      Colour kt_col;
      kt_col.r = kt;
      kt_col.g = kt;
      kt_col.b = kt;

      Colour col;
      raytrace(T_ray, objects, lights, col, depth, r_d, t_d-1);

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
  int width = 200;
  int height = 200;
  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
				       0.0f, 1.0f, 0.0f, -5.0f,
				       0.0f, 0.0f, 1.0f, 8.0f,
				       0.0f, 0.0f, 0.0f, 1.0f);

  //  Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot.ply", transform);

  // create material properties for teapot
  Phong bp1;
	bp1.ambient.r = 0.0f;
	bp1.ambient.g = 0.9f;
	bp1.ambient.b = 0.9f;
	bp1.diffuse.r = 0.0f;
	bp1.diffuse.g = 0.9f;
	bp1.diffuse.b = 0.9f;
	bp1.specular.r = 0.4f;
	bp1.specular.g = 0.4f;
	bp1.specular.b = 0.4f;
	bp1.power = 40.0f;

	pm->material = &bp1;

  pm->material->bool_reflection = false;
  pm->material->k_reflection = 0.4f;

  pm->material->bool_refraction = false;
  pm->material->index_refraction = 0.0f;
  pm->material->k_refraction = 0.0f;

  // create box for teapot to sit in
  Transform *transform2 = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
				       0.0f, 1.0f, 0.0f, 0.0f,
				       0.0f, 0.0f, 1.0f, 0.0f,
				       0.0f, 0.0f, 0.0f, 1.0f);
  PolyMesh *background_pm = new PolyMesh((char *)"background.ply", transform2);
  PolyMesh *floor_pm = new PolyMesh((char *)"floor.ply", transform2);
  PolyMesh *left_wall = new PolyMesh((char *)"left_wall.ply", transform2);
  PolyMesh *right_wall = new PolyMesh((char *)"right_wall.ply", transform2);
  PolyMesh *ceiling_pm = new PolyMesh((char *)"ceiling.ply", transform2);

  Phong bp4;
  bp4.ambient.r = 1.0f;
	bp4.ambient.g = 0.7f;
	bp4.ambient.b = 1.0f;
	bp4.diffuse.r = 1.0f;
	bp4.diffuse.g = 0.7f;
	bp4.diffuse.b = 1.0f;
	bp4.specular.r = 0.2f;
	bp4.specular.g = 0.2f;
	bp4.specular.b = 0.2f;
	bp4.power = 40.0f;

  Phong bp5;
  bp5.ambient.r = 0.5f;
	bp5.ambient.g = 0.0f;
	bp5.ambient.b = 1.0f;
	bp5.diffuse.r = 0.5f;
	bp5.diffuse.g = 0.0f;
	bp5.diffuse.b = 1.0f;
	bp5.specular.r = 0.2f;
	bp5.specular.g = 0.2f;
	bp5.specular.b = 0.2f;
	bp5.power = 40.0f;

  Phong bp6;
  bp6.ambient.r = 1.0f;
	bp6.ambient.g = 0.5f;
	bp6.ambient.b = 0.0f;
	bp6.diffuse.r = 1.0f;
	bp6.diffuse.g = 0.5f;
	bp6.diffuse.b = 0.0f;
	bp6.specular.r = 0.2f;
	bp6.specular.g = 0.2f;
	bp6.specular.b = 0.2f;
	bp6.power = 40.0f;

  // back wall
  background_pm->material = &bp6;
  background_pm->material->bool_reflection = false;
  background_pm->material->k_reflection = 0.0f;
  
  background_pm->material->bool_refraction = false;
  background_pm->material->k_refraction = 0.0f;
  background_pm->material->index_refraction = 0.0f;

  // floor
  floor_pm->material = &bp4;
  floor_pm->material->bool_reflection = false;
  floor_pm->material->k_reflection = 0.0f;
  
  floor_pm->material->bool_refraction = false;
  floor_pm->material->k_refraction = 0.0f;
  floor_pm->material->index_refraction = 0.0f;

  // left wall
  left_wall->material = &bp5;
  left_wall->material->bool_reflection = false;
  left_wall->material->k_reflection = 0.0f;
  
  left_wall->material->bool_refraction = false;
  left_wall->material->k_refraction = 0.0f;
  left_wall->material->index_refraction = 0.0f;

  // right wall
  right_wall->material = &bp5;
  right_wall->material->bool_reflection = false;
  right_wall->material->k_reflection = 0.0f;
  
  right_wall->material->bool_refraction = false;
  right_wall->material->k_refraction = 0.0f;
  right_wall->material->index_refraction = 0.0f;

  // ceiling
  ceiling_pm->material = &bp4;
  ceiling_pm->material->bool_reflection = false;
  ceiling_pm->material->k_reflection = 0.0f;
  
  ceiling_pm->material->bool_refraction = false;
  ceiling_pm->material->k_refraction = 0.0f;
  ceiling_pm->material->index_refraction = 0.0f;

  // create bubbles
  Vertex v;
  v.x = 3.0f; 
  v.y = 1.0f; 
  v.z = 8.0f; 
  
  Sphere *sphere = new Sphere(v, 0.4);
  Phong bp2;
  bp2.ambient.r = 0.8f;
	bp2.ambient.g = 0.8f;
	bp2.ambient.b = 0.8f;
	bp2.diffuse.r = 0.8f;
	bp2.diffuse.g = 0.8f;
	bp2.diffuse.b = 0.8f;
	bp2.specular.r = 0.4f;
	bp2.specular.g = 0.4f;
	bp2.specular.b = 0.4f;
	bp2.power = 40.0f;

  sphere->material->bool_reflection = true;
  sphere->material->k_reflection = 0.4f;

  sphere->material->bool_refraction = true;
  sphere->material->k_refraction = 0.6f;
  sphere->material->index_refraction = 1.33f; // water


  // create spheres and material properties of spheres
  /*
  Vertex v;
  v.x = 0.0f; // 2
  v.y = 0.1f; // 2
  v.z = 1.5f; // 4
  
  Sphere *sphere = new Sphere(v, 0.4);
  Phong bp2;
  // rgb(127,0,255) purple

  bp2.ambient.r = 0.8f;
	bp2.ambient.g = 0.8f;
	bp2.ambient.b = 0.8f;
	bp2.diffuse.r = 0.8f;
	bp2.diffuse.g = 0.8f;
	bp2.diffuse.b = 0.8f;
	bp2.specular.r = 0.4f;
	bp2.specular.g = 0.4f;
	bp2.specular.b = 0.4f;
	bp2.power = 40.0f;

	sphere->material = &bp2;

  sphere->material->bool_reflection = true;
  sphere->material->k_reflection = 0.4f;

  sphere->material->bool_refraction = true;
  sphere->material->k_refraction = 0.6f;
  sphere->material->index_refraction = 1.33f; // water
*/
/*
  Vertex v2;
  v2.x = 2.0f;
  v2.y = 0.5f;
  v2.z = 3.0f;
  
  Sphere *sphere2 = new Sphere(v2,0.25f); // bubble above spout

  Phong bp3;

  bp3.ambient.r = 1.0f;
	bp3.ambient.g = 1.0f;
	bp3.ambient.b = 1.0f;
	bp3.diffuse.r = 1.0f;
	bp3.diffuse.g = 1.0f;
	bp3.diffuse.b = 1.0f;
	bp3.specular.r = 0.5f;
	bp3.specular.g = 0.5f;
	bp3.specular.b = 0.5f;
	bp3.power = 40.0f;

 	sphere2->material = &bp3;

  sphere2->material->bool_reflection = true;
  sphere2->material->k_reflection = 0.5f;

  sphere2->material->bool_refraction = true;
  sphere2->material->k_refraction = 0.5f;
  sphere2->material->index_refraction = 1.52f; // glass 

  Vertex v3;
  v3.x = 1.80f;
  v3.y = 1.5f;
  v3.z = 3.0f;
  
  Sphere *sphere3 = new Sphere(v3,0.25f); // bubble above spout

  sphere3->material = &bp3;

  sphere3->material->bool_reflection = true;
  sphere3->material->k_reflection = 0.5f;

  sphere3->material->bool_refraction = true;
  sphere3->material->k_refraction = 0.5f;
  sphere3->material->index_refraction = 1.52f; // glass 


*/
  // plane
  
/*
  Vertex a (1, -3, 1);
  Vertex b (1,2,1);
  Vector plane_normal = a.getDirection(b);
  Vector norm (0,1,0);
  Plane *floor = new Plane(a, norm);
*/
  
/*
  floor->material = &bp4;

  floor->material->bool_reflection = false;
  floor->material->k_reflection = 0.0f;
  
  floor->material->bool_refraction = false;
  floor->material->k_refraction = 0.0f;
  floor->material->index_refraction = 0.0f;
  */

  // link objects
 // pm->next = sphere;
  //sphere->next = sphere2;
  //sphere2->next = sphere3;
  

  // link objects
  //pm->next = sphere2;
  //sphere->next = sphere2;
  //sphere2->next = floor_pm;
  
  // link objects
  pm->next = background_pm;
  background_pm->next = floor_pm;
  floor_pm->next = left_wall;
  left_wall->next = right_wall;
  right_wall->next = ceiling_pm;
  ceiling_pm->next = sphere;
  
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
      int r_d = 3;
      // transparency (refraction) recursion depth
      int t_d = 10;

      raytrace(ray, pm, dl, colour, depth, r_d, t_d);

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
