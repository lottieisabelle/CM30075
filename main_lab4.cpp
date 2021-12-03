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
#include "photon.h"
#include "point_light.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <random>

using namespace std;

void clamp(float &x){
  if (x <= -1){
    x = -1.0;
  } else if (x >= 1){
    x = 1.0;
  } 
}

float fresnel(Ray ray, Hit &best_hit){
  // kr to be returned from this function
  float kr;

  // calculate cos of angle between incident ray and normal to surface
  float cos_i = best_hit.normal.dot(ray.direction);
  clamp(cos_i); // limit between -1 and 1

  float n_out;
  float n_in;

  // assume
  n_out = best_hit.what->material->ior_surround;
  n_in = best_hit.what->material->ior_object;

  // unless
  if (cos_i < 0){
    n_out = best_hit.what->material->ior_object;
    n_in = best_hit.what->material->ior_surround;
  }

  float n = n_in/n_out;
  best_hit.what->material->index_refraction = n;

  float sin_t = n_out / n_in * sqrtf(std::max(0.f, 1-cos_i *cos_i));
  if (sin_t >= 1){
    // total internal reflection
    kr = 1;
    return kr;
  }

  float cos_t = sqrtf(std::max(0.0f, 1 - sin_t *sin_t));
  cos_i = fabsf(cos_i);

  float fr_1 = ((n_in*cos_i)-(n_out*cos_t)) / ((n_in*cos_i)+(n_out*cos_t));
  float fr_2 = ((n_out*cos_i)-(n_in*cos_t)) / ((n_out*cos_i)+(n_in*cos_t));
  kr = (fr_1*fr_1 + fr_2*fr_2) /2;

  return kr;
}

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
  if (d <= 0){
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

    if(best_hit.what->material->bool_refraction){
      // do fresnel equation
      float kr = fresnel(ray, best_hit);
      float kt = 1.0 - kr;

      best_hit.what->material->k_reflection = kr;
      best_hit.what->material->k_refraction = kt;

      //printf("%f , %f \n", kr, kt);
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
      raytrace(r_ray, objects, lights, col, depth, d-1);

      col.scale(kr_col);
      colour.add(col);
      
    }

    // compute refraction ray if material supports it.
    if(best_hit.what->material->bool_refraction && best_hit.what->material->k_refraction > 0.0)
    {
      float n = best_hit.what->material->index_refraction;
      // tray.dir = refraction(ray.dir, hit.normal, hit.ior);
      
      // cos θi = N.I
      // I = incident ray direction vector
      // N = normal to surface direction vector
      float cos_i = best_hit.normal.dot(ray.direction);
      clamp(cos_i);

      // cos θt = sqrt(1 – (1/η2) * (1 - cos2 θi) )
      float n2 = n*n;
      float cos_i2 = cos_i*cos_i;
      float cos_t2 = 1 - (1/n2) * (1 - cos_i2);

      float cos_t = sqrt(cos_t2);
      clamp(cos_t);

      // T = 1/η * I – (cos θt – (1/η)* cos θi ) * N
      Vector T_dir = 1/n * ray.direction - (cos_t - (1/n) * cos_i) * best_hit.normal;
      
      // tray.pos = hit.pos + small_e * tray.dir;
      Vertex T_pos;
      T_pos.x = best_hit.position.x + 0.00222 * T_dir.x;
      T_pos.y = best_hit.position.y + 0.00222 * T_dir.y;
      T_pos.z = best_hit.position.z + 0.00222 * T_dir.z;
      
      Ray T_ray (T_pos, T_dir);

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
  //pm->material->k_reflection = 0.4f;
  pm->material->bool_refraction = false;
  //pm->material->ior_object = 1.52f; // glass
  //pm->material->ior_surround = 1.0003f; // air


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
	bp6.ambient.g = 1.0f;
	bp6.ambient.b = 1.0f;
	bp6.diffuse.r = 1.0f;
	bp6.diffuse.g = 1.0f;
	bp6.diffuse.b = 1.0f;
	bp6.specular.r = 0.2f;
	bp6.specular.g = 0.2f;
	bp6.specular.b = 0.2f;
	bp6.power = 40.0f;

  // back wall
  background_pm->material = &bp6;
  background_pm->material->bool_reflection = false;  
  background_pm->material->bool_refraction = false; 

  // floor
  floor_pm->material = &bp4;
  floor_pm->material->bool_reflection = false;
  floor_pm->material->bool_refraction = false;

  // left wall
  left_wall->material = &bp5;
  left_wall->material->bool_reflection = false;
  left_wall->material->bool_refraction = false;

  // right wall
  right_wall->material = &bp5;
  right_wall->material->bool_reflection = false;  
  right_wall->material->bool_refraction = false;

  // ceiling
  ceiling_pm->material = &bp4;
  ceiling_pm->material->bool_reflection = false;  
  ceiling_pm->material->bool_refraction = false;

  // create bubbles
  // lower bubble
  Vertex v;
  v.x = 3.0f; 
  v.y = -1.0f; 
  v.z = 8.0f; 
  
  Sphere *sphere = new Sphere(v, 0.4);
  Phong bp2;
  bp2.ambient.r = 0.6f;
	bp2.ambient.g = 0.6f;
	bp2.ambient.b = 0.6f;
	bp2.diffuse.r = 0.6f;
	bp2.diffuse.g = 0.6f;
	bp2.diffuse.b = 0.6f;
	bp2.specular.r = 0.4f;
	bp2.specular.g = 0.4f;
	bp2.specular.b = 0.4f;
	bp2.power = 40.0f;

  sphere->material = &bp2;
  sphere->material->bool_reflection = true;
  sphere->material->bool_refraction = true;
  sphere->material->ior_object = 1.38f; // soap bubbles
  sphere->material->ior_surround = 1.0003f; // air
  //sphere->material->ior_object = 1.52f; // glass
  //sphere->material->ior_object = 1.33f; // water

  // top bubble
  Vertex v2;
  v2.x = 2.0f;
  v2.y = 1.0f;
  v2.z = 8.0f;
  
  Sphere *sphere2 = new Sphere(v2,0.6f);
  Phong bp3;
  bp3.ambient.r = 0.6f;
	bp3.ambient.g = 0.6f;
	bp3.ambient.b = 0.6f;
	bp3.diffuse.r = 0.6f;
	bp3.diffuse.g = 0.6f;
	bp3.diffuse.b = 0.6f;
	bp3.specular.r = 0.4f;
	bp3.specular.g = 0.4f;
	bp3.specular.b = 0.4f;
	bp3.power = 40.0f;

 	sphere2->material = &bp3;

  sphere2->material->bool_reflection = true;
  sphere2->material->bool_refraction = true;
  //sphere2->material->ior_object = 1.33f; // water
  sphere2->material->ior_object = 1.38f; // soap bubbles
  sphere2->material->ior_surround = 1.0003f; // air

  // big bubble
  Vertex v3;
  v3.x = 0.0f;
  v3.y = 2.0f;
  v3.z = 4.0f;
  
  Sphere *sphere3 = new Sphere(v3,0.6f);
  Phong bp7;
  bp7.ambient.r = 0.6f;
	bp7.ambient.g = 0.6f;
	bp7.ambient.b = 0.6f;
	bp7.diffuse.r = 0.6f;
	bp7.diffuse.g = 0.6f;
	bp7.diffuse.b = 0.6f;
	bp7.specular.r = 0.4f;
	bp7.specular.g = 0.4f;
	bp7.specular.b = 0.4f;
	bp7.power = 40.0f;

 	sphere3->material = &bp7;

  sphere3->material->bool_reflection = true;
  sphere3->material->bool_refraction = true;
  //sphere3->material->ior_object = 1.33f; // water
  sphere3->material->ior_object = 1.38f; // soap bubbles
  sphere3->material->ior_surround = 1.0003f; // air

  // link objects
  pm->next = background_pm;
  background_pm->next = floor_pm;
  floor_pm->next = left_wall;
  left_wall->next = right_wall;
  right_wall->next = ceiling_pm;
  ceiling_pm->next = sphere;
  sphere->next = sphere2;
  //sphere2->next = sphere3;
  
  // generate shooting ray from camera point
  Ray ray;

  ray.position.x = 0.0001f;
  ray.position.y = 0.0f;
  ray.position.z = 0.0f;

  // create directional light
  DirectionalLight *dl = new DirectionalLight(Vector(0.5f, -0.2f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f));
  // an approximation that simulates the sun
  // x how much pointing left and right
  // y how much up or down
  // z into and out of the scene, more extemely into the scene = more positive

  PointLight *pl = new PointLight(Vertex(-5.0f,2.0f,-10.0f), Colour(1.0f,1.0f,1.0f,0.0f));

  // photon mapping here

  // random number generator set up code, for random direction vectors
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> dist(0.0, 1.0); //range defined here (inclusive) TODO define range
  
  int n = 10; // number of photons TODO set number

  // for n photons
  for (int i=0; i<n; ++i){
    // create photon
    Photon p;

    // send out into picture, meaning give photon direction vector
    // create random direction vectors all across image

    float p_x = dist(mt);
    float p_y = dist(mt);
    float p_z = dist(mt);

    Vector photon_dir (p_x, p_y, p_z);
    photon_dir.normalise();

    p.set_dir(photon_dir);
    //printf("%f , %f , %f\n",p_x, p_y, p_z);

    // determine hit point (use intersection?) basically like raytrace?


    // store hit position, type of photon (d/i/s), intensity of photon? in photon object in kdtree

  } 




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

      // reflection and refraction recursion depth
      int d = 10;

      // uses point light pl, change to dl to use directional light
      raytrace(ray, pm, pl, colour, depth, d);

      fb->plotPixel(x, y, colour.r, colour.g, colour.b);
      //fb->plotDepth(x,y, depth);
    }

    cerr << "*" << flush;
  }

  // blur middle line maybe
  // for image height in pixels
  // get 2 middle pixels colour
  // x/2 and x/2 +1
  // get average of the two 
  // set pixels again
  
  // Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");
  //  fb->writeDepthFile((char *)"depth.ppm");
  return 0;
  
}
