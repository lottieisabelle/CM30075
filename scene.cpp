
#include "scene.h"


void Scene::addObject(Object *object)
{
    object_list.push_back(object);
}

void Scene::addLight(Lighting light)
{
    light_list.push_back(light);
}

void Scene::render_image()
{

}

Colour Scene::raytracer(Ray ray, Hit &hit)
{
    Colour final_colour = ambient_intensity;

    for(Object *object : object_list){
        object->intersection(ray, hit);
    }

    if (hit.flag != true){
        // if hit nothing, make background black
        return Colour (0.0, 0.0, 0.0);
    }

    for(Lighting light : light_list){
        Hit shadow_hit;
        shadow_hit.flag = false;
        shadow_hit.t = 99999999;
        
        Vertex shadow_point = Vertex (hit.position.x-0.00222, hit.position.y-0.00222, hit.position.z-0.00222);
        Vector shadow_dir = shadow_point.getDirection(light.position);
        Ray shadow_ray (shadow_point, shadow_dir);

        if (shadow_hit.flag == true){
            continue;
        }

        // do diffuse and specular here
    }

    
    // TODO : reverse y direction?
    //int h = screen_height-1-(ray_y+1)*(screen_height/2);

    // determine which lighting calculation is needed
    if (shooting_hit.flag==true){
    // calculate if shadows here
    

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


    return final_colour;
}