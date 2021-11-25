
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

Colour Scene::calculate_lighting(Lighting light, Hit &hit)
{
    Colour col;

    Vector L = hit.position.getDirection(light.position);
    L.normalise();
    Vector I = light.position.getDirection(hit.position);
    I.normalise();
    Vector R;
    hit.normal.reflection(I,R);
    R.normalise();
    Vector V = hit.position.getDirection(Vertex (0,0,0));
    V.normalise();
    int n = 40;

    col.red = light.diffuse_intensity *  ( hit.what->diffuse.red*(hit.normal.dot(L))  +  hit.what->specular.red*  pow(R.dot(V),n));
    col.green = light.diffuse_intensity *  ( hit.what->diffuse.green*(hit.normal.dot(L))  +  hit.what->specular.green*  pow(R.dot(V),n));
    col.blue = light.diffuse_intensity *  ( hit.what->diffuse.blue*(hit.normal.dot(L))  +  hit.what->specular.blue*  pow(R.dot(V),n));

    return col;
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

    final_colour.multiply(hit.what->ambient);

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

        // adds diffuse and specular light
        final_colour = calculate_lighting(light, hit);
    }

    return final_colour;
}