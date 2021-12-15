
#include "point_light.h"

PointLight::PointLight()
{
    Light();
}

PointLight::PointLight(Vertex pos, Vector dir, Colour col)
{
    Light();

    position = pos;
    intensity = col;
    direction = dir;
}

bool PointLight::get_direction(Vertex &surface, Vector &dir)
{
	dir = position.getDirection(surface);
    dir.normalise();

	return true;
}

void PointLight::get_intensity(Vertex &surface, Colour &level)
{
	level = intensity;
}

Vertex PointLight::get_position(){
    return position;
}

Vector PointLight::get_direction(){
    return direction;
}
