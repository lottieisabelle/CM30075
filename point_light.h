#pragma once
#include "light.h"

class PointLight : public Light{
public:
    Vertex position;
    Colour intensity;
    Vector direction;

    PointLight();
    PointLight(Vertex pos, Vector dir, Colour col);

    bool get_direction(Vertex &surface, Vector &dir);
	void get_intensity(Vertex &surface, Colour &intensity);
    Vertex get_position();
    Vector get_direction();

};