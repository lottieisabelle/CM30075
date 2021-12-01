#pragma once
#include "light.h"

class PointLight : public Light{
public:
    Vertex position;
    Colour intensity;

    PointLight();
    PointLight(Vertex pos, Colour col);

    bool get_direction(Vertex &surface, Vector &dir);
	void get_intensity(Vertex &surface, Colour &intensity);

};