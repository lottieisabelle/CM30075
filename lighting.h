
#include "vertex.h"

class Lighting{
public:
    float red;
    float green;
    float blue;

    float ambient_intensity;
    float diffuse_intensity;

    Vertex position;

    Lighting()
    {
        ambient_intensity = 0.5;
        diffuse_intensity = 0.5;
    }

    Lighting(float Ia, float Ii, Vertex set_position)
    {
        ambient_intensity = Ia;
        diffuse_intensity = Ii;
        position = set_position;
    }

};