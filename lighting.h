
#include "vertex.h"

class Lighting{
public:
    float red;
    float green;
    float blue;

    float diffuse_intensity;

    Vertex position;

    Lighting()
    {
        diffuse_intensity = 0.5;
        position = Vertex (-4,-1,1);
    }

    Lighting(float Ii, Vertex set_position)
    {
        diffuse_intensity = Ii;
        position = set_position;
    }

};