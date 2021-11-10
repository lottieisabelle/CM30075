

#include "ray.h"
#include "hit.h"

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

    Lighting(float Ia, float Ii)
    {
        ambient_intensity = Ia;
        diffuse_intensity = Ii;
    }

    void get_direction();

    void get_intensity();

    void get_position();


};