#pragma once

#include "vector.h"
#include "colour.h"
#include "object.h"


class Photon {
public :
    Vector direction;
    Vertex position;
    Colour intensity;
    Object *what;
    char p_type;

    Photon()
    {

    }

    void set_intensity(Colour col)
    {
        intensity = col;
    }

    void set_dir(Vector dir)
    {
        direction = dir;
        //direction.normalise();
    }

    void set_type(char type)
    {
        p_type = type;
    }

};