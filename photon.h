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

};