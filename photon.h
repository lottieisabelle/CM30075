#include "vector.h"
#include "colour.h"
#include "hit.h"


class Photon{
public :
    Vector direction;
    Colour intensity;
    Hit intersection;
    string p_type;

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
        direction.normalise();
    }

    void set_hit(Hit &hit)
    {
        intersection = hit;
    }

    void set_type(string type)
    {
        p_type = type;
    }

};