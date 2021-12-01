#include "vector.h"
#include "colour.h"
#include "hit.h"


class Photon{
public :
    Vector direction;
    Colour intensity;
    Hit intersection;
    bool direct;
    bool indirect;
    bool shadow;

    Photon()
    {
        direct = false;
        indirect = false;
        shadow = false;
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

};