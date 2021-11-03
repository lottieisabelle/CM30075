
class Lighting{
public:
    // coefficients
    float ambient;
    float specular;
    float diffuse;

    float red;
    float green;
    float blue;

    float* colour_hit(Hit &hit);
};