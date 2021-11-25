#include <vector>
#include "object.h"
#include "lighting.h"
#include "colour.h"

class Scene{
public:
    Colour ambient_intensity;
    std::vector<Object*> object_list;
    std::vector<Lighting> light_list;

    Scene()
    {
        ambient_intensity = Colour (0.5, 0.5, 0.5);
    }

    Scene(float Ia)
    {
        ambient_intensity = Colour (Ia, Ia, Ia);
    }

    void addObject(Object *object);
    void addLight(Lighting light);

    void render_image();

    Colour raytracer(Ray ray, Hit &hit);

    
};