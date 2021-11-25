

class Scene{
public:
    float ambient_intensity;

    Scene()
    {
        ambient_intensity = 0.5;
    }

    Scene(float Ia)
    {
        ambient_intensity = Ia;
    }
};