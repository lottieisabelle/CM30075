
class Colour{
    public:
        float red;
        float green;
        float blue;

        Colour(){

        }

        Colour(float r, float g, float b){
            red = r;
            green = g;
            blue = b;
        }

        void multiply(Colour scale){
            red = red * scale.red;
            green = green * scale.green;
            blue = blue * scale.blue;
        }
};