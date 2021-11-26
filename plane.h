#include "vertex.h"
#include "vector.h"
#include "object.h"

class Plane : public Object {
	Vertex point;
    Vector normal;
public:
	Plane(Vertex p, Vector n);
	void intersection(Ray ray, Hit &hit);
};