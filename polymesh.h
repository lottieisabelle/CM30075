/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include "vertex.h"
#include "transform.h"
#include <math.h>

#include "ray.h"
#include "hit.h"

#include "lighting.h"

typedef int TriangleIndex[3];

class PolyMesh : public Object {
public:
	int vertex_count;
	int triangle_count;
    Vertex *vertex;
	TriangleIndex *triangle;

	void do_construct(char *file, Transform *transform, int x);

	void intersection(Ray ray, Hit &hit);

	float* colour_hit(Hit &hit);

	float* colour_no_hit(Ray ray);
	
	PolyMesh(char *file, int x);
	PolyMesh(char *file, Transform *transform, int x);
};
