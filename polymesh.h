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

class PolyMesh {
public:
	int vertex_count;
	int triangle_count;
    Vertex *vertex;
	TriangleIndex *triangle;

	// coefficients
	float *ambient;
    float *diffuse;
	float *specular;

	void do_construct(char *file, Transform *transform);

	Vector getDirection(Vertex a, Vertex b);

	void intersection(Ray ray, Hit &hit);

	float* colour_hit(Hit &hit);

	float* colour_no_hit(Ray ray);

	float* calculate_lighting(Hit &hit, Lighting light, int flag);

	void set_coeffs(float ar, float ag, float ab, float dr, float dg, float db, float sr, float sg, float sb);
	
	PolyMesh(char *file);
	PolyMesh(char *file, Transform *transform);
};
