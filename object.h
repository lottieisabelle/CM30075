/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Object is the base class for objects.
#ifndef _OBJECT_H_
#define _OBJECT_H_

#include "ray.h"
#include "hit.h"

class Object {
public:

	Object *next;
	// coefficients
	float *ambient;
    float *diffuse;
	float *specular;

	Object()
	{
		next = (Object *)0;
	}
	
	virtual void intersection(Ray ray, Hit &hit)
	{

	}

	void set_coeffs(float ar, float ag, float ab, float dr, float dg, float db, float sr, float sg, float sb)
	{

		ambient[0] = ar;
		ambient[1] = ag;
		ambient[2] = ab; 

		diffuse[0] = dr;
		diffuse[1] = dg;
		diffuse[2] = db;

		specular[0] = sr;
		specular[1] = sg;
		specular[2] = sb; 
	}
};

#endif
