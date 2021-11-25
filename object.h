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

#include "colour.h"

class Object {
public:

	Object *next;
	// coefficients
	Colour ambient;
	Colour diffuse;
	Colour specular;

	Object()
	{
		next = (Object *)0;
	}
	
	virtual void intersection(Ray ray, Hit &hit)
	{

	}

	void set_coeffs(float ar, float ag, float ab, float dr, float dg, float db, float sr, float sg, float sb)
	{
		ambient.red = ar;
		ambient.green = ag;
		ambient.blue = ab;

		diffuse.red = dr;
		diffuse.green = dg;
		diffuse.blue = db;

		specular.red = sr;
		specular.green = sg;
		specular.blue = sb; 
	}
};

#endif
