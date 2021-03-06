/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Material is the base class for materials.

#pragma once

#include "vector.h"
#include "colour.h"

class Material {
public:
	float k_reflection;
	float k_refraction;
	float index_refraction;
	float ior_object;
	float ior_surround;
	bool bool_reflection;
	bool bool_refraction;
	bool bool_specular;

	virtual void compute_base_colour(Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}
	virtual void compute_light_colour(Vector &viewer, Vector &normal, Vector &ldir, Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}

	virtual void compute_diffuse(Vector &viewer, Vector &normal, Vector &ldir, Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}

	virtual void compute_specular(Vector &viewer, Vector &normal, Vector &ldir, Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}

	virtual float prob_diff(){
		return 0.0f;
	}

	virtual float prob_ref(){
		return 0.0f;
	}

	virtual float prob_spec(){
		return 0.0f;
	}

	virtual void get_diffuse(Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}

	virtual void get_specular(Colour &result)
	{
		result.r = 0.0f;
		result.g = 0.0f;
		result.b = 0.0f;
	}

};
