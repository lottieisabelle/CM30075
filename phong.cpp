/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include "phong.h"

#include <math.h>
#include <algorithm> 

// A simple Phong based lighting model

void Phong::compute_base_colour(Colour &result)
{
	result.r = ambient.r;
	result.g = ambient.g;
	result.b = ambient.b;
}

void Phong::compute_light_colour(Vector &viewer, Vector &normal, Vector &ldir, Colour &result)
{

	float diff;

	Vector tolight;
	Vector toviewer;

	result.r=0.0f;
	result.g=0.0f;
	result.b=0.0f;

	tolight = ldir;
	tolight.negate();

	toviewer = viewer;
	toviewer.negate();

	diff = normal.dot(tolight);
	
	if (diff < 0.0f) // light is behind surface
	{
		return;
	}

	// diffuse

	result.r += diffuse.r * diff;
	result.g += diffuse.g * diff;
	result.b += diffuse.b * diff;

	// the specular component

	Vector r;
	
	normal.reflection(tolight, r);
	r.normalise();

	float h;

	h = r.dot(toviewer);

	if (h > 0.0f)
	{
		float p = (float)pow(h, power);

		result.r += specular.r * p;
		result.g += specular.g * p;
		result.b += specular.b * p;
	}
}

void Phong::compute_diffuse(Vector &viewer, Vector &normal, Vector &ldir, Colour &result)
{

	float diff;

	Vector tolight;
	Vector toviewer;

	result.r=0.0f;
	result.g=0.0f;
	result.b=0.0f;

	tolight = ldir;
	tolight.negate();

	toviewer = viewer;
	toviewer.negate();

	diff = normal.dot(tolight);
	
	if (diff < 0.0f) // light is behind surface
	{
		return;
	}

	// diffuse

	result.r += diffuse.r * diff;
	result.g += diffuse.g * diff;
	result.b += diffuse.b * diff;
}

void Phong::compute_specular(Vector &viewer, Vector &normal, Vector &ldir, Colour &result)
{
	// the specular component

	Vector r;

	Vector tolight;
	Vector toviewer;

	result.r=0.0f;
	result.g=0.0f;
	result.b=0.0f;

	tolight = ldir;
	tolight.negate();

	toviewer = viewer;
	toviewer.negate();
	
	normal.reflection(tolight, r);
	r.normalise();

	float h;

	h = r.dot(toviewer);

	if (h > 0.0f)
	{
		float p = (float)pow(h, power);

		result.r += specular.r * p;
		result.g += specular.g * p;
		result.b += specular.b * p;
	}
}

float Phong::prob_diff(){
	float prob_diff = std::max(diffuse.r, diffuse.g);
	prob_diff = std::max(prob_diff, diffuse.b);

	//float prob_ref = std::max((diffuse.r +specular.r), (diffuse.g + specular.g));
	//prob_ref = std::max(prob_ref, (diffuse.b + specular.b));
	//float prob_diff = ((diffuse.r + diffuse.g + diffuse.b) / (diffuse.r + diffuse.g + diffuse.b + specular.r + specular.g + specular.b)) * prob_ref;
	return prob_diff;	
}

float Phong::prob_spec(){
	float prob_spec = std::max(specular.r, specular.g);
	prob_spec = std::max(prob_spec, specular.b);
	return prob_spec;
}

float Phong::prob_ref(){
	float prob_ref = std::max((diffuse.r +specular.r), (diffuse.g + specular.g));
	prob_ref = std::max(prob_ref, (diffuse.b + specular.b));
	return prob_ref;
}

void Phong::get_diffuse(Colour &result){
	result.r = diffuse.r;
	result.g = diffuse.g;
	result.b = diffuse.b;
}

void Phong::get_specular(Colour &result){
	result.r = specular.r;
	result.g = specular.g;
	result.b = specular.b;
}
