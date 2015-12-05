/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	Vector3D lightDirection(_pos[0] - ray.intersection.point[0], _pos[1] - ray.intersection.point[1], _pos[2] - ray.intersection.point[2]);
	lightDirection.normalize();

	// calculate the ambient component
	Colour ambientIntensity = _col_ambient * ray.intersection.mat->ambient;

	// calculate the diffuse component
	double diffusePara = fmax(0, ray.intersection.normal.dot(lightDirection));
	Colour diffuseIntensity = diffusePara * _col_diffuse * ray.intersection.mat->diffuse;

	// calculate the specular component
	Vector3D mirrorDir = ((2 * (ray.intersection.normal.dot(-lightDirection))) * ray.intersection.normal) - (-lightDirection);
	mirrorDir.normalize();
	ray.dir.normalize();
	double specularPara = pow(fmax(0, mirrorDir.dot(ray.dir)), ray.intersection.mat->specular_exp);
	Colour specularIntensity = specularPara * _col_specular * ray.intersection.mat->specular;



	// uncomment below line to show scene signature without using Phong model
	ray.col = ray.intersection.mat->ambient + ray.intersection.mat->diffuse;

	// uncomment below line to shade the scene using Phong model without specular
	ray.col = ambientIntensity + diffuseIntensity;

	// shade the scene using Phong model
	ray.col = ambientIntensity + diffuseIntensity + specularIntensity;
	
	// clamp color values
	ray.col.clamp();

}

void PointLight::shadeShadow( Ray3D& ray ) {
	Colour ambientIntensity = _col_ambient * ray.intersection.mat->ambient;
	ray.col = ambientIntensity;
	ray.col.clamp();

}

