/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[])
{	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	int width = 320; 
	int height = 240; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, 0, 10);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2, 0, NULL );
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0, NULL );
	Material blue( Colour(0, 0, 0), Colour(0.3, 0.3, 0.3), 
			Colour(0.2, 0.2, 0.7), 
			5.0, 0, NULL );
	Material glass( Colour(0.1, 0.1, 0.1), Colour(0.1, 0.1, 0.1), 
			Colour(0.9, 0.9, 0.9), 
			100, 1.5, NULL );

	// Define texture for shading
	unsigned char *rarray, *garray, *barray;
	unsigned long int bmpwidth;
	long int bmpheight;
	bool err = bmp_read("LightningTexture.bmp", &bmpwidth, &bmpheight, &rarray, &garray, &barray);
	if (err) { return -1; }

	Texture lightning(bmpwidth, bmpheight, rarray, garray, barray);

	Material wall( Colour(0.1, 0.1, 0.1), Colour(0.1, 0.1, 0.1), 
		Colour(0.1, 0.1, 0.1), 
		5.0, 1.8, &lightning );


	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 5), 
				Colour(0.9, 0.9, 0.9) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* cylinder1 = raytracer.addObject(new UnitCylinder(), &gold );
	SceneDagNode* cylinder2 = raytracer.addObject(new UnitCylinder(), &gold );
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &glass );
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &blue );
	SceneDagNode* wall1 = raytracer.addObject( new UnitSquare(), &wall );
	SceneDagNode* wall2 = raytracer.addObject( new UnitSquare(), &wall );
	SceneDagNode* wall3 = raytracer.addObject( new UnitSquare(), &wall );
	SceneDagNode* wall4 = raytracer.addObject( new UnitSquare(), &wall );

	// Apply some transformations to the unit square.
	double factor1[3] = { 2.0, 2.0, 2.0 };
	double factor2[3] = { 14.0, 14.0, 14.0 };
	double factor3[3] = { 0.75, 0.75, 1.0 };

	raytracer.translate(sphere, Vector3D(0, 0, -5));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);
	
	raytracer.translate(plane, Vector3D(0, 0, -7));	
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	raytracer.translate(wall1, Vector3D(7, 0, 0));	
	raytracer.rotate(wall1, 'y', -90); 
	raytracer.scale(wall1, Point3D(0, 0, 0), factor2);

	raytracer.translate(wall2, Vector3D(-7, 0, 0));	
	raytracer.rotate(wall2, 'y', 90); 
	raytracer.scale(wall2, Point3D(0, 0, 0), factor2);

	raytracer.translate(wall3, Vector3D(0, 7, 0));
	raytracer.rotate(wall3, 'x', 90); 	
	raytracer.scale(wall3, Point3D(0, 0, 0), factor2);

	raytracer.translate(wall4, Vector3D(0, -7, 0));	
	raytracer.rotate(wall4, 'x', -90); 
	raytracer.scale(wall4, Point3D(0, 0, 0), factor2);

    raytracer.translate(cylinder1, Vector3D(0, 4, -6));
    raytracer.scale(cylinder1, Point3D(0, 0, 0), factor3);

    raytracer.translate(cylinder2, Vector3D(0, -4, -6));
    raytracer.scale(cylinder2, Point3D(0, 0, 0), factor3);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	
	// Render it from a different point of view.
	Point3D eye2(5, 4, 1);
	Vector3D view2(-4, -2, -6);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");

	Point3D eye3(6, 5, -4);
	Vector3D view3(-4, -2, 0);
	raytracer.render(width, height, eye3, view3, up, fov, "view3.bmp");
	
	return 0;
}

