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

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel; 
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	Colour sum_shade(0.0, 0.0, 0.0);
	int light_num = 1;
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.
		
	   // Implement shadows here if needed.

		Point3D origin = ray.intersection.point;
		Point3D light_position = curLight->light->get_position();
		Vector3D dir(light_position[0] - origin[0], light_position[1] - origin[1], light_position[2] - origin[2]);
		dir.normalize();
		Ray3D shadowRay = Ray3D(origin + 0.01 * dir, dir); 
		traverseScene(_root, shadowRay); 
		if (shadowRay.intersection.none) {
			curLight->light->shade(ray);
		}
		
		// Entended light on each light source to generate soft shadows
		sum_shade = sum_shade + 0.5*ray.col;
		for (int k=0; k < light_num; k++) {
			double varx = ((((double)rand() / (double)RAND_MAX)*2) - 1) * 0.2;
			double vary = ((((double)rand() / (double)RAND_MAX)*2) - 1) * 0.2;
			double varz = ((((double)rand() / (double)RAND_MAX)*2) - 1) * 0.2;
			Point3D tempShadowOrigin = origin;
			tempShadowOrigin[0] = tempShadowOrigin[0] + varx;
			tempShadowOrigin[1] = tempShadowOrigin[1] + vary;
			tempShadowOrigin[2] = tempShadowOrigin[2] + varz;			
			tempShadowOrigin = tempShadowOrigin + 0.2*dir; 	
			shadowRay = Ray3D(tempShadowOrigin, dir);
			traverseScene(_root, shadowRay);
			if (shadowRay.intersection.none) {
				curLight->light->shade(ray);				
				sum_shade = sum_shade + 0.1*ray.col;
			}
			
		}
		
		curLight->light->shade(ray);
		curLight = curLight->next;
		
	}
	sum_shade = (1.0/( (0.5) + (0.1)*(light_num-1)))*sum_shade;
	sum_shade.clamp();
	ray.col = sum_shade;	
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray, int recursiveDepth ) {
	Colour col(0.0, 0.0, 0.0); 
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray);
		col = ray.col; 
 
		// You'll want to call shadeRay recursively (with a different ray, 
		// of course) here to implement reflection/refraction effects.
		if (recursiveDepth > 0){
			// implement reflection
			if (ray.intersection.mat->refraction == 0){
				ray.dir.normalize();
				Vector3D mirrorDir = ((2 * (ray.intersection.normal.dot(-ray.dir))) * ray.intersection.normal) - (-ray.dir);
				Ray3D reflectRay(ray.intersection.point + 0.01 * mirrorDir, mirrorDir);
				shadeRay(reflectRay, recursiveDepth - 1);
				if (!reflectRay.intersection.none)
					col = col + ray.intersection.mat->specular * reflectRay.col;
			}
			
			// implement refraction - transparent object has refraction > 0
			if (ray.intersection.mat->refraction > 0) {
				Vector3D surfaceNormal = ray.intersection.normal;
				double cosAngle = surfaceNormal.dot(-ray.dir);
				double index = 1.0;
				// ray into material
				if (cosAngle > 0){
					index = 1.0 / ray.intersection.mat->refraction;
					// need to change surface normal
					surfaceNormal = -surfaceNormal;
				}				
				else{
					index = ray.intersection.mat->refraction;
				}
				double sinTheta = index * index * (1.0 - cosAngle * cosAngle);
				// not total internal refraction
				if (sinTheta <= 1.0){
					Vector3D refractionDir = index * -ray.dir - (index * cosAngle + sqrt(1.0 - sinTheta)) * surfaceNormal;
					refractionDir.normalize();
					Ray3D refractionRay(ray.intersection.point + 0.01 * refractionDir, refractionDir);			
					shadeRay(refractionRay, recursiveDepth - 1);
					if (!refractionRay.intersection.none)
						col = refractionRay.col;
				}			
			}
		}
	}
	col.clamp();

	return col; 
}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);
	bool antiAlias = false;
	int recursiveDepth = 5;

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);

	// Construct a ray for each pixel.
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			// Sets up ray origin and direction in view space, 
			// image plane is at z = -1.
			Point3D origin(0, 0, 0);
			Point3D imagePlane;
			
			//perform anti-Alias
			if (antiAlias){
				for (float frag1 = i;frag1 < i + 1.0f; frag1 += 0.5f){
					for (float frag2= j;frag2 < j + 1.0f; frag2 += 0.5f){
						imagePlane[0] = (-double(width)/2 + 0.5 + frag2)/factor;
                        imagePlane[1] = (-double(height)/2 + 0.5 + frag1)/factor;
                        imagePlane[2] = -1;
                        
                        Point3D rayOriginWorld = viewToWorld * imagePlane;
                        Vector3D rayDirWorld = rayOriginWorld - eye;
                        rayDirWorld.normalize();
                        
                        Ray3D ray(rayOriginWorld, rayDirWorld);
                        Colour col = shadeRay(ray, recursiveDepth);
                        _rbuffer[i*width+j] += int(col[0]*255*0.25f);
                        _gbuffer[i*width+j] += int(col[1]*255*0.25f);
                        _bbuffer[i*width+j] += int(col[2]*255*0.25f);
                    }
				}			
			}else{			
				imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
				imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
				imagePlane[2] = -1;

				// TODO: Convert ray to world space and call 
				// shadeRay(ray) to generate pixel colour. 				

				Point3D worldOrigin = viewToWorld * origin;
				Point3D worldImagePlane = viewToWorld * imagePlane;

				Vector3D direction(worldImagePlane[0] - worldOrigin[0], worldImagePlane[1] - worldOrigin[1], worldImagePlane[2] - worldOrigin[2]);

				Ray3D ray(worldImagePlane, direction);

				Colour col = shadeRay(ray, recursiveDepth); 

				_rbuffer[i*width+j] = int(col[0]*255);
				_gbuffer[i*width+j] = int(col[1]*255);
				_bbuffer[i*width+j] = int(col[2]*255);
			}
		}
	}
	flushPixelBuffer(fileName);
}


