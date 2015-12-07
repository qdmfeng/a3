/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"


bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, Material *material ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	Vector3D rayDir = worldToModel * ray.dir;
	Point3D rayOrigin = worldToModel * ray.origin;
	Vector3D normal(0.0, 0.0, 1.0);
    Point3D topRight(0.5, 0.5, 0.0);

	double dn = rayDir.dot(normal);

	if (dn > -0.00001 && 0.00001 > dn){ 
		// if the ray lies on the plane
		return false;
	}else{
		double t = (normal.dot(topRight - rayOrigin)) / dn;

		if (t < 0.00001)
			return false;

        //find intersection point
		Point3D intersect = rayOrigin + t * rayDir;
		double x = intersect[0];
		double y = intersect[1];
		double z = intersect[2];

		bool in = (z < 0.00001 && z > -0.00001 && x <= 0.5 &&
		            x >= -0.5 && y <= 0.5 && y >= -0.5);
		//update ray values
		if (in && (ray.intersection.none ||
		   (t < ray.intersection.t_value && !ray.intersection.none))){
			ray.intersection.none = false;
			ray.intersection.point = modelToWorld * intersect;
			ray.intersection.t_value = t;
			ray.intersection.normal = transNorm(worldToModel, normal);
			ray.intersection.normal.normalize();

			// when material has texture associated with it
			if (material->texture != NULL) {
				int width = floor(material->texture->width * (x + 0.5));
		        int height = floor(material->texture->height * (y + 0.5));
		        // index of the array
		        int index = height * material->texture->width + width;
		        // get texture value
		        int r_value = (int) *(material->texture->rarray + index);
		        double R = (double) r_value / 255;
		        int g_value = (int) *(material->texture->barray + index);
		        double G = (double) g_value / 255;
		        int b_value = (int) *(material->texture->garray + index);
		        double B = (double) b_value / 255;
				Colour textureCol = Colour(R, G, B);
		        //material->ambient = textureCol;
		        material->diffuse = textureCol;
    		}

			return true;
		}else{
			return false;
		}
	}
}


bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, Material * material ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	Point3D rayOrigin = worldToModel * ray.origin;
	Vector3D rayDir = worldToModel * ray.dir;
	Point3D center(0,0,0);
	
	double t, t1, t2;

    //quadratic formula coefficient
    double A = rayDir.dot(rayDir);
    double B = rayDir.dot(rayOrigin - center);
    double C = (rayOrigin - center).dot(rayOrigin - center) - 1;

    //Find delta
    double delta = 4*B*B-4*A*C;

    //If delta is negative there is no solution i.e., no intersection
    if ( delta < 0 )
        return false;
    
    // find solutions
    t1 = (-2*B + sqrt(delta))/(2.0 * A);
    t2 = (-2*B - sqrt(delta))/(2.0 * A);
    
    if ( t1 < 0 && t2 < 0 )
    //Both intersections are behind the view-plane, and are not visible
        return false;
    else if ( t1 < 0 )
    //p(t2) is a visible intersection but p(t1) is not
        t = t2;
    else if ( t2 < 0)
    //p(t1) is a visible intersection but p(t2) is not
        t = t1;
    else //Both intersections are in front of the view-plane.
		//Pick the closest intersection.
        t = fmin(t2,t1);

    //find intersection point
    Point3D intersect = rayOrigin + t * rayDir;
	double x = intersect[0];
	double y = intersect[1];
	double z = intersect[2];

    //update ray values
    if (ray.intersection.none || t < ray.intersection.t_value) {
	    Vector3D norm(x,y,z);
		norm.normalize();
        ray.intersection.t_value = t;
        ray.intersection.point = modelToWorld*intersect;
        ray.intersection.normal = worldToModel.transpose()*norm;
        ray.intersection.normal.normalize();
        ray.intersection.none = false;
        return true;
    }else
		return false;
}

bool UnitCylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, Material * material ) {
        // The intersection for a unit cylinder with the top and 
        // the base are unit circles with radius=1 centered at (0,0,z),
        // where z is -0.5, 0.5, respectively. 
        // So the height of the cylinder is 1.

	Point3D rayOrigin = worldToModel*ray.origin;
    Vector3D rayDir = worldToModel*ray.dir;
	Point3D sphereOrigin(0,0,0);
	double t1;
	double t2;
	double lambdaStar_1;
	double lambdaStar_2;

    //quadratic formula coefficient
	double A = rayDir[0]*rayDir[0] + rayDir[1]*rayDir[1];
	double B = (rayOrigin[0]*rayDir[0] + rayOrigin[1] * rayDir[1]);
	double C = rayOrigin[0]*rayOrigin[0] + rayOrigin[1]*rayOrigin[1] - 1;
	
    //Find delta
	double delta = B*B-A*C;
	
	Point3D intersectionPoint;
	Vector3D normal_1;
	Vector3D normal_2;
	Vector3D normal;

    //If delta is negative there is no solution i.e., no intersection
	if (delta<0)
		return false;

    // find solutions
	t1 = -B/A + sqrt(delta) / A;
	t2 = -B/A - sqrt(delta) / A;
	if (t1 < 0 && t2 < 0) 
	//Both intersections are behind the view-plane, and are not visible
		return false;
	else if (t1 > 0 && t2 < 0)
	//p(t1) is a visible intersection but p(t2) is not
		lambdaStar_2 = t1;
	else
		lambdaStar_2 = t2;
		
	//pick the closer solution
	t1 = (-0.5-rayOrigin[2])/rayDir[2];
	t2 = (0.5-rayOrigin[2])/rayDir[2];
	if (t1 < t2){
		lambdaStar_1 = t1;
		Point3D normal_temp(0,0,-1);
		normal_1 = normal_temp - sphereOrigin;
		normal_1.normalize();
	}
	else{
		lambdaStar_1 = t2;
		Point3D normal_temp(0,0,1);
		normal_1 = normal_temp - sphereOrigin;
		normal_1.normalize();
	}

	
	intersectionPoint = rayOrigin + lambdaStar_1 * rayDir;
	if (lambdaStar_1* lambdaStar_1 < 0.0001){
		return false;
	}
	
	//check if it intersects with the top or bottom by the first lambda
	if (intersectionPoint[0]*intersectionPoint[0] + intersectionPoint[1] * intersectionPoint[1] <= 1)
	{
		if (!ray.intersection.none > ray.intersection.t_value){
			return false;
		}
		ray.intersection.point = intersectionPoint;
		ray.intersection.normal = normal_1;
		ray.intersection.t_value = lambdaStar_1;
		ray.intersection.none = false;
		return true;
	}


	 //if not intersected with the base, check if intersects with 
	 //the side by second lambda
	intersectionPoint = rayOrigin + lambdaStar_2 * rayDir;
	if (lambdaStar_2 * lambdaStar_2 < 0.001)
		return false;
	
	normal_2[0] = intersectionPoint[0];
	normal_2[1] = intersectionPoint[1];
	normal_2[2] = 0;
	normal_2.normalize();



	if (intersectionPoint[2] < 0.5 && intersectionPoint[2] > -0.5)
	{
		if (!ray.intersection.none > ray.intersection.t_value)
			return false;
		
		ray.intersection.point = modelToWorld * intersectionPoint;
		Point3D normalll;
		normalll[0] = intersectionPoint[0];
		normalll[1] = intersectionPoint[1];
		normalll[2] = 0;
		ray.intersection.normal = modelToWorld * (normalll - sphereOrigin);
		ray.intersection.t_value = lambdaStar_2;
		ray.intersection.none = false;
		return true;
	}
	else
		return false;

}
