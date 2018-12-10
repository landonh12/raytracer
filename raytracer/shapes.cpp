//
//  sphere.cpp
//  raytracer
//
//  Created by Landon Haugh on 12/6/18.
//  Copyright © 2018 Landon Haugh. All rights reserved.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include "vec3.cpp"

typedef Vec3<float> Vec3f;

class Sphere {
    
public:
    Vec3f center;                           /// position of the sphere
    float radius, radius2;                  /// sphere radius and radius^2
    Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
    float transparency, reflection;         /// surface transparency and reflectivity
    Sphere(
           const Vec3f &c,
           const float &r,
           const Vec3f &sc,
           const float &refl = 0,
           const float &transp = 0,
           const Vec3f &ec = 0) :
    center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
    transparency(transp), reflection(refl)
    { /* empty */ }
    
    // Returns t0 and t1 through pointers
    // Returns true if intersection
    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
    {
        Vec3f l = center - rayorig;
        float tca = l.dot(raydir);
        if (tca < 0) return false;
        float d2 = l.dot(l) - tca * tca;
        if (d2 > radius2) return false;
        float thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;
        
        return true;
    }
};

class Ray {
    
public:
    Ray(const Vec3f &orig, const Vec3f &dir) : orig(orig), dir(dir)
    {
        invdir.x = 1 / dir.x;
        invdir.y = 1 / dir.y;
        invdir.z = 1 / dir.z;
        sign[0] = (invdir.x < 0);
        sign[1] = (invdir.y < 0);
        sign[2] = (invdir.z < 0);
    }
    Vec3f orig, dir; // ray orig and dir
    Vec3f invdir;
    int sign[3];
};

class Cube {
    
public:
    
    Vec3f bounds[2];
    Vec3f c;
    Cube(const Vec3f &b0, const Vec3f &b1, const Vec3f &color) {
        bounds[0] = b0;
        bounds[1] = b1;
        c = color;
    }
    
    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t) const {
        
        Ray r = Ray(rayorig, raydir);
        
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        
        tmin = (bounds[r.sign[0]].x - r.orig.x) * r.invdir.x;
        tmax = (bounds[1-r.sign[0]].x - r.orig.x) * r.invdir.x;
        tymin = (bounds[r.sign[1]].y - r.orig.y) * r.invdir.y;
        tymax = (bounds[1-r.sign[1]].y - r.orig.y) * r.invdir.y;
        
        if ((tmin > tymax) || (tymin > tmax))
            return false;
        
        if (tymin > tmin)
            tmin = tymin;
        if (tymax < tmax)
            tmax = tymax;
        
        tzmin = (bounds[r.sign[2]].z - r.orig.z) * r.invdir.z;
        tzmax = (bounds[1-r.sign[2]].z - r.orig.z) * r.invdir.z;
        
        if ((tmin > tzmax) || (tzmin > tmax))
            return false;
        
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;
        
        t = tmin;
        
        if (t < 0) {
            t = tmax;
            if (t < 0) return false;
        }
        
        return true;
    }
};
