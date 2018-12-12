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

class Triangle {
    
public:
    Vec3f v0;
    Vec3f v1;
    Vec3f v2;
    Vec3f surfaceColor;
    
    Triangle(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, const Vec3f &surfaceColor) :
    v0(v0), v1(v1), v2(v2), surfaceColor(surfaceColor) {}
    
    bool rayTriangleIntersect(
                              const Vec3f &orig, const Vec3f &dir,
                              const Vec3f &v0, const Vec3f &v1, const Vec3f &v2,
                              float &t, Vec3f &N)
    {
        // compute plane's normal
        Vec3f v0v1 = v1 - v0;
        Vec3f v0v2 = v2 - v0;
        // no need to normalize
        N = v0v1 * v0v2; // N
        float area2 = N.length();
        
        // Step 1: finding P
        
        // check if ray and plane are parallel ?
        float NdotRayDirection = N.dot(dir);
        if (fabs(NdotRayDirection) < 0) // almost 0
            return false; // they are parallel so they don't intersect !
        
        // compute d parameter using equation 2
        float d = N.dot(v0);
        
        // compute t (equation 3)
        t = (N.dot(orig) + d) / NdotRayDirection;
        // check if the triangle is in behind the ray
        if (t < 0) return false; // the triangle is behind
        
        // compute the intersection point using equation 1
        Vec3f dir2 = dir;
        dir2.x = t * dir.x;
        dir2.y = t * dir.y;
        dir2.z = t * dir.z;
        Vec3f P = orig + dir2;
        
        // Step 2: inside-outside test
        Vec3f C; // vector perpendicular to triangle's plane
        
        // edge 0
        Vec3f edge0 = v1 - v0;
        Vec3f vp0 = P - v0;
        C = edge0 * vp0;
        if (N.dot(C) < 0) return false; // P is on the right side
        
        // edge 1
        Vec3f edge1 = v2 - v1;
        Vec3f vp1 = P - v1;
        C = edge1 * vp1;
        if (N.dot(C) < 0)  return false; // P is on the right side
        
        // edge 2
        Vec3f edge2 = v0 - v2;
        Vec3f vp2 = P - v2;
        C = edge2 * vp2;
        if (N.dot(C) < 0) return false; // P is on the right side;
        
        return true; // this ray hits the triangle
    }
};
