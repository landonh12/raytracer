//
//  main.cpp
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
#include "shapes.cpp"

#if defined __linux__ || defined __APPLE__
// "Compiled for Linux
#else
// Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif

#define MAX_RAY_DEPTH 2

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

//TODO
Vec3f trace(const Vec3f &rayorig, const Vec3f &raydir, const std::vector<Sphere> &spheres, const int &depth) {
    
    float dist = INFINITY;       // Distance to intersection
    const Sphere *sphere = NULL; // Create sphere object
    float t0;                    // t0 sphere intersection
    float t1;                    // t1 sphere intersection
    
    // Check each object for intersection
    for(int i = 0; i < spheres.size(); i++) {
        // t0, t1 are intersection points
        t0 = INFINITY;
        t1 = INFINITY;
        // Gives t0 and t1 back from intersect function
        if(spheres[i].intersect(rayorig, raydir, t0, t1)) {
            // If we are already inside the sphere, set t0 to t1.
            if(t0 < 0) {
                t0 = t1;
            }
            // if t0 is less than infinity, set distance to t0 and assign sphere object
            if(t0 < dist) {
                dist = t0;
                sphere = &spheres[i];
            }
        }
    }
    
    //if(tc < t0) printf("t0: %f, t1: %f, tc: %f", t0,t1,tc);
    
    // If there's no cube intersections, then there's no hit!
    if(!sphere) {
        return Vec3f(0.0, 0.0, 0.0);
    }
    
    Vec3f sColor = Vec3f(0.0, 0.0, 0.0);  // Surface color
    Vec3f pHit = rayorig + raydir * dist; // Point of intersection with the sphere
    Vec3f nHit = pHit - sphere->center;   // Normal of the hit
    nHit.normalize(); // Normalize normal dir
    
    // START REFLECTION AND REFRACTION
    
    // Biasing
    float bias = 1e-2;
    
    // Shouldn't be inside the object at first
    bool inside = false;
    // If the ray direction * normal is > 0, then we are inside the object
    if(raydir.dot(nHit) > 0) {
        nHit = -nHit;
        inside = true;
    }
    // If there's a sphere object, and it has transparecny or reflection, start a recursive ray tracing function
    if(sphere && (sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
        // This is just some fresnel effect stuff, got from the example.
        float facingratio = -raydir.dot(nHit);
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
        
        // Create reflection vector
        Vec3f reflection = Vec3f(0.0,0.0,0.0);
        if(sphere->reflection) {
            // Reflection direction is ray direction - normal * 2 * ray direction dot normal
            Vec3f refldir = raydir - nHit * 2 * raydir.dot(nHit);
            // Normalize reflection direction
            refldir.normalize();
            // Trace reflection ray. Ray origin is pHit + nHit * bias
            reflection = trace(pHit + nHit * bias, refldir, spheres, depth + 1);
        }
        // Create refraction vector
        Vec3f refraction = 0;
        if(sphere->transparency) {
            // Index of reflection and entering/exiting (inverse of index or reflection)
            float ior = 1.345, eta = (inside) ? ior : 1 / ior;
            // cosi is inverse normal dot ray direction
            float cosi = -nHit.dot(raydir);
            float k = 1 - eta * eta * (1 - cosi * cosi);
            Vec3f refrdir = raydir * eta + nHit * (eta * cosi - sqrt(k));
            refrdir.normalize();
            // Trace refraction ray
            refraction = trace(pHit - nHit * bias, refrdir, spheres, depth + 1);
        }
        
        // Create surface color based on refraction/reflection rays
        sColor = ((reflection * 0.7) * 0.5) + ((refraction * (1/fresneleffect) * 0.8 * sphere->transparency) * 0.5);
        
    }
    
    for(int i = 0; i < spheres.size(); ++i) {
        if(spheres[i].emissionColor.x > 0) {
            // If we have a sphere that isn't reflecting or refracting, then shoot rays from sphere to light
            Vec3f trans = 1;
            // Calculate direction towards light
            Vec3f lightDir = spheres[i].center - pHit;
            lightDir.normalize();
            // For each sphere that isn't a light
            for(int j = 0; j < spheres.size(); ++j) {
                if(i != j) {
                    float t0, t1;
                    // Check to see if the ray towards the light intersects any sphere
                    if(spheres[j].intersect(pHit + nHit * bias, lightDir, t0, t1)) {
                        // if so, then scale color to a factor of 0.3
                        trans = 0.3;
                        break;
                    }
                }
            }
            sColor += (sphere->surfaceColor) * (trans) * std::max(float(0), nHit.dot(lightDir)) * (spheres[i].emissionColor);
            Vec3f h = (-raydir) + lightDir;
            Vec3f g;
            if(h.x < 0) g.x = h.x * -1;
            else g.x = h.x;
            if(h.y < 0) g.y = h.y * -1;
            else g.y = h.y;
            if(h.z < 0) g.z = h.z * -1;
            else g.z = h.z;
            Vec3f u;
            u.x = h.x / g.x;
            u.y = h.y / g.y;
            u.z = h.z / g.z;
            //printf("%f %f %f\n", u.x, u.y, u.z);
            sColor += 0.2 * std::max(float(0), nHit.dot(u));
        }
    }
    
    return sColor + sphere->emissionColor;
    
}

//TODO
void render(const std::vector<Sphere> &spheres)
{
    unsigned width = 1920, height = 1080;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 90, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    // Trace rays
    Vec3f pix;
    for (unsigned y = 0; y < height; ++y) {
        for (unsigned x = 0; x < width; ++x, ++pixel) {
            for(int i = -2; i < 3; i++) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx + (0.0001 * i), yy + (0.0001 * i), -1);
            raydir.normalize();
            pix = trace(Vec3f(0), raydir, spheres, 0);
            pix.x = pix.x / 5;
            pix.y = pix.y / 5;
            pix.z = pix.z / 5;
            *pixel += pix;
            }
        }
    }
    
    // Save result to a PPM image
    std::ofstream ofs("./tracedimage.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
        (unsigned char)(std::min(float(1), image[i].y) * 255) <<
        (unsigned char)(std::min(float(1), image[i].z) * 255);
    }
    ofs.close();
    delete [] image;
}


// Create Spheres and render
int main(int argc, char **argv)
{
    srand48(13);
    std::vector<Sphere> spheres;
    // position, radius, surface color, reflectivity, transparency, emission color
    spheres.push_back(Sphere(Vec3f( 0.0, -5, -20), 2, Vec3f(0.30, 0.30, 0.30), 0.01, 0.1));
    spheres.push_back(Sphere(Vec3f( 4.0, -5, -10), 2, Vec3f(0.30, 0.30, 0.30), 0.0, 0.0));
    spheres.push_back(Sphere(Vec3f( 0.0, -20007, -20), 20000, Vec3f(0.5, 0.2, 0.3), 0.0, 0.0));
    spheres.push_back(Sphere(Vec3f( 0.0, 0.0, -20060), 20000, Vec3f(0.10, 0.30, 0.10), 0.0, 0.0));
    // light
    spheres.push_back(Sphere(Vec3f(0, 60, 10), 3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(2)));
    render(spheres);
    
    return 0;
}
