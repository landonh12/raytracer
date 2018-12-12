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
Vec3f trace(const Vec3f &rayorig, const Vec3f &raydir, const std::vector<Sphere> &spheres, const std::vector<Triangle> &triangles, const int &depth) {
    
    float dist = INFINITY;       // Distance to intersection
    float distT = INFINITY;
    const Sphere *sphere = NULL; // Create sphere object
    const Triangle *triangle = NULL;
    float t0;                    // t0 sphere intersection
    float t1;                    // t1 sphere intersection
    float t;
    Vec3f tN;
    bool isTriangle;
    
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
    
    for(int i = 0; i < triangles.size(); i++) {
        t = INFINITY;
        if(triangles[i].intersect(rayorig, raydir, t, tN)) {
            printf("Tint!");
            if(t < distT) {
                distT = t;
                triangle = &triangles[i];
            }
        }
    }
    
    // If there's no cube intersections, then there's no hit!
    if(!sphere && !triangle) {
        return Vec3f(0.0, 0.0, 0.0);
    }
    
    if(triangle && !sphere)
        isTriangle = true;
    
    if(sphere && !triangle)
        isTriangle = false;
    
    if(sphere && triangle) {
        //printf("Sphere & Triangle!");
        if(t < t0)
            isTriangle = true;
        else
            isTriangle = false;
    }
    
    
    if(isTriangle) {
        //printf("tN: %f %f %f\n", tN.x, tN.y, tN.z);
        Vec3f sColor = Vec3f(0.0, 0.0, 0.0);
        Vec3f pHit = rayorig + raydir * distT;
        //printf("pHit: %f %f %f\n", pHit.x, pHit.y, pHit.z);
        //tN.normalize();
        
        float bias = 1e-2;
        
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
                        if(spheres[j].intersect(pHit, lightDir, t0, t1)) {
                            // if so, then scale color to a factor of 0.3
                            trans = 0.3;
                            break;
                        }
                    }
                }
                printf("tN: %f, %f, %f\n", tN.x, tN.y, tN.z);
                printf("tN.dot(lightDir): %f\n", tN.dot(lightDir));
                sColor = (triangle->surfaceColor) * tN.dot(lightDir);
            }
        }
        
        //printf("sColor: %f %f %f\n", sColor.x, sColor.y, sColor.z);
        return sColor;
        
    } else {
        Vec3f sColor = Vec3f(0.0, 0.0, 0.0);  // Surface color
        Vec3f pHit = rayorig + raydir * dist; // Point of intersection with the sphere
        Vec3f nHit;
        if(sphere) nHit = pHit - sphere->center;   // Normal of the hit
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
                reflection = trace(pHit + nHit * bias, refldir, spheres, triangles, depth + 1);
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
                refraction = trace(pHit - nHit * bias, refrdir, spheres, triangles, depth + 1);
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
            }
        }
        
        return sColor + sphere->emissionColor;
    }
    
}

//TODO
void render(const std::vector<Sphere> &spheres, const std::vector<Triangle> &triangles)
{
    unsigned width = 1920, height = 1080;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    // Trace rays
    Vec3f pix;
    for (unsigned y = 0; y < height; ++y) {
        for (unsigned x = 0; x < width; ++x, ++pixel) {
            for(int i = 0; i < 1; i++) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx + (0.0001 * i), yy + (0.0001 * i), -1);
            raydir.normalize();
            pix = trace(Vec3f(0), raydir, spheres, triangles, 0);
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
    std::vector<Triangle> triangles;
    // position, radius, surface color, reflectivity, transparency, emission color
    spheres.push_back(Sphere(Vec3f( 0.0, -5, -40), 2, Vec3f(0.30, 0.30, 0.30), 0.0, 0.0));
    spheres.push_back(Sphere(Vec3f( 6.0, -5, -30), 2, Vec3f(0.30, 0.30, 0.30), 0.1, 0.0));
    triangles.push_back(Triangle(Vec3f(-7, -7, -30), Vec3f(-7, -4, -30), Vec3f(-4, -7, -30), Vec3f(0.0, 1.0, 0.0)));
    triangles.push_back(Triangle(Vec3f(-4, -4, -30), Vec3f(-7, -4, -30), Vec3f(-4, -7, -30), Vec3f(0.0, 0.0, 1.0)));
    triangles.push_back(Triangle(Vec3f(-4, -4, -30), Vec3f(-4, -7, -30), Vec3f(-4, -7, -34), Vec3f(1.0, 0.0, 0.0)));
    triangles.push_back(Triangle(Vec3f(-4, -4, -34), Vec3f(-4, -7, -34), Vec3f(-4, -4, -30), Vec3f(0.0, 0.0, 1.0)));
    triangles.push_back(Triangle(Vec3f(-7, -4, -30), Vec3f(-7, -4, -34), Vec3f(-4, -4, -30), Vec3f(0.0, 1.0, 0.0)));
    triangles.push_back(Triangle(Vec3f(-4, -4, -34), Vec3f(-7, -4, -34), Vec3f(-4, -4, -30), Vec3f(1.0, 0.0, 0.0)));
    triangles.push_back(Triangle(Vec3f(-7, -7, -34), Vec3f(-7, -4, -34), Vec3f(-4, -7, -34), Vec3f(0.0, 1.0, 0.0)));
    triangles.push_back(Triangle(Vec3f(-4, -4, -34), Vec3f(-7, -4, -34), Vec3f(-4, -7, -34), Vec3f(0.0, 0.0, 1.0)));
    triangles.push_back(Triangle(Vec3f(-7, -4, -30), Vec3f(-7, -7, -30), Vec3f(-7, -7, -34), Vec3f(1.0, 0.0, 0.0)));
    triangles.push_back(Triangle(Vec3f(-7, -4, -34), Vec3f(-7, -7, -34), Vec3f(-7, -4, -30), Vec3f(0.0, 0.0, 1.0)));
    triangles.push_back(Triangle(Vec3f(-7, -7, -30), Vec3f(-7, -7, -34), Vec3f(-4, -7, -30), Vec3f(0.0, 1.0, 0.0)));
    triangles.push_back(Triangle(Vec3f(-4, -7, -34), Vec3f(-7, -7, -34), Vec3f(-4, -7, -30), Vec3f(1.0, 0.0, 0.0)));
    spheres.push_back(Sphere(Vec3f( 0.0, -20007, -20), 20000, Vec3f(0.5, 0.2, 0.3), 0.0, 0.0));
    spheres.push_back(Sphere(Vec3f( 0.0, 0.0, -20060), 20000, Vec3f(0.10, 0.30, 0.10), 0.0, 0.0));
    // light
    spheres.push_back(Sphere(Vec3f(0, 50, -20), 3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(10)));
    render(spheres, triangles);
    
    return 0;
}
