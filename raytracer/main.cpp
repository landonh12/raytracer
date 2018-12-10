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

#define MAX_RAY_DEPTH 5

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

//TODO
Vec3f trace(const Vec3f &rayorig, const Vec3f &raydir, const std::vector<Sphere> &spheres, const std::vector<Cube> &cubes, const int &depth) {
    
    float dist = INFINITY; // Dist to intersection
    const Sphere *sphere = NULL; // Create sphere object
    const Cube *cube = NULL;
    float t0;
    float t1;
    float tc;
    // Check each object for intersection
    for(int i = 0; i < spheres.size(); i++) {
        // t0, t1 are intersection points
        t0 = INFINITY;
        t1 = INFINITY;
        // Gives t0 and t1 back from intersect function
        if(spheres[i].intersect(rayorig, raydir, t0, t1)) {
            if(t0 < 0) {
                t0 = t1;
            }
            if(t0 < dist) {
                dist = t0;
                sphere = &spheres[i];
            }
        }
    }
    // If there's no intersection, then sphere will be NULL. Return a bogus Vec3f.
    if(!sphere){
        // Check each object for intersection
        for(int i = 0; i < cubes.size(); i++) {
            // t0, t1 are intersection points
            t0 = INFINITY;
            t1 = INFINITY;
            // Gives t0 and t1 back from intersect function
            if(cubes[i].intersect(rayorig, raydir, tc)) {
                    dist = tc;
                    cube = &cubes[i];
            }
        }
        if(!cube) {
            return Vec3f(0.0, 0.0, 0.0);
        }
    }
    Vec3f sColor = 0; // Surface color of the sphere
    Vec3f pHit = rayorig + raydir * dist; // Point of intersection with the sphere
    Vec3f nHit;
    if(!sphere) {
        nHit = pHit - cube->bounds[0];
    } else {
        nHit = pHit - sphere->center; // Normal of the surface of intersection
    }
    nHit.normalize(); // Normalize normal dir
    // START REFLECTION AND REFRACTION
    
    float bias = 1e-2;
    
    bool inside = false;
    if(raydir.dot(nHit) > 0) {
        nHit = -nHit;
        inside = true;
    }
    if(sphere && (sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
        float facingratio = -raydir.dot(nHit);
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
        Vec3f reflection = Vec3f(0.0,0.0,0.0);
        if(sphere->reflection) {
            Vec3f refldir = raydir - nHit * 2 * raydir.dot(nHit);
            refldir.normalize();
            reflection = trace(pHit + nHit * bias, refldir, spheres, cubes, depth + 1);
        }
        Vec3f refraction = 0;
        if(sphere->transparency) {
            float ior = 1.345, eta = (inside) ? ior : 1 / ior;
            float cosi = -nHit.dot(raydir);
            float k = 1 - eta * eta * (1 - cosi * cosi);
            Vec3f refrdir = raydir * eta + nHit * (eta * cosi - sqrt(k));
            refrdir.normalize();
            refraction = trace(pHit - nHit * bias, refrdir, spheres, cubes, depth + 1);
        }
        
        sColor = (reflection * 0.7) + (refraction * (1/fresneleffect) * 0.8 * sphere->transparency);
        
    } else {
    
        // LIGHTS
        if(sphere) {
            for(int i = 0; i < spheres.size(); ++i) {
                if(spheres[i].emissionColor.x > 0) {
                    Vec3f trans = 1;
                    Vec3f lightDir = spheres[i].center - pHit;
                    lightDir.normalize();
                    for(int j = 0; j < spheres.size(); ++j) {
                        if(i != j) {
                            float t0, t1;
                            if(spheres[j].intersect(pHit + nHit * bias, lightDir, t0, t1)) {
                                trans = 0.3;
                                break;
                            }
                        }
                    }
                    sColor += (sphere->surfaceColor) * (trans) * std::max(float(0), nHit.dot(lightDir)) * (spheres[i].emissionColor);
                }
            }
        } else if(cube) {
            for(int i = 0; i < spheres.size(); ++i) {
                if(spheres[i].emissionColor.x > 0) {
                    Vec3f trans = 1;
                    Vec3f lightDir = spheres[i].center - pHit;
                    lightDir.normalize();
                    for(int j = 0; j < cubes.size(); ++j) {
                        float t0;
                        if(cubes[j].intersect(pHit + nHit * bias, lightDir, t0)) {
                            trans = 0;
                            break;
                        }
                    }
                    printf("%f, ", nHit.x);
                    sColor += (cube->c) * (trans) * std::max(float(0), nHit.dot(lightDir)) * (spheres[i].emissionColor);
                }
            }
        }
    }
    
    if(sphere) {
        return sColor + sphere->emissionColor;
    } else {
        return cube->c;
    }
    
    
}

//TODO
void render(const std::vector<Sphere> &spheres, const std::vector<Cube> &cubes)
{
    unsigned width = 1920, height = 1080;
    Vec3f *image = new Vec3f[width * height], *pixel = image;
    float invWidth = 1 / float(width), invHeight = 1 / float(height);
    float fov = 30, aspectratio = width / float(height);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    // Trace rays
    for (unsigned y = 0; y < height; ++y) {
        for (unsigned x = 0; x < width; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -1);
            raydir.normalize();
            *pixel = trace(Vec3f(0), raydir, spheres, cubes, 0);
        }
    }
    
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
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
    std::vector<Cube> cubes;
    // position, radius, surface color, reflectivity, transparency, emission color
    spheres.push_back(Sphere(Vec3f( -6.0, 0.0, -30), 2, Vec3f(0.50, 0.50, 0.50), 0.0, 0.0));
    spheres.push_back(Sphere(Vec3f( 0.0, 0.0, -40), 2, Vec3f(0.50, 0.50, 0.50), 0.5, 0.0));
    spheres.push_back(Sphere(Vec3f( 4.0, 0.0, -23), 2, Vec3f(0.90, 0.90, 0.90), 0.0, 0.1));
    spheres.push_back(Sphere(Vec3f( 0.0, -20003, -20), 20000, Vec3f(0.5, 0.2, 0.3), 0, 0.0));
    spheres.push_back(Sphere(Vec3f( 0.0, 0.0, -20060), 20000, Vec3f(0.10, 0.30, 0.10), 0, 0.0));
    //cubes.push_back(Cube(Vec3f(-6.0, -3.0, -37), Vec3f(6.0, 3.0, -31), Vec3f(0.50, 0.0, 0.00)));
    // light
    spheres.push_back(Sphere(Vec3f(3.0,  40, -10),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(2)));
    render(spheres, cubes);
    
    return 0;
}
