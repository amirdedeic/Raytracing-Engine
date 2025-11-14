#include <SDL2/SDL.h>
#include <cmath>
#include <vector>
#include <limits>
#include <string>

using namespace std;

// Constants
const int Cw = 900;
const int Ch = 900;
const float Vw = 1.5f;
const float Vh = 1.5f;
const float d = 1.0f;

// Vec3 struct
struct Vec3 {
    float x, y, z;
    
    Vec3 operator-(const Vec3& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }
    
    Vec3 operator*(float t) const {
        return {x * t, y * t, z * t};
    }
    
    Vec3 operator+(const Vec3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }
    
    Vec3 operator/(float t) const {
        return {x / t, y / t, z / t};
    }

    Vec3 operator-() const {
        return {-x, -y, -z};
    }
};

float dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

float length(const Vec3& v) {
    return sqrt(dot(v, v));
}

// Color struct
struct Color {
    int r, g, b;
    
    Color operator*(float intensity) const {
        return {
            min(255, max(0, (int)(r * intensity))),
            min(255, max(0, (int)(g * intensity))),
            min(255, max(0, (int)(b * intensity)))
        };
    }
};

const Color BACKGROUND_COLOR = {255, 255, 255};

// Sphere struct
struct Sphere {
    Vec3 center;
    float radius;
    Color color;
    float specular;
};

// Light struct
struct Light {
    string type;  // "ambient", "point", "directional"
    float intensity;
    Vec3 position;   // for point lights
    Vec3 direction;  // for directional lights
};

// Scene
vector<Sphere> spheres = {
    {{0, -1, 3}, 1, {255, 0, 0}, 500},      // Red
    {{2, 0, 4}, 1, {0, 0, 255}, 500},       // Blue
    {{-2, 0, 4}, 1, {0, 255, 0}, 300},      // Green
    {{0, -5001, 0}, 5000, {255, 255, 0}, 1000} // Yellow floor
};

// Lights
vector<Light> lights = {
    {"ambient", 0.2f, {0,0,0}, {0,0,0}},
    {"point", 0.6f, {-2, 1.5, 0}, {0,0,0}},
    {"directional", 0.2f, {0,0,0}, {1, -4, 4}}
};

Vec3 CanvasToViewport(int x, int y) {
    return {x * Vw / Cw, y * Vh / Ch, d};
}

pair<float, float> IntersectRaySphere(Vec3 O, Vec3 D, const Sphere& sphere) {
    float r = sphere.radius;
    Vec3 CO = O - sphere.center;
    
    float a = dot(D, D);
    float b = 2 * dot(CO, D);
    float c = dot(CO, CO) - r * r;
    
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {INFINITY, INFINITY};
    }
    
    float t1 = (-b + sqrt(discriminant)) / (2 * a);
    float t2 = (-b - sqrt(discriminant)) / (2 * a);
    return {t1, t2};
}

float ComputeLighting(Vec3 P, Vec3 N, Vec3 V, float s) {
    float i = 0.0f;
    
    //# check if there is any intersection with another primative

    for (auto& light : lights) {
        if (light.type == "ambient") {
            i += light.intensity;
        } else {
            Vec3 L;
            if (light.type == "point") {
                L = light.position - P;
            } else {  // directional
                L = light.direction;
            }
            
            // point
            float n_dot_l = dot(N, L);
            if (n_dot_l > 0) {
                i += light.intensity * n_dot_l / (length(N) * length(L));
            }

            // specular
            if (s != -1) {
                Vec3 R = N * (2.0f * dot(N, L)) - L;
                float r_dot_v = dot(R, V);

                if (r_dot_v > 0) {
                    i += light.intensity * pow(r_dot_v/(length(R) * length(V)), s);
                }
            }
        }
    }
    
    return i;
}

pair<Sphere*, float> ClosestIntersection(Vec3 O, Vec3 D, float t_min, float t_max) {
    float closest_t = INFINITY;
    Sphere* closest_sphere = nullptr;
    
    for (auto& sphere : spheres) {
        auto [t1, t2] = IntersectRaySphere(O, D, sphere);
        
        if (t1 >= t_min && t1 <= t_max && t1 < closest_t) {
            closest_t = t1;
            closest_sphere = &sphere;
        }
        if (t2 >= t_min && t2 <= t_max && t2 < closest_t) {
            closest_t = t2;
            closest_sphere = &sphere;
        }
    }
    
    return {closest_sphere, closest_t};
}

Color TraceRay(Vec3 O, Vec3 D, float t_min, float t_max) {
    float closest_t = INFINITY;
    Sphere* closest_sphere = nullptr;
    
    for (auto& sphere : spheres) {
        auto [t1, t2] = IntersectRaySphere(O, D, sphere);
        
        if (t1 >= t_min && t1 <= t_max && t1 < closest_t) {
            closest_t = t1;
            closest_sphere = &sphere;
        }
        if (t2 >= t_min && t2 <= t_max && t2 < closest_t) {
            closest_t = t2;
            closest_sphere = &sphere;
        }
    }
    
    if (closest_sphere == nullptr) {
        return BACKGROUND_COLOR; // White background
    }

    
    // Compute intersection point and normal
    Vec3 P = O + D * closest_t;
    Vec3 N = P - closest_sphere->center;
    N = N / length(N);
    
    // Compute lighting
    float intensity = ComputeLighting(P, N, -D, closest_sphere->specular);
    
    return closest_sphere->color * intensity;
}

void PutPixel(SDL_Renderer* renderer, int x, int y, Color color) {
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
    SDL_RenderDrawPoint(renderer, Cw/2 + x, Ch/2 - y);
}

int main(int argc, char* argv[]) {
    SDL_Init(SDL_INIT_VIDEO);
    
    SDL_Window* window = SDL_CreateWindow("Raytracer", 
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, Cw, Ch, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    
    // Clear to black
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    
    // Raytrace
    Vec3 O = {0, 0, 0};
    for (int x = -Cw/2; x < Cw/2; x++) {
        for (int y = -Ch/2; y < Ch/2; y++) {
            Vec3 D = CanvasToViewport(x, y);
            Color color = TraceRay(O, D, 1, INFINITY);
            PutPixel(renderer, x, y, color);
        }
    }
    
    SDL_RenderPresent(renderer);
    
    // Event loop
    bool running = true;
    SDL_Event event;
    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }
    }
    
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}