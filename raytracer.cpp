#include <SDL2/SDL.h>
#include <cmath>
#include <vector>
#include <limits>
#include <string>
#include <chrono>
#include <iostream>

using namespace std;

// Constants
const int Cw = 280;
const int Ch = 280;
const float Vw = 1.5f;
const float Vh = 1.5f;
const float d = 2.0f;
const int RECURSION_DEPTH = 2;
const Uint8* keys = SDL_GetKeyboardState(NULL);

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

struct RotationMatrix {
    double m[3][3];
    
    RotationMatrix(double yaw, double pitch) {
        double cy = cos(yaw);
        double sy = sin(yaw);
        double cp = cos(pitch);
        double sp = sin(pitch);
        
        m[0][0] = cy;       m[0][1] = sy*sp;   m[0][2] = sy*cp;
        m[1][0] = 0;        m[1][1] = cp;      m[1][2] = -sp;
        m[2][0] = -sy;      m[2][1] = cy*sp;   m[2][2] = cy*cp;
    }
    
    Vec3 operator*(const Vec3& v) const { 
        return {
            m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
            m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
            m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z
        };
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

    Color operator+(const Color& other) const {
        return {
            min(255, r + other.r),
            min(255, g + other.g),
            min(255, b + other.b)
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
    float reflective;
};

// Light struct
struct Light {
    string type;
    float intensity;
    Vec3 position;
    Vec3 direction;
};

// Scene
vector<Sphere> spheres = {
    {{0, 0, 3}, 1, {255, 0, 0}, 200, 0.2},
    {{3, 0, 6}, 1, {0, 0, 255}, 500, 0.3},
    {{-3, 0, 6}, 1, {0, 255, 0}, 300, 0.4},
    {{0, -5001, 0}, 5000, {255, 255, 0}, 1000, 0.0}
};

// Lights
vector<Light> lights = {
    {"directional", 0.9f, {0,0,0}, {1, 4, -4}},
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

float ComputeLighting(Vec3 P, Vec3 N, Vec3 V, float s) {
    float i = 0.0f;

    for (auto& light : lights) {
        if (light.type == "ambient") {
            i += light.intensity;
        } else {
            Vec3 L;
            float t_max;
            if (light.type == "point") {
                L = light.position - P;
                t_max = 1;
            } else {
                L = light.direction;
                t_max = INFINITY;
            }

            auto [shadow_sphere, shadow_t] = ClosestIntersection(P, L, 0.001, t_max);
            if (shadow_sphere != nullptr) {
                continue;
            }

            float n_dot_l = dot(N, L);
            if (n_dot_l > 0) {
                i += light.intensity * n_dot_l / (length(N) * length(L));
            }

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

Vec3 ReflectRay(Vec3 R, Vec3 N) {
    return N * 2.0f * dot(N, R) - R;
}

Color TraceRay(Vec3 O, Vec3 D, float t_min, float t_max, int RECURSION_DEPTH) {
    auto [closest_sphere, closest_t] = ClosestIntersection(O, D, t_min, t_max);
    if (closest_sphere == nullptr) {
        return BACKGROUND_COLOR;
    }

    Vec3 P = O + D * closest_t;
    Vec3 N = P - closest_sphere->center;
    N = N / length(N);
    
    float intensity = ComputeLighting(P, N, -D, closest_sphere->specular);
    Color local_color = closest_sphere->color * intensity;

    float r = closest_sphere->reflective;
    if (RECURSION_DEPTH <= 0 || r <= 0) {
        return local_color;
    }

    Vec3 R = ReflectRay(-D, N);
    Color reflected_color = TraceRay(P, R, 0.001, INFINITY, RECURSION_DEPTH - 1);

    return local_color * (1 - r) + reflected_color * r;
}

void PutPixel(SDL_Renderer* renderer, int x, int y, Color color) {
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
    SDL_RenderDrawPoint(renderer, Cw/2 + x, Ch/2 - y);
}

int main() {
    SDL_Init(SDL_INIT_VIDEO);
    // omp_set_num_threads(8);
    
    SDL_Window* window = SDL_CreateWindow("Raytracer", 
    SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, Cw, Ch, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    
    auto start = chrono::high_resolution_clock::now();
    
    Vec3 camera_pos = {0, 0, -0.1};
    double yaw = 0.0;
    double pitch = 0.0;
    bool needs_redraw = true;
    bool running = true;
    SDL_Event event;

    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }
        if (keys[SDL_SCANCODE_W]) {
            Vec3 forward = {sin(yaw), 0, cos(yaw)};
            camera_pos = camera_pos + forward * 0.2;
            needs_redraw = true;
        }
        if (keys[SDL_SCANCODE_S]) {
            Vec3 forward = {sin(yaw), 0, cos(yaw)};
            camera_pos = camera_pos + forward * (-0.2);
            needs_redraw = true;
        }
        if (keys[SDL_SCANCODE_A]) {
            Vec3 right = {cos(yaw), 0, -sin(yaw)};
            camera_pos = camera_pos + right * (-0.2);
            needs_redraw = true;
        }
        if (keys[SDL_SCANCODE_D]) {
            Vec3 right = {cos(yaw), 0, -sin(yaw)};
            camera_pos = camera_pos + right * 0.2;
            needs_redraw = true;
        }
        if (keys[SDL_SCANCODE_UP]) {
            pitch -= 0.15;
            needs_redraw = true;
        }
        if (keys[SDL_SCANCODE_DOWN]) {
            pitch += 0.15;
            needs_redraw = true;
        }
        if (keys[SDL_SCANCODE_LEFT]) {
            yaw -= 0.15;
            needs_redraw = true;
        }
        if (keys[SDL_SCANCODE_RIGHT]) {
            yaw += 0.15;
            needs_redraw = true;
        }
    
        if (needs_redraw) {
            RotationMatrix rotation(yaw, pitch);
            vector<Color> pixels(Cw * Ch);
            
            #pragma omp parallel for collapse(2)
            for (int x = -Cw/2; x < Cw/2; x++) {
                for (int y = -Ch/2; y < Ch/2; y++) {
                    Vec3 D = rotation * CanvasToViewport(x, y);
                    Color color = TraceRay(camera_pos, D, 1, INFINITY, RECURSION_DEPTH);
                    int idx = (x + Cw/2) + (y + Ch/2) * Cw;
                    pixels[idx] = color;
                }
            }
            
            for (int x = -Cw/2; x < Cw/2; x++) {
                for (int y = -Ch/2; y < Ch/2; y++) {
                    int idx = (x + Cw/2) + (y + Ch/2) * Cw;
                    PutPixel(renderer, x, y, pixels[idx]);
                }
            }
            SDL_RenderPresent(renderer);
            needs_redraw = false;
            
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
            cout << "Render time: " << duration.count() << " ms" << endl;
            start = chrono::high_resolution_clock::now();
        }
    }
    
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}