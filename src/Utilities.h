#pragma once
//#include <cstdint>
#include <cmath>
#include <string>

struct Color {
    double r, g, b;

    Color& operator+(const Color&);
    Color& operator*(const Color&);
    Color& operator/(const Color&);
    Color& operator/(unsigned int x);
};

struct Point {
    double x, y, z;

    Point& operator+(const Point&);
    Point& operator*(const Point&);
    bool operator==(const Point&) const;
};

struct Sphere {
    double r;
    double x0, y0, z0;
    double KsR, KsG, KsB;
    double KdR, KdG, KdB;
    double KaR, KaG, KaB;
    double n;
};

struct LightSource {
    double xs, ys, zs;
    double IsR, IsG, IsB;
    double IdR, IdG, IdB;
    double IaR, IaG, IaB;
};

struct Parameters {
    std::size_t x_input_file_resolution, y_input_file_resolution, ray_tracing_reflections;
    std::string input_file_path, output_file_name;
    bool help_mode = false;
    bool is_input_set = false;
};