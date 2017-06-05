#pragma once

#include <iostream>

template<typename T>
struct vec2_base {
    T x, y;

    vec2_base operator-() const {
        return { -x, -y };
    }
};

template<typename T>
vec2_base<T> operator+(const vec2_base<T>& a, const vec2_base<T>& b) {
    return { a.x + b.x, a.y + b.y };
}

template<typename T>
vec2_base<T> operator-(const vec2_base<T>& a, const vec2_base<T>& b) {
    return a + -b;
}

template<typename T>
vec2_base<T> operator*(const vec2_base<T>& a, const double k) {
    return { a.x * k, a.y * k };
}

template<typename T>
vec2_base<T> operator*(const double k, const vec2_base<T>& a) {
    return a * k;
}

template<typename T>
vec2_base<T> operator*(const vec2_base<T>& a, const T& k) {
    return { a.x * k, a.y * k };
}

template<typename T>
vec2_base<T> operator*(const T& k, const vec2_base<T>& a) {
    return a * k;
}

template<typename T>
vec2_base<T> operator/(const vec2_base<T>& a, const T& k) {
    return { a.x / k, a.y / k };
}

//  Good ole' dot product
template<typename T>
T dot(const vec2_base<T>& a, const vec2_base<T>& b) {
    return a.x * b.x + a.y * b.y;
}

//  Module of cross product, let's call it cross for convenience
template<typename T>
T cross(const vec2_base<T>& a, const vec2_base<T>& b) {
    return a.x * b.y - a.y * b.x;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const vec2_base<T>& v) {
    return os << "{ " << v.x << ", " << v.y << " }";
}
