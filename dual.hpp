#pragma once

#include <cmath>

enum class dual_value {
    f,
    df
};

struct dual {
    double x, dx;

    dual() = default;
    dual(double x) : x{x} {}
    dual(double x, double dx) : x{x}, dx{dx} {}

    dual operator-() const {
        return { -x, -dx };
    }
};

template<dual_value>
double get(const dual& d) { return d.x; }

template<>
double get<dual_value::df>(const dual& d) { return d.dx; }

dual operator+(const dual& a, const dual& b) {
    return { a.x + b.x, a.dx + b.dx };
}

dual operator-(const dual& a, const dual& b) {
    return a + -b;
}

dual operator*(const dual& a, const dual& b) {
    return { a.x * b.x, a.dx * b.x + a.x * b.dx };
}

dual operator*(const double a, const dual& b) {
    return dual{a, 0} * b;
}

dual operator*(const dual& a, const double b) {
    return b * a;
}

dual operator/(const dual& a, const dual& b) {
    return { a.x / b.x, (a.dx * b.x - a.x * b.dx) / (b.x * b.x) };
}

dual operator/(const dual& a, const double b) {
    return a / dual{b, 0};
}

dual operator/(const double a, const dual& b) {
    return dual{a, 0} / b;
}

dual sqrt(const dual& a) {
    const auto val = std::sqrt(a.x);
    return { val, a.dx / (2 * val) };
}

dual pow(const dual& a, double b) {
    const auto f = std::pow(a.x, b);
    return { f, f * b * a.dx / a.x };
}

dual pow(const dual& a, const dual& b) {
    const auto f = std::pow(a.x, b.x);
    return { f, f * (b.dx * std::log(a.x) + b.x * a.dx / a.x) };
}
