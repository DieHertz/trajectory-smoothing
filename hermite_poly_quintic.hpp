#pragma once

double h0(const double t) {
    // 1 - 10t^3 + 15t^4 - 6t^5
    return 1 + t * t * t * (-10 + t * (15 - 6 * t));
}

double h1(const double t) {
    // t - 6t^3 + 8t^4 - 3t^5
    return t * (1 + t * t * (-6 + t * (8 - 3 * t)));
}

double h2(const double t) {
    // 0.5t^2 - 1.5t^3 + 1.5t^4 - 0.5t^5
    return t * t * (0.5 + t * (-1.5 + t * (1.5 - 0.5 * t)));
}

double h3(const double t) {
    // 0.5t^3 - t^4 + 0.5t^5
    return t * t * t * (0.5 + t * (-1 + 0.5 * t));
}

double h4(const double t) {
    // -4t^3 + 7t^4 - 3t^5
    return t * t * t * (-4 + t * (7 - 3 * t));
}

double h5(const double t) {
    // 10t^3 - 15t^4 + 6t^5
    return t * t * t * (10 + t * (-15 + 6 * t));
}

// dC/dt
double dh0(const double t) {
    // -30t^2 + 60t^3 - 30t^4
    return t * t * (-30 + t * (60 - 30 * t));
}

double dh1(const double t) {
    // 1 - 18t^2 + 32t^3 - 15t^4
    return 1 + t * t * (-18 + t * (32 - 15 * t));
}

double dh2(const double t) {
    // t - 4.5t^2 + 6t^3 - 2.5t^4
    return t * (1 + t * (-4.5 + t * (6 - 2.5 * t)));
}

double dh3(const double t) {
    // 1.5t^2 - 4t^3 + 2.5t^4
    return t * t * (1.5 + t * (-4 + 2.5 * t));
}

double dh4(const double t) {
    // -12t^2 + 28t^3 - 15t^4
    return t * t * (-12 + t * (28 - 15 * t));
}

double dh5(const double t) {
    // 30t^2 - 60t^3 + 30t^4
    return t * t * (30 + t * (-60 + 30 * t));
}

// d^2C/dt^2
double d2h0(const double t) {
    // -60t + 180t^2 - 120t^3
    return t * (-60 + t * (180 - 120 * t));
}

double d2h1(const double t) {
    // -36t + 96t^2 - 60t^3
    return t * (-36 + t * (96 - 60 * t));
}

double d2h2(const double t) {
    // 1 - 9t + 18t^2 - 10t^3
    return 1 + t * (-9 + t * (18 - 10 * t));
}

double d2h3(const double t) {
    // 3t - 12t^2 + 10t^3
    return t * (3 + t * (-12 + 10 * t));
}

double d2h4(const double t) {
    // -24t + 84t^2 - 60t^3
    return t * (-24 + t * (84 - 60 * t));
}

double d2h5(const double t) {
    // 60t - 180t^2 + 120t^3
    return t * (60 + t * (-180 + 120 * t));
}

// d^3C/dt^3
double d3h0(const double t) {
    // -60 + 360t - 360t^2
    return -60 + t * (360 - 360 * t);
}

double d3h1(const double t) {
    // -36 + 192t - 180t^2
    return -36 + t * (192 - 180 * t);
}

double d3h2(const double t) {
    // -9 + 36t - 30t^2
    return -9 + t * (36 - 30 * t);
}

double d3h3(const double t) {
    // 3 - 24t + 30t^2
    return 3 + t * (-24 + 30 * t);
}

double d3h4(const double t) {
    // -24 + 168t - 180t^2
    return -24 + t * (168 - 180 * t);
}

double d3h5(const double t) {
    // 60 - 360t + 360^2
    return 60 + t * (-360 + 360 * t);
}