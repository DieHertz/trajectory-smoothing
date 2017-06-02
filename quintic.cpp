#include <gsl/gsl_integration.h>
#include <nlopt.hpp>
#include <vector>
#include <iostream>
#include <fstream>

struct vec2 {
    double x, y;

    vec2 operator-() const {
        return { -x, -y };
    }
};

vec2 operator+(const vec2& a, const vec2& b) {
    return { a.x + b.x, a.y + b.y };
}

vec2 operator-(const vec2& a, const vec2& b) {
    return a + -b;
}

vec2 operator*(const vec2& a, const double k) {
    return { a.x * k, a.y * k };
}

vec2 operator*(const double k, const vec2& a) {
    return a * k;
}

//  Good ole' dot product
double dot(const vec2& a, const vec2& b) {
    return a.x * b.x + a.y * b.y;
}

//  Module of cross product, let's call it cross for convenience
double cross(const vec2& a, const vec2& b) {
    return a.x * b.y - a.y * b.x;
}

struct curve_vertex {
    vec2 pos;
    vec2 tangent;
    vec2 curvature;
    double t;
};

struct curve_segment {
    double magnitude_start;
    double magnitude_end;
    double alpha_start;
    double alpha_end;
};

std::vector<vec2> originalVertices{
    { 0, 0 },
    { 1, 1 },
    { 2, 2 },
    { 3, 2 }
};

std::vector<curve_vertex> vertices{
    { { 0, 0 }, { 1, 1 }, { 0, 0 }, 0 },
    { { 1, 1 }, { 1, 1 }, { 0, 0 }, 1 },
    { { 2, 2 }, { 1, 1 }, { 0, 0 }, 2 },
    { { 3, 2 }, { 1, 0 }, { 0, 0 }, 3 }
};

std::vector<curve_segment> segments{
    { 1, 1, 1, 1 },
    { 1, 1, 1, 1 },
    { 1, 1, 1, 1 }
};

void packArguments(double* out) {
    for (size_t i = 1; i < vertices.size() - 1; ++i) {
        *out++ = vertices[i].pos.x;
        *out++ = vertices[i].pos.y;
    }

    for (auto&& v : vertices) {
        *out++ = v.tangent.x;
        *out++ = v.tangent.y;
        *out++ = v.curvature.x;
        *out++ = v.curvature.y;
    }

    for (auto&& s : segments) {
        *out++ = s.magnitude_start;
        *out++ = s.magnitude_end;
        *out++ = s.alpha_start;
        *out++ = s.alpha_end;
    }
}

void unpackArguments(const double* out) {
    for (size_t i = 1; i < vertices.size() - 1; ++i) {
        vertices[i].pos.x = *out++;
        vertices[i].pos.y = *out++;
    }

    for (auto&& v : vertices) {
       v.tangent.x = *out++;
       v.tangent.y = *out++;
       v.curvature.x = *out++;
       v.curvature.y = *out++;
    }

    for (auto&& s : segments) {
        s.magnitude_start = *out++;
        s.magnitude_end = *out++;
        s.alpha_start = *out++;
        s.alpha_end = *out++;
    }
}

double hermiteQuintic0(const double t) {
    // 1 - 10t^3 + 15t^4 - 6t^5
    return 1 + t * t * t * (-10 + t * (15 - 6 * t));
}

double hermiteQuintic1(const double t) {
    // t - 6t^3 + 8t^4 - 3t^5
    return t * (1 + t * t * (-6 + t * (8 - 3 * t)));
}

double hermiteQuintic2(const double t) {
    // 0.5t^2 - 1.5t^3 + 1.5t^4 - 0.5t^5
    return t * t * (0.5 + t * (-1.5 + t * (1.5 - 0.5 * t)));
}

double hermiteQuintic3(const double t) {
    // 0.5t^3 - t^4 + 0.5t^5
    return t * t * t * (0.5 + t * (-1 + 0.5 * t));
}

double hermiteQuintic4(const double t) {
    // -4t^3 + 7t^4 - 3t^5
    return t * t * t * (-4 + t * (7 - 3 * t));
}

double hermiteQuintic5(const double t) {
    // 10t^3 - 15t^4 + 6t^5
    return t * t * t * (10 + t * (-15 + 6 * t));
}

// dC/dt
double dhermiteQuintic0(const double t) {
    // -30t^2 + 60t^3 - 30t^4
    return t * t * (-30 + t * (60 - 30 * t));
}

double dhermiteQuintic1(const double t) {
    // 1 - 18t^2 + 32t^3 - 15t^4
    return 1 + t * t * (-18 + t * (32 - 15 * t));
}

double dhermiteQuintic2(const double t) {
    // t - 4.5t^2 + 6t^3 - 2.5t^4
    return t * (1 + t * (-4.5 + t * (6 - 2.5 * t)));
}

double dhermiteQuintic3(const double t) {
    // 1.5t^2 - 4t^3 + 2.5t^4
    return t * t * (1.5 + t * (-4 + 2.5 * t));
}

double dhermiteQuintic4(const double t) {
    // -12t^2 + 28t^3 - 15t^4
    return t * t * (-12 + t * (28 - 15 * t));
}

double dhermiteQuintic5(const double t) {
    // 30t^2 - 60t^3 + 30t^4
    return t * t * (30 + t * (-60 + 30 * t));
}

// d^2C/dt^2
double d2hermiteQuintic0(const double t) {
    // -60t + 180t^2 - 120t^3
    return t * (-60 + t * (180 - 120 * t));
}

double d2hermiteQuintic1(const double t) {
    // -36t + 96t^2 - 60t^3
    return t * (-36 + t * (96 - 60 * t));
}

double d2hermiteQuintic2(const double t) {
    // 1 - 9t + 18t^2 - 10t^3
    return 1 + t * (-9 + t * (18 - 10 * t));
}

double d2hermiteQuintic3(const double t) {
    // 3t - 12t^2 + 10t^3
    return t * (3 + t * (-12 + 10 * t));
}

double d2hermiteQuintic4(const double t) {
    // -24t + 84t^2 - 60t^3
    return t * (-24 + t * (84 - 60 * t));
}

double d2hermiteQuintic5(const double t) {
    // 60t - 180t^2 + 120t^3
    return t * (60 + t * (-180 + 120 * t));
}

// d^3C/dt^3
double d3hermiteQuintic0(const double t) {
    // -60 + 360t - 360t^2
    return -60 + t * (360 - 360 * t);
}

double d3hermiteQuintic1(const double t) {
    // -36 + 192t - 180t^2
    return -36 + t * (192 - 180 * t);
}

double d3hermiteQuintic2(const double t) {
    // -9 + 36t - 30t^2
    return -9 + t * (36 - 30 * t);
}

double d3hermiteQuintic3(const double t) {
    // 3 - 24t + 30t^2
    return 3 + t * (-24 + 30 * t);
}

double d3hermiteQuintic4(const double t) {
    // -24 + 168t - 180t^2
    return -24 + t * (168 - 180 * t);
}

double d3hermiteQuintic5(const double t) {
    // 60 - 360t + 360^2
    return 60 + t * (-360 + 360 * t);
}

using fn = double(double);

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 spline(const curve_vertex& start, const curve_vertex& end, const double t) {
    return h0(t) * start.pos + h1(t) * start.tangent + h2(t) * start.curvature + h3(t) * end.curvature + h4(t) * end.tangent + h5(t) * end.pos;
}

vec2 hermiteQuintic(const curve_vertex& start, const curve_vertex& end, const double t) {
    return spline<hermiteQuintic0, hermiteQuintic1, hermiteQuintic2, hermiteQuintic3, hermiteQuintic4, hermiteQuintic5>(start, end, t);
}

vec2 dhermiteQuintic(const curve_vertex& start, const curve_vertex& end, const double t) {
    return spline<dhermiteQuintic0, dhermiteQuintic1, dhermiteQuintic2, dhermiteQuintic3, dhermiteQuintic4, dhermiteQuintic5>(start, end, t);
}

vec2 d2hermiteQuintic(const curve_vertex& start, const curve_vertex& end, const double t) {
    return spline<d2hermiteQuintic0, d2hermiteQuintic1, d2hermiteQuintic2, d2hermiteQuintic3, d2hermiteQuintic4, d2hermiteQuintic5>(start, end, t);
}

vec2 d3hermiteQuintic(const curve_vertex& start, const curve_vertex& end, const double t) {
    return spline<d3hermiteQuintic0, d3hermiteQuintic1, d3hermiteQuintic2, d3hermiteQuintic3, d3hermiteQuintic4, d3hermiteQuintic5>(start, end, t);
}

double curvature(const curve_vertex& start, const curve_vertex& end, const double t) {
    const auto dcurve = dhermiteQuintic(start, end, t);
    const auto d2curve = d2hermiteQuintic(start, end, t);

    return cross(dcurve, d2curve) / std::pow(dot(dcurve, dcurve), 1.5);
}

double dcurvature(const curve_vertex& start, const curve_vertex& end, const double t) {
    const auto dcurve = dhermiteQuintic(start, end, t);
    const auto d2curve = d2hermiteQuintic(start, end, t);
    const auto d3curve = d3hermiteQuintic(start, end, t);

    const auto velocitySq = dot(dcurve, dcurve);

    return (cross(dcurve, d3curve) * std::pow(velocitySq, 1.5) - 3 * cross(dcurve, d2curve) * std::sqrt(velocitySq) * dot(dcurve, d2curve)) / std::pow(velocitySq, 3);
}

const auto DISTANCE_WEIGHT = 1.0;
const auto VARIATION_WEIGHT = 1.0;

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dposx_start(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return { h0(t), 0 };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dposy_start(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return { 0, h0(t) };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dtangentx_start(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return {
        (h1(t) + segment.alpha_start * h2(t)) * segment.magnitude_start * segment.magnitude_start,
        0
    };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dtangenty_start(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return {
        0,
        (h1(t) + segment.alpha_start * h2(t)) * segment.magnitude_start * segment.magnitude_start
    };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dcurvaturex_start(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return {
        h2(t) * segment.magnitude_start * segment.magnitude_start * segment.magnitude_start * segment.magnitude_start,
        0
    };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dcurvaturey_start(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return {
        0,
        h2(t) * segment.magnitude_start * segment.magnitude_start * segment.magnitude_start * segment.magnitude_start
    };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dposx_end(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return { h5(t), 0 };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dposy_end(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return { 0, h5(t) };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dtangentx_end(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return {
        (h4(t) + segment.alpha_end * h3(t)) * segment.magnitude_end * segment.magnitude_end,
        0
    };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dtangenty_end(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return {
        0,
        (h4(t) + segment.alpha_end * h3(t)) * segment.magnitude_end * segment.magnitude_end
    };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dcurvaturex_end(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return {
        h3(t) * segment.magnitude_end * segment.magnitude_end * segment.magnitude_end * segment.magnitude_end,
        0
    };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dcurvaturey_end(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return {
        0,
        h3(t) * segment.magnitude_end * segment.magnitude_end * segment.magnitude_end * segment.magnitude_end
    };
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dmagnitude_start(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return h1(t) * 2 * start.tangent * segment.magnitude_start +
        h2(t) * (4 * start.curvature * segment.magnitude_start * segment.magnitude_start * segment.magnitude_start + 2 * start.tangent * segment.magnitude_start * segment.alpha_start);
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dmagnitude_end(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return h4(t) * 2 * end.tangent * segment.magnitude_end +
        h3(t) * (4 * end.curvature * segment.magnitude_end * segment.magnitude_end * segment.magnitude_end + 2 * end.tangent * segment.magnitude_end * segment.alpha_end);
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dalpha_start(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return h2(t) * segment.magnitude_start * segment.magnitude_start * start.tangent;
}

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 df_dalpha_end(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end, const double t) {
    return h3(t) * segment.magnitude_end * segment.magnitude_end * end.tangent;
}

using partial_derivative_f = vec2(const curve_vertex&, const curve_segment&, const curve_vertex&, const double);

template<partial_derivative_f f, partial_derivative_f df, partial_derivative_f d2f, partial_derivative_f d3f>
double integrate(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end) {
    return 2*(((double) rand()) / RAND_MAX - 0.5);
}

const auto INTEGRATION_WORKSPACE_SIZE = 1000;
const auto workspace = gsl_integration_workspace_alloc(INTEGRATION_WORKSPACE_SIZE);

double objective(unsigned n, const double* x, double* grad, void*) {
    unpackArguments(x);

    // Distance objective, only applied for inner points
    auto distanceSq = 0.;
    for (size_t i = 1; i < vertices.size() - 1; ++i) {
        distanceSq += dot(originalVertices[i] - vertices[i].pos, originalVertices[i] - vertices[i].pos);
    }

    // Variation constraint integral
    // @todo replace euler integration with something smarter
    auto variationIntegralSum = 0.;

    if (grad) {
        // for all points
        // integrate dF/dtangentx over segment to the left and segment to the right
        // integrate dF/dtangenty over segment to the left and segment to the right
        // integrate dF/dcurvaturex over segment to the left and segment to the right
        // integrate dF/dcurvaturey over segment to the left and segment to the right
        *grad++ = integrate<dcurve_dtangentx_start, ddcurve_dtangentx_start, dd2curve_dtangentx_start, dd3curve_dtangentx_start>(vertices.front(), segments.front(), vertices[1]);
        *grad++ = integrate<dcurve_dtangenty_start, ddcurve_dtangenty_start, dd2curve_dtangenty_start, dd3curve_dtangenty_start>(vertices.front(), segments.front(), vertices[1]);
        *grad++ = integrate<dcurve_dcurvaturex_start, ddcurve_dcurvaturex_start, dd2curve_dcurvaturex_start, dd3curve_dcurvaturex_start>(vertices.front(), segments.front(), vertices[1]);
        *grad++ = integrate<dcurve_dcurvaturey_start, ddcurve_dcurvaturey_start, dd2curve_dcurvaturey_start, dd3curve_dcurvaturey_start>(vertices.front(), segments.front(), vertices[1]);

        for (size_t i = 1; i < vertices.size() - 1; ++i) {
            *grad++ = integrate<dcurve_dtangentx_end, ddcurve_dtangentx_end, dd2curve_dtangentx_end, dd3curve_dtangentx_end>(vertices[i - 1], segments[i - 1], vertices[i]) + integrate<dcurve_dtangentx_start, ddcurve_dtangentx_start, dd2curve_dtangentx_start, dd3curve_dtangentx_start>(vertices[i], segments[i], vertices[i + 1]);
            *grad++ = integrate<dcurve_dtangenty_end, ddcurve_dtangenty_end, dd2curve_dtangenty_end, dd3curve_dtangenty_end>(vertices[i - 1], segments[i - 1], vertices[i]) + integrate<dcurve_dtangenty_start, ddcurve_dtangenty_start, dd2curve_dtangenty_start, dd3curve_dtangenty_start>(vertices[i], segments[i], vertices[i + 1]);
            *grad++ = integrate<dcurve_dcurvaturex_end, ddcurve_dcurvaturex_end, dd2curve_dcurvaturex_end, dd3curve_dcurvaturex_end>(vertices[i - 1], segments[i - 1], vertices[i]) + integrate<dcurve_dcurvaturex_start, ddcurve_dcurvaturex_start, dd2curve_dcurvaturex_start, dd3curve_dcurvaturex_start>(vertices[i], segments[i], vertices[i + 1]);
            *grad++ = integrate<dcurve_dcurvaturey_end, ddcurve_dcurvaturey_end, dd2curve_dcurvaturey_end, dd3curve_dcurvaturey_end>(vertices[i - 1], segments[i - 1], vertices[i]) + integrate<dcurve_dcurvaturey_start, ddcurve_dcurvaturey_start, dd2curve_dcurvaturey_start, dd3curve_dcurvaturey_start>(vertices[i], segments[i], vertices[i + 1]);
        }

        *grad++ = integrate<dcurve_dtangentx_end, ddcurve_dtangentx_end, dd2curve_dtangentx_end, dd3curve_dtangentx_end>(vertices[vertices.size() - 1], segments.back(), vertices.back());
        *grad++ = integrate<dcurve_dtangenty_end, ddcurve_dtangenty_end, dd2curve_dtangenty_end, dd3curve_dtangenty_end>(vertices[vertices.size() - 1], segments.back(), vertices.back());
        *grad++ = integrate<dcurve_dcurvaturex_end, ddcurve_dcurvaturex_end, dd2curve_dcurvaturex_end, dd3curve_dcurvaturex_end>(vertices[vertices.size() - 1], segments.back(), vertices.back());
        *grad++ = integrate<dcurve_dcurvaturey_end, ddcurve_dcurvaturey_end, dd2curve_dcurvaturey_end, dd3curve_dcurvaturey_end>(vertices[vertices.size() - 1], segments.back(), vertices.back());

        // for all segments
        // integrate dF/dmagnitude_start over current segment
        // integrate dF/dmagnitude_end over current segment
        // integrate dF/dalpha_start over current segment
        // integrate dF/dalpha_end over current segment
        for (size_t i = 0; i < segments.size(); ++i) {
            *grad++ = integrate<dcurve_dmagnitude_start, ddcurve_dmagnitude_start, dd2curve_dmagnitude_start, dd3curve_dmagnitude_start>(vertices[i], segments[i], vertices[i + 1]);
            *grad++ = integrate<dcurve_dmagnitude_end, ddcurve_dmagnitude_end, dd2curve_dmagnitude_end, dd3curve_dmagnitude_end>(vertices[i], segments[i], vertices[i + 1]);
            *grad++ = integrate<dcurve_dalpha_start, ddcurve_dalpha_start, dd2curve_dalpha_start, dd3curve_dalpha_start>(vertices[i], segments[i], vertices[i + 1]);
            *grad++ = integrate<dcurve_dalpha_end, ddcurve_dalpha_end, dd2curve_dalpha_end, dd3curve_dalpha_end>(vertices[i], segments[i], vertices[i + 1]);
        }
    }

    for (size_t i = 0; i < segments.size(); ++i) {
        auto& seg = segments[i];
        auto& v = vertices[i];
        auto& nv = vertices[i + 1];

        const auto magStart = seg.magnitude_start * seg.magnitude_start;
        const auto magEnd = seg.magnitude_end * seg.magnitude_end;

        const curve_vertex start{
            v.pos,
            v.tangent * magStart,
            v.curvature * magStart * magStart + v.tangent * seg.alpha_start * magStart
        };

        const curve_vertex end{
            nv.pos,
            nv.tangent * magEnd,
            nv.curvature * magEnd * magEnd + nv.tangent * seg.alpha_end * magEnd,
        };

        using vertex_pair = std::pair<const curve_vertex&, const curve_vertex&>;
        vertex_pair vertices{start, end};

        gsl_function integrand{
            [] (double t, void* params) {
                auto& vertices = *static_cast<vertex_pair*>(params);

                const auto dc = dcurvature(vertices.first, vertices.second, t);
                const auto speed = dhermiteQuintic(vertices.first, vertices.second, t);

                return dc * dc / std::sqrt(dot(speed, speed));
            },
            &vertices
        };

        double result, absError;
        gsl_integration_qag(&integrand, 0., 1., 1e-3, 1e-3, INTEGRATION_WORKSPACE_SIZE, GSL_INTEG_GAUSS15, workspace, &result, &absError);

        variationIntegralSum += result;
    }

    const auto objective = DISTANCE_WEIGHT * distanceSq + VARIATION_WEIGHT * variationIntegralSum;

    return objective;
}

int main() {
    const unsigned int taskSize = segments.size() * 4 + vertices.size() * 4 + (vertices.size() - 2) * 2;
    nlopt::opt optimizer(nlopt::LD_SLSQP/*LN_COBYLA*/, taskSize);
    optimizer.set_min_objective(objective, nullptr);
    optimizer.set_xtol_rel(1e-5);
    optimizer.set_ftol_abs(1e-5);
    optimizer.set_ftol_rel(1e-5);

    std::vector<double> x(taskSize);
    packArguments(&x[0]);

    double objectiveVal;
    try {
        optimizer.optimize(x.data(), x.size(), objectiveVal);
    } catch (const std::exception& exception) {
        std::cout << "e.what(): " << exception.what() << std::endl;
    }

    std::cout << "objective=" << objectiveVal << std::endl;
    std::copy(std::begin(x), std::end(x), std::ostream_iterator<double>{std::cout, " "});
    std::cout << std::endl;

    std::ofstream file{"spline.js"};

    file << "var smoothedSpline = [\n";
    for (size_t i = 0; i < segments.size(); ++i) {
        auto& seg = segments[i];
        auto& v = vertices[i];
        auto& nv = vertices[i + 1];

        const auto magStart = seg.magnitude_start * seg.magnitude_start;
        const auto magEnd = seg.magnitude_end * seg.magnitude_end;

        const curve_vertex start{
            v.pos,
            v.tangent * magStart,
            v.curvature * magStart * magStart + v.tangent * seg.alpha_start * magStart
        };

        const curve_vertex end{
            nv.pos,
            nv.tangent * magEnd,
            nv.curvature * magEnd * magEnd + nv.tangent * seg.alpha_end * magEnd,
        };

        const auto dt = 5.0e-3;
        for (double t = 0; t <= 1.0; t += dt) {
            auto&& curvePoint = hermiteQuintic(start, end, t);
            file << curvePoint.x << ", " << curvePoint.y << ", " << i + t << ",\n";
        }
    }
    file << "];" << std::endl;

    file << "var originalSpline = [\n";
    for (auto&& vertex : originalVertices) {
        file << vertex.x << ", " << vertex.y << ", " << vertex.x << ",\n";
    }
    file << "];" << std::endl;

    file << "var curvaturePlot = [\n";
    for (size_t i = 0; i < segments.size(); ++i) {
        auto& seg = segments[i];
        auto& v = vertices[i];
        auto& nv = vertices[i + 1];

        const auto magStart = seg.magnitude_start * seg.magnitude_start;
        const auto magEnd = seg.magnitude_end * seg.magnitude_end;

        const curve_vertex start{
            v.pos,
            v.tangent * magStart,
            v.curvature * magStart * magStart + v.tangent * seg.alpha_start * magStart
        };

        const curve_vertex end{
            nv.pos,
            nv.tangent * magEnd,
            nv.curvature * magEnd * magEnd + nv.tangent * seg.alpha_end * magEnd,
        };

        const auto dt = 5.0e-3;
        for (double t = 0; t <= 1.0; t += dt) {
            auto&& c = curvature(start, end, t);
            file << i + t << ", " << c << ", " << i + t << ",\n";
        }
    }
    file << "];" << std::endl;

    gsl_integration_workspace_free(workspace);
}

