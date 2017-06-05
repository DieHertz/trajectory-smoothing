#include "vec2.hpp"
#include "dual.hpp"
#include "hermite_poly_quintic.hpp"
#include <gsl/gsl_integration.h>
#include <nlopt.hpp>
#include <random>
#include <vector>
#include <iostream>
#include <fstream>

using vec2 = vec2_base<dual>;

struct curve_vertex {
    vec2 pos;
    vec2 tangent;
    vec2 curvature;
    double t;
};

struct curve_segment {
    dual magnitude_start;
    dual magnitude_end;
    dual alpha_start;
    dual alpha_end;
};

struct vertex {
    vec2 pos;
    double t;
};

std::vector<vertex> originalVertices{
    { { 0, 0 }, 0 },
    { { 1, 1 }, 1 },
    { { 2, 2 }, 2 },
    { { 3, 2 }, 3 }
};

std::vector<curve_vertex> vertices;
std::vector<curve_segment> segments;

void packParameters(double* out) {
    for (size_t i = 1; i < vertices.size() - 1; ++i) {
        *out++ = vertices[i].pos.x.x;
        *out++ = vertices[i].pos.y.x;
    }

    for (auto&& v : vertices) {
        *out++ = v.tangent.x.x;
        *out++ = v.tangent.y.x;
        *out++ = v.curvature.x.x;
        *out++ = v.curvature.y.x;
    }

    for (auto&& s : segments) {
        *out++ = s.magnitude_start.x;
        *out++ = s.magnitude_end.x;
        *out++ = s.alpha_start.x;
        *out++ = s.alpha_end.x;
    }
}

void unpackParameters(const double* out) {
    for (size_t i = 1; i < vertices.size() - 1; ++i) {
        vertices[i].pos.x.x = *out++;
        vertices[i].pos.y.x = *out++;
    }

    for (auto&& v : vertices) {
       v.tangent.x.x = *out++;
       v.tangent.y.x = *out++;
       v.curvature.x.x = *out++;
       v.curvature.y.x = *out++;
    }

    for (auto&& s : segments) {
        s.magnitude_start.x = *out++;
        s.magnitude_end.x = *out++;
        s.alpha_start.x = *out++;
        s.alpha_end.x = *out++;
    }
}

using fn = double(double);

template<fn h0, fn h1, fn h2, fn h3, fn h4, fn h5>
vec2 spline(const curve_vertex& start, const curve_vertex& end, const double t) {
    return h0(t) * start.pos + h1(t) * start.tangent + h2(t) * start.curvature + h3(t) * end.curvature + h4(t) * end.tangent + h5(t) * end.pos;
}

vec2 hermite(const curve_vertex& start, const curve_vertex& end, const double t) {
    return spline<h0, h1, h2, h3, h4, h5>(start, end, t);
}

vec2 dhermite(const curve_vertex& start, const curve_vertex& end, const double t) {
    return spline<dh0, dh1, dh2, dh3, dh4, dh5>(start, end, t);
}

vec2 d2hermite(const curve_vertex& start, const curve_vertex& end, const double t) {
    return spline<d2h0, d2h1, d2h2, d2h3, d2h4, d2h5>(start, end, t);
}

vec2 d3hermite(const curve_vertex& start, const curve_vertex& end, const double t) {
    return spline<d3h0, d3h1, d3h2, d3h3, d3h4, d3h5>(start, end, t);
}

dual curvature(const curve_vertex& start, const curve_vertex& end, const double t) {
    const auto dcurve = dhermite(start, end, t);
    const auto d2curve = d2hermite(start, end, t);

    return cross(dcurve, d2curve) / pow(dot(dcurve, dcurve), 1.5);
}

dual dcurvature(const curve_vertex& start, const curve_vertex& end, const double t) {
    const auto dcurve = dhermite(start, end, t);
    const auto d2curve = d2hermite(start, end, t);
    const auto d3curve = d3hermite(start, end, t);

    const auto velocitySq = dot(dcurve, dcurve);

    return (cross(dcurve, d3curve) * pow(velocitySq, 1.5) - 3 * cross(dcurve, d2curve) * sqrt(velocitySq) * dot(dcurve, d2curve)) / pow(velocitySq, 3);
}

const auto DISTANCE_WEIGHT = 1.0;
const auto VARIATION_WEIGHT = 1.0;

dual d(const dual& d) {
    return { d.x, 1 };
}

vec2 dx(const vec2& v) {
    return { d(v.x), v.y };
}

vec2 dy(const vec2& v) {
    return { v.x, d(v.y) };
}

curve_vertex dtangent_x(const curve_vertex& v) {
    return {
        v.pos,
        dx(v.tangent),
        v.curvature
    };
}

curve_vertex dtangent_y(const curve_vertex& v) {
    return {
        v.pos,
        dy(v.tangent),
        v.curvature
    };
}

curve_vertex dcurvature_x(const curve_vertex& v) {
    return {
        v.pos,
        v.tangent,
        dx(v.curvature)
    };
}

curve_vertex dcurvature_y(const curve_vertex& v) {
    return {
        v.pos,
        v.curvature,
        dy(v.tangent)
    };
}

curve_segment dmagnitude_start(const curve_segment& s) {
    return {
        d(s.magnitude_start),
        s.magnitude_end,
        s.alpha_start,
        s.alpha_end
    };
}

curve_segment dmagnitude_end(const curve_segment& s) {
    return {
        s.magnitude_start,
        d(s.magnitude_end),
        s.alpha_start,
        s.alpha_end
    };
}

curve_segment dalpha_start(const curve_segment& s) {
    return {
        s.magnitude_start,
        s.magnitude_end,
        d(s.alpha_start),
        s.alpha_end
    };
}

curve_segment dalpha_end(const curve_segment& s) {
    return {
        s.magnitude_start,
        s.magnitude_end,
        s.alpha_start,
        d(s.alpha_end)
    };
}

const auto INTEGRATION_WORKSPACE_SIZE = 1000;
const auto workspace = gsl_integration_workspace_alloc(INTEGRATION_WORKSPACE_SIZE);

template<dual_value which>
double integrate(const curve_vertex& start, const curve_segment& segment, const curve_vertex& end) {
    const auto magStart = segment.magnitude_start * segment.magnitude_start;
    const auto magEnd = segment.magnitude_end * segment.magnitude_end;

    using vertex_pair = std::pair<const curve_vertex&, const curve_vertex&>;
    vertex_pair vertices{
        {
            start.pos,
            start.tangent * magStart,
            start.curvature * magStart * magStart + start.tangent * segment.alpha_start * magStart
        },
        {
            end.pos,
            end.tangent * magEnd,
            end.curvature * magEnd * magEnd + end.tangent * segment.alpha_end * magEnd,
        }
    };

    /*gsl_function integrand{
        [] (double t, void* params) -> double {
            auto& vertices = *static_cast<vertex_pair*>(params);

            const auto dc = dcurvature(vertices.first, vertices.second, t);
            const auto speed = dhermite(vertices.first, vertices.second, t);
            const auto f = get<which>(dc * dc / sqrt(dot(speed, speed)));

            return f;
        },
        &vertices
    };

    double result, absError;
    size_t neval;

    gsl_integration_qag(&integrand, 0., 1., 1e-3, 1e-3, INTEGRATION_WORKSPACE_SIZE, GSL_INTEG_GAUSS15, workspace, &result, &absError);
    //gsl_integration_qng(&integrand, 0., 1., 1e-3, 10, &result, &absError, &neval);*/

    const auto steps = 1000;
    const auto dt = 1. / steps;
    auto result = 0.;
    for (auto i = 0; i < steps; ++i) {
        const auto t = i * dt;

        const auto dc = dcurvature(vertices.first, vertices.second, t);
        const auto speed = dhermite(vertices.first, vertices.second, t);
        const auto f = dc * dc / sqrt(dot(speed, speed));

        result += get<which>(f) * dt;
    }

    if (std::isnan(result)) {
        throw 1;
    }

    return result;
}

double objective(unsigned n, const double* x, double* grad, void*) {
    unpackParameters(x);

    // Distance objective, only applied for inner points
    auto distanceSq = 0.;
    for (size_t i = 1; i < vertices.size() - 1; ++i) {
        distanceSq += dot(vertices[i].pos - originalVertices[i].pos, vertices[i].pos - originalVertices[i].pos).x;
    }

    // Variation constraint integral
    // @todo replace euler integration with something smarter
    auto variationIntegralSum = 0.;

    if (grad) {
        // for movable points
        for (size_t i = 1; i < vertices.size() - 1; ++i) {
            *grad++ = DISTANCE_WEIGHT * dot(dx(vertices[i].pos) - originalVertices[i].pos, dx(vertices[i].pos) - originalVertices[i].pos).dx;
            *grad++ = DISTANCE_WEIGHT * dot(dy(vertices[i].pos) - originalVertices[i].pos, dy(vertices[i].pos) - originalVertices[i].pos).dx;
        }

        // for all points
        *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(dtangent_x(vertices.front()), segments.front(), vertices[1]);
        *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(dtangent_y(vertices.front()), segments.front(), vertices[1]);
        *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(dcurvature_x(vertices.front()), segments.front(), vertices[1]);
        *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(dcurvature_y(vertices.front()), segments.front(), vertices[1]);

        for (size_t i = 1; i < vertices.size() - 1; ++i) {
            *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[i - 1], segments[i - 1], dtangent_x(vertices[i])) +
                VARIATION_WEIGHT * integrate<dual_value::df>(dtangent_x(vertices[i]), segments[i], vertices[i + 1]);
            *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[i - 1], segments[i - 1], dtangent_y(vertices[i])) +
                VARIATION_WEIGHT * integrate<dual_value::df>(dtangent_y(vertices[i]), segments[i], vertices[i + 1]);
            *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[i - 1], segments[i - 1], dcurvature_x(vertices[i])) +
                VARIATION_WEIGHT * integrate<dual_value::df>(dcurvature_x(vertices[i]), segments[i], vertices[i + 1]);
            *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[i - 1], segments[i - 1], dcurvature_y(vertices[i])) +
                VARIATION_WEIGHT * integrate<dual_value::df>(dcurvature_y(vertices[i]), segments[i], vertices[i + 1]);
        }

        *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[vertices.size() - 1], segments.back(), dtangent_x(vertices.back()));
        *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[vertices.size() - 1], segments.back(), dtangent_y(vertices.back()));
        *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[vertices.size() - 1], segments.back(), dcurvature_x(vertices.back()));
        *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[vertices.size() - 1], segments.back(), dcurvature_y(vertices.back()));

        // for all segments
        for (size_t i = 0; i < segments.size(); ++i) {
            *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[i], dmagnitude_start(segments[i]), vertices[i + 1]);
            *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[i], dmagnitude_end(segments[i]), vertices[i + 1]);
            *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[i], dalpha_start(segments[i]), vertices[i + 1]);
            *grad++ = VARIATION_WEIGHT * integrate<dual_value::df>(vertices[i], dalpha_end(segments[i]), vertices[i + 1]);
        }
    }

    for (size_t i = 0; i < segments.size(); ++i) {
        variationIntegralSum += integrate<dual_value::f>(vertices[i], segments[i], vertices[i + 1]);
    }

    const auto objective = DISTANCE_WEIGHT * distanceSq + VARIATION_WEIGHT * variationIntegralSum;
    std::cout << "objective=" << objective << '\n';

    return objective;
}

void initialize_interpolation() {
    for (size_t i = 0; i < originalVertices.size(); ++i) {
        vertices.push_back({ originalVertices[i].pos });
        vertices.back().t = originalVertices[i].t;
    }

    for (size_t i = 1; i < vertices.size() - 1; ++i) {
        vertices[i].tangent = (vertices[i - 1].pos - vertices[i].pos) / dot(vertices[i - 1].pos - vertices[i].pos, vertices[i - 1].pos - vertices[i].pos) +
            (vertices[i].pos - vertices[i + 1].pos) / dot(vertices[i].pos - vertices[i + 1].pos, vertices[i].pos - vertices[i + 1].pos);
        vertices[i].tangent = -vertices[i].tangent / sqrt(dot(vertices[i].tangent, vertices[i].tangent));

        vertices[i].curvature = { 1, 1 };
    }

    vertices.front().tangent = vertices[1].tangent;
    vertices.back().tangent = vertices[vertices.size() - 2].tangent;

    for (size_t i = 0; i < vertices.size() - 1; ++i) {
        segments.push_back({ 
            sqrt(dot(vertices[i].pos - vertices[i + 1].pos, vertices[i].pos - vertices[i + 1].pos)),
            sqrt(dot(vertices[i].pos - vertices[i + 1].pos, vertices[i].pos - vertices[i + 1].pos)),
            0,
            0
        });
    }
}

int main() {
    initialize_interpolation();

    const unsigned int taskSize = segments.size() * 4 + vertices.size() * 4 + (vertices.size() - 2) * 2;
    std::cout << "taskSize=" << taskSize << std::endl;
    nlopt::opt optimizer(nlopt::/*LD_SLSQP*/LN_COBYLA, taskSize);
    optimizer.set_min_objective(objective, nullptr);
    optimizer.set_xtol_rel(1e-4);
    optimizer.set_ftol_abs(1e-4);
    optimizer.set_ftol_rel(1e-4);

    std::vector<double> x(taskSize);
    packParameters(&x[0]);

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
            auto&& curvePoint = hermite(start, end, t);
            file << curvePoint.x.x << ", " << curvePoint.y.x << ", " << i + t << ",\n";
        }
    }
    file << "];" << std::endl;

    file << "var originalSpline = [\n";
    for (auto&& vertex : originalVertices) {
        file << vertex.pos.x.x << ", " << vertex.pos.y.x << ", " << vertex.t << ",\n";
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
            file << i + t << ", " << c.x << ", " << i + t << ",\n";
        }
    }
    file << "];" << std::endl;

    gsl_integration_workspace_free(workspace);
}

