#include <nlopt.hpp>
#include <vector>
#include <iostream>
#include <fstream>

struct plane_curve_point {
    double x, y, t;
};

const auto DISTANCE_WEIGHT = 1.0;
const auto VARIATION_WEIGHT = 1.0;

std::vector<plane_curve_point> curve{
    { 10.0, 0.0, 0.0 },
    { 11.6981, 0.0601598, 0.425 },
    { 13.3567, 0.418023, 0.85 },
    { 14.9026, 1.11765, 1.275 },
    { 16.2661, 2.12757, 1.7 },
    { 17.4755, 3.32199, 2.125 },
    { 18.839, 4.3319, 2.55 },
    { 20.3612, 5.08687, 2.975 },
    { 21.8835, 5.84183, 3.4 },
    { 23.247, 6.85175, 3.825 },
    { 24.3668, 8.12652, 4.25 },
    { 25.1927, 9.60881, 4.675 },
    { 25.6873, 1.2319, 5.1 },
    { 25.8285, 2.9228, 5.525 },
    { 25.7293, 4.6197, 5.95 },
    { 25.8106, 6.3169, 6.375 },
    { 25.8919, 8.0141, 6.8 },
    { 25.6733, 9.6968, 7.225 },
    { 25.1047, 1.2955, 7.65 },
    { 24.2118, 2.7383, 8.075 },
    { 23.0346, 3.9604, 8.5 },
    { 21.6261, 4.9066, 8.925 },
    { 20.0498, 5.5346, 9.35 },
    { 18.3765, 5.816, 9.775 },
    { 17.4114, 5.8354, 10.0163705 },
    { 15.7133, 5.7752, 10.4413705 },
    { 14.0305, 5.5352, 10.8663705 },
    { 12.3325, 5.4751, 11.2913705 },
    { 10.6344, 5.4149, 11.7163705 },
    { 8.97573, 5.057, 12.1413705 },
    { 7.42988, 4.3574, 12.5663705 },
    { 6.06635, 3.3475, 12.9913705 },
    { 4.94649, 2.0727, 13.4163705 },
    { 4.12067, 0.5904, 13.8413705 },
    { 3.62605, 8.9673, 14.2663705 },
    { 3.48488, 7.2764, 14.6913705 },
    { 3.7035, 5.5937, 15.1163705 },
    { 4.21584, 3.9737, 15.5413705 },
    { 4.72817, 2.3536, 15.9663705 },
    { 4.94679, 0.6709, 16.3913705 },
    { 4.92561, 8.97128, 16.8163705 },
    { 5.02477, 7.27439, 17.2413705 },
    { 5.00358, 5.57474, 17.6663705 },
};

struct spline_coefficients {
    double x, dx, dx2, dx3, dx4;
};

struct spline_segment_info {
    spline_coefficients x_coeffs, y_coeffs;
};

std::vector<spline_segment_info> segments(curve.size() - 1);

using dvec = std::vector<double>;

double spline(const spline_coefficients& c, const double t) {
    return c.x + c.dx * t + c.dx2 * t * t / 2 + c.dx3 * t * t * t / 6 + c.dx4 * t * t * t * t / 24;
}

double dspline(const spline_coefficients& c, const double t) {
    return c.dx + c.dx2 * t + c.dx3 * t * t / 2 + c.dx4 * t * t * t / 6;
}

double dspline2(const spline_coefficients& c, const double t) {
    return c.dx2 + c.dx3 * t + c.dx4 * t * t / 2;
}

double dspline3(const spline_coefficients& c, const double t) {
    return c.dx3 + c.dx4 * t;
}

double curvature(const spline_segment_info& spline, const double t) {
    const auto dx = dspline(spline.x_coeffs, t);
    const auto dy = dspline(spline.y_coeffs, t);
    const auto dx2 = dspline2(spline.x_coeffs, t);
    const auto dy2 = dspline2(spline.y_coeffs, t);

    return (dx * dy2 - dx2 * dy) / std::pow(dx * dx + dy * dy, 1.5);
}

double dcurvature(const spline_segment_info& spline, const double t) {
    const auto dx = dspline(spline.x_coeffs, t);
    const auto dy = dspline(spline.y_coeffs, t);
    const auto dx2 = dspline2(spline.x_coeffs, t);
    const auto dy2 = dspline2(spline.y_coeffs, t);
    const auto dx3 = dspline3(spline.x_coeffs, t);
    const auto dy3 = dspline3(spline.y_coeffs, t);

    const auto velocitySq = dx * dx + dy * dy;

    return ((dx * dy3 - dx3 * dy) * std::pow(velocitySq, 1.5) - 3 * (dx * dy2 - dx2 * dy) * std::sqrt(velocitySq) * (dx * dx2 + dy * dy2)) / std::pow(velocitySq, 3);
}

double objective(unsigned n, const double*, double*, void*) {
    // Distance objective, only applied for inner points
    auto distanceSq = 0.;

    // @todo start and end points are hard constraints
    for (auto i = 1; i < curve.size() - 1; ++i) {
        const auto distanceX = spline(segments[i].x_coeffs, curve[i].t) - curve[i].x;
        const auto distanceY = spline(segments[i].y_coeffs, curve[i].t) - curve[i].y;
        distanceSq += distanceX * distanceX + distanceY * distanceY;
    }

    // Variation constraint integral
    // @todo replace euler integration with something smarter
    auto variationIntegralSum = 0.;

    for (auto i = 0; i < segments.size(); ++i) {
        const auto t_min = curve[i].t;
        const auto t_max = curve[i + 1].t;
        const auto DT_REFERENCE = 1.0e-3;
        const auto steps = std::ceil((t_max - t_min) / DT_REFERENCE);
        const auto dt = (t_max - t_min) / steps;

        for (auto step = 0; step < steps; ++step) {
            const auto t = t_min + dt * step;
            const auto dc = dcurvature(segments[i], t);

            variationIntegralSum += dc * dc * dt;
        }
    }

    return DISTANCE_WEIGHT * distanceSq + VARIATION_WEIGHT * variationIntegralSum;
}

const auto CONSTRAINTS_PER_KNOT = 8;
const auto CONTINUITY_CONSTRAINTS_NUM = CONSTRAINTS_PER_KNOT * (segments.size() - 1);
void continuity_constraints(unsigned m, double* result, unsigned n, const double* x, double*, void*) {
    auto ptr = result;
    for (auto i = 1; i < segments.size(); ++i, ptr += CONSTRAINTS_PER_KNOT) {
        auto&& prevSeg = segments[i - 1];
        auto&& seg = segments[i];

        // x C3 continuity
        ptr[0] = spline(prevSeg.x_coeffs, curve[i].t) - spline(seg.x_coeffs, curve[i].t);
        ptr[1] = dspline(prevSeg.x_coeffs, curve[i].t) - dspline(seg.x_coeffs, curve[i].t);
        ptr[2] = dspline2(prevSeg.x_coeffs, curve[i].t) - dspline2(seg.x_coeffs, curve[i].t);
        ptr[3] = dspline3(prevSeg.x_coeffs, curve[i].t) - dspline3(seg.x_coeffs, curve[i].t);

        // y C3 continuity
        ptr[4] = spline(prevSeg.y_coeffs, curve[i].t) - spline(seg.y_coeffs, curve[i].t);
        ptr[5] = dspline(prevSeg.y_coeffs, curve[i].t) - dspline(seg.y_coeffs, curve[i].t);
        ptr[6] = dspline2(prevSeg.y_coeffs, curve[i].t) - dspline2(seg.y_coeffs, curve[i].t);
        ptr[7] = dspline3(prevSeg.y_coeffs, curve[i].t) - dspline3(seg.y_coeffs, curve[i].t);
    }
}

const auto START_END_CONSTRAINTS_NUM = 5;
void start_end_constraints(unsigned m, double* result, unsigned n, const double* x, double*, void*) {
    // Start and end must match
    result[0] = spline(segments.front().x_coeffs, curve.front().t) - curve.front().x;
    result[1] = spline(segments.front().y_coeffs, curve.front().t) - curve.front().y;
    result[2] = spline(segments.back().x_coeffs, curve.back().t) - curve.back().x;
    result[3] = spline(segments.back().y_coeffs, curve.back().t) - curve.back().y;

    // @todo
    // second derivative at start
    // second derivative at end
    result[4] = curvature(segments.front(), curve.front().t);
}

const auto CURVATURE_CONSTRAINTS_NUM = segments.size();
const auto CURVATURE_MAX = 0.6;
void curvature_constraints(unsigned m, double* result, unsigned n, const double* x, double*, void*) {
    for (auto i = 0; i < segments.size(); ++i) {
        result[i] = std::abs(curvature(segments[i], curve[i + 1].t)) - CURVATURE_MAX;
    }
}

int main() {
    const auto TASK_SIZE = static_cast<unsigned int>(segments.size() * sizeof(segments[0]) / sizeof(double));
    std::cout << "task size is " << TASK_SIZE << std::endl;
    std::cout << "number of continuity constraints " << CONTINUITY_CONSTRAINTS_NUM << std::endl;
    std::cout << "number of start/end constraints " << START_END_CONSTRAINTS_NUM << std::endl;

    nlopt::opt optimizer{nlopt::LN_COBYLA, TASK_SIZE};
    optimizer.set_xtol_rel(1e-3);
    optimizer.set_ftol_rel(1e-3);
    optimizer.set_ftol_rel(1e-3);
    optimizer.set_min_objective(objective, nullptr);
    //optimizer.add_equality_mconstraint(continuity_constraints, nullptr, dvec(CONTINUITY_CONSTRAINTS_NUM));
    optimizer.add_equality_mconstraint(start_end_constraints, nullptr, dvec(START_END_CONSTRAINTS_NUM));
    //optimizer.add_inequality_mconstraint(curvature_constraints, nullptr, dvec(CURVATURE_CONSTRAINTS_NUM));

    double minFunctional;
    try {
        optimizer.optimize(reinterpret_cast<double*>(segments.data()), TASK_SIZE, minFunctional);
    } catch (const std::exception& exception) {
        std::cout << "e.what(): " << exception.what() << std::endl;
    }

    std::cout << "objective = " << minFunctional << std::endl;

    for (auto i = 1; i < curve.size(); ++i) {
        auto&& seg = segments[i - 1];
        std::cout << "x=" << spline(seg.x_coeffs, curve[i].t) << "\t\ty=" << spline(seg.y_coeffs, curve[i].t) << "\t\tcurvature=" << curvature(seg, curve[i].t) << "\t\tt=" << curve[i].t << '\n';
    }

    std::ofstream file("spline.js");

    file << "var splineTMin = " << curve.front().t << ", splineTMax = " << curve.back().t << std::endl;
    // print original spline
    file << "var originalSpline = [\n";
    for (auto i = 0; i < curve.size(); ++i) {
        file << '\t' << curve[i].x << ", " << curve[i].y << ", " << curve[i].t << ",\n";
    }
    file << ']' << std::endl;

    // print x(t)
    file << "var xSegments = [\n";
    for (auto i = 0; i < segments.size(); ++i) {
        file << '\t' << curve[i].t << ", " << curve[i + 1].t << ", " << segments[i].x_coeffs.x << ", " << segments[i].x_coeffs.dx << ", " << segments[i].x_coeffs.dx2 << ", " << segments[i].x_coeffs.dx3 << ", " << segments[i].x_coeffs.dx4 << ",\n";
    }
    file << ']' << std::endl;

    // print y(t)
    file << "var ySegments = [\n";
    for (auto i = 0; i < segments.size(); ++i) {
        file << '\t' << curve[i].t << ", " << curve[i + 1].t << ", " << segments[i].y_coeffs.x << ", " << segments[i].y_coeffs.dx << ", " << segments[i].y_coeffs.dx2 << ", " << segments[i].y_coeffs.dx3 << ", " << segments[i].y_coeffs.dx4 << ",\n";
    }
    file << ']' << std::endl;
}
