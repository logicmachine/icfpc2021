/////////////////////////////////////////////////////////////////////
// Brain Wall Solver  by foota, 2021/07/09
// g++ brain_wall_solver.cpp -O2 -Wall -std=c++17

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "json.hpp"

class BrainWallSolver {
private:
    int num_iters;
    int num_vertices_of_hole;
    int num_vertices;
    int num_edges;
    std::vector<std::vector<int>> hole;
    std::vector<std::vector<int>> vertices;
    std::vector<std::vector<double>> forces;
    std::vector<std::vector<int>> edges;
    std::vector<double> ene;
    double eps;

public:
    BrainWallSolver(int n) : num_iters(n) {}
    ~BrainWallSolver(){};
    void read(const std::string);
    void calc();
    void move();
    void print();
    void run();
    bool crossing(const std::vector<int>& a, const std::vector<int>& b, const std::vector<int>& c, const std::vector<int>& d);
    bool crossing_hole(const std::vector<int>& a, const std::vector<int>& b);
    bool point_in_area(const std::vector<int>& p);
    std::vector<double> nearest_point_on_line(const std::vector<int>& p1, const std::vector<int>& p2, const std::vector<int>& point);
    std::vector<int> nearest_point_on_hole(const std::vector<int>& point);
    //double dist_point_to_line(const std::vector<int>& p1, const std::vector<int>& p2, const std::vector<int> point);

    double energy() {
        double sum = 0.0;
        for (int i = 0; i < num_vertices; i++) sum += ene[i];
        return sum;
    }
    template <typename U, typename T>
    double dist(const std::vector<U>& p, const std::vector<T>& q) {
        double dx = p[0] - q[0];
        double dy = p[1] - q[1];
        return dx * dx + dy * dy;
    }
};

void BrainWallSolver::read(const std::string fname) {
    std::fstream fs(fname, std::ios_base::in);

    if (!fs.is_open()) {
        std::cerr << "File open error: " << fname << std::endl;
        exit(1);
    }

    fs >> num_vertices_of_hole;
    int x, y;
    for (int i = 0; i < num_vertices_of_hole; i++) {
        fs >> x >> y;
        hole.push_back(std::vector<int>({x, y}));
    }
    fs >> num_vertices >> num_edges;
    for (int i = 0; i < num_vertices; i++) {
        fs >> x >> y;
        vertices.push_back(std::vector<int>({x, y}));
    }
    for (int i = 0; i < num_edges; i++) {
        fs >> x >> y;
        edges.push_back(std::vector<int>({x, y}));
    }
    fs >> eps;
    eps /= 1000000.0;

    fs.close();

    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 2; j++) forces.push_back(std::vector<double>({0.0, 0.0}));
        ene.push_back(0.0);
    }
}

bool BrainWallSolver::crossing(const std::vector<int>& a, const std::vector<int>& b, const std::vector<int>& c, const std::vector<int>& d) {
    double s, t;

    s = (a[0] - b[0]) * (c[1] - a[1]) - (a[1] - b[1]) * (c[0] - a[0]);
    t = (a[0] - b[0]) * (d[1] - a[1]) - (a[1] - b[1]) * (d[0] - a[0]);
    if (s * t >= 0.0) return false;

    s = (c[0] - d[0]) * (a[1] - c[1]) - (c[1] - d[1]) * (a[0] - c[0]);
    t = (c[0] - d[0]) * (b[1] - c[1]) - (c[1] - d[1]) * (b[0] - c[0]);
    if (s * t >= 0.0) return false;

    return true;
}

bool BrainWallSolver::crossing_hole(const std::vector<int>& a, const std::vector<int>& b) {
    for (int i = 0; i < num_vertices_of_hole - 1; i++) {
        if (crossing(hole[i], hole[i + 1], a, b)) return true;
    }
    if (crossing(hole[num_vertices_of_hole - 1], hole[0], a, b)) return true;

    return false;
}

#if 0
bool BrainWallSolver::point_in_area(const std::vector<int>& p) {
    int cn = 0;
    
    for (int i = 0; i < num_vertices_of_hole - 1; i++) {
        if (((hole[i][1] <= p[1]) && (hole[i+1][1] > p[1])) || ((hole[i][1] > p[1]) && (vertices[i+1][1] <= p[1]))) {
            //double vt = (p[1] - hole[i][1]) / (hole[i+1][1] - hole[i][1]);
            double vt = (hole[i+1][1] - hole[i][1]) != 0.0 ? (p[1] - hole[i][1]) / (hole[i+1][1] - hole[i][1]) : 0.0;
            if(p[0] < (hole[i][0] + (vt * (hole[i+1][0] - hole[i][0])))) {
                cn++;
            }
        }
    }
    
    if ((cn & 1) == 1) return true;
    return false;
}
#endif

bool BrainWallSolver::point_in_area(const std::vector<int>& p) {
    int cnt = 0;

    std::vector<int> p0 = hole[0];
    bool f0x = (p[0] <= p0[0]);
    bool f0y = (p[1] <= p0[1]);

    for (int i = 1; i < num_vertices_of_hole + 1; i++) {
        std::vector<int> p1(hole[i % num_vertices_of_hole]);
        bool f1x = (p[0] <= p1[0]);
        bool f1y = (p[1] <= p1[1]);
        if (f0y != f1y) {
            if (f0x == f1x) {
                if (f0x) cnt += (f0y ? -1 : 1);
            } else {
                if (p[0] <= (p0[0] + (p1[0] - p0[0]) * (p[1] - p0[1]) / (p1[1] - p0[1]))) cnt += (f0y ? -1 : 1);
            }
        }
        p0 = p1;
        f0x = f1x;
        f0y = f1y;
    }

    if (cnt == 0) {
        for (int i = 0; i < num_vertices_of_hole; i++) {
            if (p == hole[i]) return true;
        }
    }

    return cnt != 0;
}

#if 0
double BrainWallSolver::dist_point_to_line(const std::vector<int>& p1, const std::vector<int>& p2, const std::vector<int> point) {
    double x0 = point[0];
    double y0 = point[1];
    double x1 = p1[0];
    double y1 = p1[1];
    double x2 = p2[0];
    double y2 = p2[1];

    double a = x2 - x1;
    double b = y2 - y1;
    double a2 = a * a;
    double b2 = b * b;
    double r2 = a2 + b2;
    double tt = -(a * (x1 - x0) + b * (y1 - y0));

    if (tt < 0) return std::sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
    else if (tt > r2) return std::sqrt((x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0));

    double f1 = a * (y1 - y0) - b * (x1 - x0);
    return std::sqrt((f1 * f1) / r2);
}
#endif

std::vector<double> BrainWallSolver::nearest_point_on_line(const std::vector<int>& p1, const std::vector<int>& p2, const std::vector<int>& point) {
    double x0 = point[0];
    double y0 = point[1];
    double x1 = p1[0];
    double y1 = p1[1];
    double x2 = p2[0];
    double y2 = p2[1];

    double ax = x2 - x1;
    double ay = y2 - y1;
    double bx = x0 - x1;
    double by = y0 - y1;

    double r = (ax * bx + ay * by) / (ax * ax + ay * ay);

    if (r <= 0)
        return {(x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0), x1, y1};
    else if (r >= 1)
        return {(x2 - x0) * (x2 - x0) + (y2 - y0) * (y2 - y0), x2, y2};
    else {
        double x3 = p1[0] + r * ax;
        double y3 = p1[1] + r * ay;
        return {(x3 - x0) * (x3 - x0) + (y3 - y0) * (y3 - y0), x3, y3};
    }
}

std::vector<int> BrainWallSolver::nearest_point_on_hole(const std::vector<int>& point) {
    std::vector<int> ans(hole[0]);
    double min_d = dist(hole[0], point);

    for (int i = 1; i < num_vertices_of_hole; i++) {
        double d = dist(hole[i], point);
        if (d < min_d) {
            min_d = d;
            ans = hole[i];
        }
    }

    return ans;
}

// Calculate energies and forces
// Total energy = vdW + Coulomb
// vdW
//  U = eps * [(Rij / r[i])^12 - 2 * (Rij / r[i])^6]
//  F = -12 * eps / Rij * [(Rij / r[i])^13 - (Rij / r[i])^7] * r_xyz / r[i]
// Coulomb
//  U = SUM_i>j qiqj / r[i]
//  F = SUM_j qiqj / r[i]^3 * r_xyz
void BrainWallSolver::calc() {
    //std::vector<std::vector<double>> nearest_points;
    std::vector<std::vector<int>> nearest_points;
    std::vector<bool> is_point_in_area;
    for (int i = 0; i < num_vertices; i++) {
        is_point_in_area.push_back(point_in_area(vertices[i]));
#if 0
        std::vector<std::vector<double>> nearest_points_;
        for (int j = 0; j < num_vertices_of_hole - 1; j++) {
            nearest_points_.push_back(nearest_point_on_line(hole[j], hole[j+1], vertices[i]));
        }
        nearest_points_.push_back(nearest_point_on_line(hole[num_vertices_of_hole-1], hole[0], vertices[i]));
        //std::sort(nearest_points_.begin(), nearest_points_.end(), [](auto const& a, auto const& b){ return a[0] < b[0]; });
        std::sort(nearest_points_.begin(), nearest_points_.end());
        nearest_points.push_back(nearest_points_.front());
        //for (auto x : nearest_points) std::cerr << x[0] << std::endl; ///// debug
#endif
        nearest_points.push_back(nearest_point_on_hole(vertices[i]));
    }

    /// debug
    std::cerr << "point in area / nearest point:" << std::endl;
    for (int i = 0; i < num_vertices; i++) {
        //std::cerr << is_point_in_area[i] << ", (" << nearest_points[i][1] << ", " << nearest_points[i][2] << ") " << nearest_points[i][0] << std::endl;
        std::cerr << is_point_in_area[i] << ", (" << nearest_points[i][0] << ", " << nearest_points[i][1] << ")" << std::endl;
    }

    std::vector<bool> is_crossing_edges;
    for (int i = 0; i < num_edges; i++) {
        int a = edges[i][0];
        int b = edges[i][1];
        if (is_point_in_area[a] && is_point_in_area[b]) {
            is_crossing_edges.push_back(crossing_hole(vertices[a], vertices[b]));
            //is_crossing_edges.push_back(false);
        } else
            is_crossing_edges.push_back(true);
    }

    /// debug
    std::cerr << "crossing edges:" << std::endl;
    for (int i = 0; i < num_edges; i++) {
        std::cerr << "(" << edges[i][0] << ", " << edges[i][1] << ") : " << is_crossing_edges[i] << std::endl;
    }

#if 0
    const double eps = 0.2f;
    const double Rij = 2.5f;
    double r, r_xyz, df, r12, r6, q;

    for (int i = 0; i < num_vertices; i++) {
        ene[i] = 0.0f;
        for (int j = 0; j < num_vertices; j++) {
            if (i == j) continue;
            r = dist(i, j);
            q = chg[i] * chg[j];
            r12 = pow(Rij / r, 12);
            r6 = pow(Rij / r, 6);
            ene[i] += eps * (r12 - 2.0f * r6) + q / r;
            for (int k = 0; k < 3; k++) {
                r_xyz = dist(i, j, k);
//                df = -(12 * eps / Rij) * (r12 * Rij / r - r6 * Rij / r) * r_xyz / r;
                df = (12 * eps / Rij) * (r12 * Rij / r - r6 * Rij / r) * r_xyz / r;
                df += q / pow(r, 3) * r_xyz;
                f[i * 3 + k] += df;
            }
        }
    }
#endif
}

void BrainWallSolver::move() {
    double pos[2];
    double rate = 0.01;
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 2; j++) {
            pos[j] = vertices[i][j] + forces[i][j] * rate;
            /// TODO: verification
            vertices[i][j] = static_cast<int>(pos[j] + 0.5);
            forces[i][j] = 0.0;
        }
#if 0
    double r, dr;
        r = 0.0;
        for (int j = 0; j < 2; j++) r += vertices[i][j] * vertices[i][j];
        dr = r - cap_range * cap_range;
        if (dr > 0.0) {
            for (int j = 0; j < 2; j++)
                forces[i][j] -= crd[i][j] / cap_range * dr * 0.01;
        }
#endif
    }
}

void BrainWallSolver::print() {
    std::cout << std::endl;
    std::cout << num_vertices_of_hole << std::endl;
    std::cout << num_vertices << std::endl;
    std::cout << "num_vertices: " << num_vertices << std::endl;
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 2; j++) std::cout << " " << std::fixed << std::setw(10) << std::setprecision(6) << vertices[i][j];
        std::cout << std::endl;
    }
#if 0
    std::cout << "num_edges: " << num_edges << std::endl;
    for (int i = 0; i < num_edges; i++) {
        for (int j = 0; j < 2; j++) std::cout << " " << std::fixed << std::setw(10) << std::setprecision(6) << edges[i][j];
        std::cout << std::endl;
    }
#endif
    std::cout << "eps: " << eps << std::endl;
}

void BrainWallSolver::run() {
    for (int i = 0; i < num_iters; i++) {
        calc();
        move();
        std::cout << "Energy (" << std::setw(7) << i + 1 << "): ";
        std::cout << std::fixed << std::setw(15) << std::setprecision(5) << energy() << std::endl;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input_file [num_iters=1]" << std::endl;
        return 1;
    }

    std::string fname(argv[1]);
    int num_iters = argc > 2 ? std::stol(argv[2]) : 1;

    BrainWallSolver solver(num_iters);

    solver.read(fname);
    solver.print();
    solver.run();
    solver.print();

    return 0;
}
