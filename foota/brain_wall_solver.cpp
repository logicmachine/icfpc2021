/////////////////////////////////////////////////////////////////////
// Brain Wall Solver  by foota, 2021/07/09
// g++ brain_wall_solver.cpp -O2 -Wall -std=c++17

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <random>
#include "json.hpp"

std::random_device rd;
std::mt19937 mt(rd());

class BrainWallSolver {
private:
    int num_iters;
    int num_vertices_of_hole;
    int num_vertices;
    int num_edges;
    std::vector<std::vector<int>> hole;
    std::vector<int> hole_params; // (min_x, max_x, min_y, max_y, center_x, center_y, centroid_x, centroid_y)
    std::vector<std::vector<int>> vertices;
    std::vector<int> figure_params; // (min_x, max_x, min_y, max_y, center_x, center_y, centroid_x, centroid_y)
    std::vector<std::vector<double>> forces;
    std::vector<std::vector<int>> edges;
    std::vector<int> figure_distances;
    std::vector<double> ene;
    double eps;

public:
    BrainWallSolver(int n) : num_iters(n) {}
    ~BrainWallSolver(){};
    //void read(const std::string);
    void read();
    double dislike();
    double energy();
    void first_step();
    void second_step();
    void calc();
    void move();
    bool validate_distance(int e);
    int validation();
    void print();
    void run();
    bool crossing(const std::vector<int>& a, const std::vector<int>& b, const std::vector<int>& c, const std::vector<int>& d);
    bool crossing_hole(const std::vector<int>& a, const std::vector<int>& b);
    bool point_in_area(const std::vector<int>& p);
    std::vector<double> nearest_point_on_line(const std::vector<int>& p1, const std::vector<int>& p2, const std::vector<int>& point);
    //std::vector<int> nearest_point_on_hole(const std::vector<int>& point);
    double nearest_point_on_hole(const std::vector<int>& point, std::vector<int>& ans);
    //double dist_point_to_line(const std::vector<int>& p1, const std::vector<int>& p2, const std::vector<int> point);
    void output();
    void output_json(std::string fname);

    template <typename U, typename T>
    double dist(const std::vector<U>& p, const std::vector<T>& q) {
        double dx = p[0] - q[0];
        double dy = p[1] - q[1];
        return dx * dx + dy * dy;
    }
};

void BrainWallSolver::read() {
    std::cin >> num_vertices_of_hole;
    int x, y;
    for (int i = 0; i < num_vertices_of_hole; i++) {
        std::cin >> x >> y;
        hole.push_back(std::vector<int>({x, y}));
    }
    std::cin >> num_vertices >> num_edges;
    for (int i = 0; i < num_vertices; i++) {
        std::cin >> x >> y;
        vertices.push_back(std::vector<int>({x, y}));
    }
    for (int i = 0; i < num_edges; i++) {
        std::cin >> x >> y;
        edges.push_back(std::vector<int>({x, y}));
        figure_distances.push_back(static_cast<int>(std::round(dist(vertices[x], vertices[y]))));
    }
    std::cin >> eps;
    eps /= 1000000.0;
    if (eps == 0.0) eps = 1e-6;

    int min_x, max_x, min_y, max_y;
    double center_x, center_y, centroid_x, centroid_y;
    min_x = max_x = hole[0][0];
    min_y = max_y = hole[0][1];
    centroid_x = hole[0][0];
    centroid_y = hole[0][1];
    for (int i = 1; i < num_vertices_of_hole; i++) {
        min_x = std::min(hole[i][0], min_x);
        max_x = std::max(hole[i][0], max_x);
        min_y = std::min(hole[i][1], min_y);
        max_y = std::max(hole[i][1], max_y);
        centroid_x += hole[i][0];
        centroid_y += hole[i][1];
    }
    centroid_x /= num_vertices_of_hole;
    centroid_y /= num_vertices_of_hole;
    center_x = (min_x + max_x) / 2.0;
    center_y = (min_y + max_y) / 2.0;
    hole_params = {min_x, max_x, min_y, max_y,
                   static_cast<int>(std::round(center_x)),
                   static_cast<int>(std::round(center_y)),
                   static_cast<int>(std::round(centroid_x)),
                   static_cast<int>(std::round(centroid_y))};

    min_x = max_x = vertices[0][0];
    min_y = max_y = vertices[0][1];
    centroid_x = vertices[0][0];
    centroid_y = vertices[0][1];
    for (int i = 1; i < num_vertices; i++) {
        min_x = std::min(vertices[i][0], min_x);
        max_x = std::max(vertices[i][0], max_x);
        min_y = std::min(vertices[i][1], min_y);
        max_y = std::max(vertices[i][1], max_y);
        centroid_x += vertices[i][0];
        centroid_y += vertices[i][1];
    }
    centroid_x /= num_vertices;
    centroid_y /= num_vertices;
    center_x = (min_x + max_x) / 2.0;
    center_y = (min_y + max_y) / 2.0;
    figure_params = {min_x, max_x, min_y, max_y,
                     static_cast<int>(std::round(center_x)),
                     static_cast<int>(std::round(center_y)),
                     static_cast<int>(std::round(centroid_x)),
                     static_cast<int>(std::round(centroid_y))};

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

double BrainWallSolver::nearest_point_on_hole(const std::vector<int>& point, std::vector<int>& ans) {
    //std::vector<int> ans(hole[0]);
    ans = hole[0];
    double min_d = dist(hole[0], point);

    for (int i = 1; i < num_vertices_of_hole; i++) {
        double d = dist(hole[i], point);
        if (d < min_d) {
            min_d = d;
            ans = hole[i];
        }
    }

    return min_d;
}

double BrainWallSolver::dislike() {
    double energy = 0.0;
#if 0
    std::vector<int> ans;
    for (int i = 0; i < num_vertices; i++) {
        energy += nearest_point_on_hole(vertices[i], ans);
    }
#endif
    for (int i = 1; i < num_vertices_of_hole; i++) {
        //std::vector<int> ans(vertices[0]);
        double min_d = dist(vertices[0], hole[i]);
        for (int j = 0; j < num_vertices; j++) {
            double d = dist(vertices[j], hole[i]);
            if (d < min_d) {
                min_d = d;
                //ans = vertices[j];
            }
        }
        energy += min_d;
    }
    return energy;
}

 double BrainWallSolver::energy() {
     double sum = 0.0;
     //for (int i = 0; i < num_vertices; i++) sum += ene[i];
     sum += dislike();
     for (int e = 0; e < num_edges; e++) {
         int i = edges[e][0];
         int j = edges[e][1];
         double d = dist(vertices[i], vertices[j]);
         double d0 = figure_distances[e];
         double d1 = std::abs(d - d0);
         sum += d1 * 1.0;
         //if (std::abs(d / d0 - 1.0) > eps) sum += d1 * 10.0;
     }

#if 0
    int cnt_in_area = 0;
    std::vector<bool> is_point_in_area;
    for (int i = 0; i < num_vertices; i++) {
        bool is_in = point_in_area(vertices[i]);
        is_point_in_area.push_back(is_in);
        if (!is_in) cnt_in_area++;
    }
    sum += cnt_in_area * 10.0;

    int cnt_crossing = 0;
    std::vector<bool> is_crossing_edges;
    for (int i = 0; i < num_edges; i++) {
        int a = edges[i][0];
        int b = edges[i][1];
        int x = (vertices[a][0] + vertices[b][0]) / 2;
        int y = (vertices[a][1] + vertices[b][1]) / 2;
        if (is_point_in_area[a] && is_point_in_area[b] && point_in_area(std::vector<int>({x, y}))) {
        //if (is_point_in_area[a] && is_point_in_area[b]) {
            bool is_crossing = crossing_hole(vertices[a], vertices[b]);
            is_crossing_edges.push_back(is_crossing);
            if (is_crossing) cnt_crossing++;
        } else {
            is_crossing_edges.push_back(true);
            cnt_crossing++;
        }
    }
    sum += cnt_crossing * 10.0;
#endif

    return sum;
 }

void BrainWallSolver::first_step() {
    int offset_x = hole_params[4] - figure_params[4]; // center_x
    int offset_y = hole_params[5] - figure_params[5]; // center_y
    for (int i = 0; i < num_vertices; i++) {
        vertices[i][0] += offset_x;
        vertices[i][1] += offset_y;
    }
}

void BrainWallSolver::second_step() {
    std::vector<std::vector<int>> nearest_points;
    //std::vector<bool> is_point_in_area;
#if 0
    for (int i = 0; i < num_vertices; i++) {
        if (!point_in_area(vertices[i])) {
            std::vector<int> ans;
            nearest_point_on_hole(vertices[i], ans);
            vertices[i] = ans;
        }
    }
#endif
    for (int i = 0; i < num_vertices; i++) {
        bool in_area = point_in_area(vertices[i]);
        //is_point_in_area.push_back(in_area);
        if (in_area) continue;
        //is_point_in_area.push_back(point_in_area(vertices[i]));
        std::vector<std::vector<double>> nearest_points_;
        for (int j = 0; j < num_vertices_of_hole - 1; j++) {
            nearest_points_.push_back(nearest_point_on_line(hole[j], hole[j+1], vertices[i]));
        }
        nearest_points_.push_back(nearest_point_on_line(hole[num_vertices_of_hole-1], hole[0], vertices[i]));
        //std::sort(nearest_points_.begin(), nearest_points_.end(), [](auto const& a, auto const& b){ return a[0] < b[0]; });
        std::sort(nearest_points_.begin(), nearest_points_.end());
        auto np = nearest_points_.front();
        int x = static_cast<int>(std::round(np[1]));
        int y = static_cast<int>(std::round(np[2]));
        //nearest_points.push_back({x, y});
        //vertices[i] = nearest_points.front();
        vertices[i] = std::vector<int>({x, y});
    }

#if 0
    //for (;;) {
    //bool is_check = true;
    for (int e = 0; e < num_edges; e++) {
        int i = edges[e][0];
        int j = edges[e][1];
#if 0
        if (vertices[i] == vertices[j]) {
            for (const auto& pos : hole) {
                if (vertices[j] != pos) vertices[j] = pos;
            }
        }
#endif
#if 0
        if (vertices[i] == vertices[j]) {
            //is_check = false;
            bool is_found = false;
            for (int x = std::max(vertices[i][0] - 5, hole_params[0]); x < std::min(vertices[i][0] + 5, hole_params[1]); x++) {
                for (int y = std::max(vertices[i][1] - 5, hole_params[2]); y < std::min(vertices[i][1] + 5, hole_params[3]); y++) {
                    auto pos = std::vector<int>({x, y});
                    if (point_in_area(pos)) {
                        vertices[j] = pos;
                        is_found = true;
                        break;
                    }
                }
                if (is_found) break;
            }
        }
#endif
    }
    //if (is_check) break;
    //}
#endif
#if 0
    std::vector<std::vector<int>> nearest_points;
    std::vector<bool> is_point_in_area;
    for (int i = 0; i < num_vertices; i++) {
        is_point_in_area.push_back(point_in_area(vertices[i]));
        std::vector<int> ans;
        nearest_point_on_hole(vertices[i], ans);
        nearest_points.push_back(ans);
    }

    std::vector<bool> is_crossing_edges;
    for (int i = 0; i < num_edges; i++) {
        int a = edges[i][0];
        int b = edges[i][1];
        if (is_point_in_area[a] && is_point_in_area[b]) {
            is_crossing_edges.push_back(crossing_hole(vertices[a], vertices[b]));
        } else
            is_crossing_edges.push_back(true);
    }
#endif
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
        //nearest_points.push_back(nearest_point_on_hole(vertices[i]));
        std::vector<int> ans;
        nearest_point_on_hole(vertices[i], ans);
        nearest_points.push_back(ans);
    }

#if 0
    /// debug
    std::cerr << "point in area / nearest point:" << std::endl;
    for (int i = 0; i < num_vertices; i++) {
        //std::cerr << is_point_in_area[i] << ", (" << nearest_points[i][1] << ", " << nearest_points[i][2] << ") " << nearest_points[i][0] << std::endl;
        std::cerr << is_point_in_area[i] << ", (" << nearest_points[i][0] << ", " << nearest_points[i][1] << ")" << std::endl;
    }
#endif
    std::vector<bool> is_crossing_edges;
    for (int i = 0; i < num_edges; i++) {
        int a = edges[i][0];
        int b = edges[i][1];
        int x = (vertices[a][0] + vertices[b][0]) / 2;
        int y = (vertices[a][1] + vertices[b][1]) / 2;
        if (is_point_in_area[a] && is_point_in_area[b] && point_in_area(std::vector<int>({x, y}))) {
            is_crossing_edges.push_back(crossing_hole(vertices[a], vertices[b]));
            //is_crossing_edges.push_back(false);
        } else
            is_crossing_edges.push_back(true);
    }

#if 0
    /// debug
    std::cerr << "crossing edges:" << std::endl;
    for (int i = 0; i < num_edges; i++) {
        std::cerr << "(" << edges[i][0] << ", " << edges[i][1] << ") : " << is_crossing_edges[i] << std::endl;
    }
#endif
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
                //df = -(12 * eps / Rij) * (r12 * Rij / r - r6 * Rij / r) * r_xyz / r;
                df = (12 * eps / Rij) * (r12 * Rij / r - r6 * Rij / r) * r_xyz / r;
                df += q / pow(r, 3) * r_xyz;
                f[i * 3 + k] += df;
            }
        }
    }
#endif
}

void BrainWallSolver::move() {
    std::uniform_int_distribution<> randint(0, num_vertices - 1);
    std::uniform_real_distribution<> rand(0.0, 1.0);
    for (int e = 0; e < num_edges; e++) {
        int i = edges[e][0];
        int j = edges[e][1];
        double d0 = figure_distances[e];
#if 0
        if (rand(mt) < 0.1) {
            int k = randint(mt);
            std::swap(vertices[j][0], vertices[k][0]);
            std::swap(vertices[j][1], vertices[k][1]);
            //std::swap(vertices[j], vertices[k]);
        }
#endif
        std::vector<std::pair<double, std::vector<int>>> candidates;
        double d = dist(vertices[i], vertices[j]);
        //double add_ene = 0.0;
        if (std::abs(d / d0 - 1.0) > eps) {
            //bool is_found = false;
            for (int x = hole_params[0]; x <= hole_params[1]; x++) {
            //for (int x = vertices[j][0] - 10; x <= vertices[j][0] + 10; x++) {
                for (int y = hole_params[2]; y <= hole_params[3]; y++) {
                //for (int y = vertices[j][1] - 10; y <= vertices[j][1] + 10; y++) {
                    auto pos = std::vector<int>({x, y});
                    auto pos2 = std::vector<int>({(vertices[i][0] + x) / 2, (vertices[i][1] + y) / 2});
                    d = dist(vertices[i], pos);
                    
                    if (std::abs(d / d0 - 1.0) > eps || !point_in_area(pos) || crossing_hole(vertices[i], pos) || !point_in_area(pos2)) continue;
                    
                    /*
                    if (std::abs(d / d0 - 1.0) > eps) add_ene += std::abs(d - d0) * 10.0;
                    if (!point_in_area(pos)) add_ene += 100.0;
                    if (crossing_hole(vertices[i], pos) || !point_in_area(pos2)) add_ene += 100.0;
                    */
                    auto prev = vertices[j];
                    vertices[j] = pos;
                    //candidates.push_back(std::make_pair(energy() + add_ene, pos));
                    candidates.push_back(std::make_pair(energy(), pos));
                    vertices[j] = prev;
                    //is_found = true;
                }
                //if (is_found) break;
            }
            if (!candidates.empty()) {
                std::uniform_int_distribution<> randcd(0, candidates.size() - 1);
                std::sort(candidates.begin(), candidates.end());
                //std::reverse(candidates.begin(), candidates.end());
                if (rand(mt) < 0.1) {
                    vertices[j] = candidates[randcd(mt)].second;
                } else {
                    vertices[j] = candidates.front().second;
                }
                /*
                if (rand(mt) < 0.1) {
                    vertices[j][0] += randint(mt) - num_vertices / 2;
                    vertices[j][1] += randint(mt) - num_vertices / 2;
                }
                */
            }
            /*
            if (!is_found) {
                std::swap(vertices[i][0], vertices[j][0]);
                std::swap(vertices[i][1], vertices[j][1]);
            }
            */
            if (candidates.empty()) {
                std::cerr << "not found" << std::endl;
                int k = randint(mt);
                int l = randint(mt);
                //std::swap(vertices[l][0], vertices[k][0]);
                //std::swap(vertices[l][1], vertices[k][1]);
                std::swap(vertices[l], vertices[k]);
            }
            //if (vertices[i] == vertices[j]) vertices[j][0] += 1;
            /*
            if (!is_found) {
                vertices[i][0] += 1;
                vertices[i][1] += 1;
            }
            */
        }
    }
#if 0
    double pos[2];
    double rate = 0.01;
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 2; j++) {
            pos[j] = vertices[i][j] + forces[i][j] * rate;
            /// TODO: validation
            vertices[i][j] = static_cast<int>(std::round(pos[j]));
            forces[i][j] = 0.0;
        }
    double r, dr;
        r = 0.0;
        for (int j = 0; j < 2; j++) r += vertices[i][j] * vertices[i][j];
        dr = r - cap_range * cap_range;
        if (dr > 0.0) {
            for (int j = 0; j < 2; j++)
                forces[i][j] -= crd[i][j] / cap_range * dr * 0.01;
        }
    }
#endif
}

bool BrainWallSolver::validate_distance(int e) {
    int i = edges[e][0];
    int j = edges[e][1];
    double d0 = figure_distances[e];
    double d = dist(vertices[i], vertices[j]);
    return std::abs(d / d0 - 1.0) <= eps;
}

int BrainWallSolver::validation() {
    int cnt_in_area = 0;
    int cnt_crossing = 0;
    int cnt_distance = 0;

    std::vector<bool> is_point_in_area;
    for (int i = 0; i < num_vertices; i++) {
        bool is_in = point_in_area(vertices[i]);
        is_point_in_area.push_back(is_in);
        if (!is_in) cnt_in_area++;
    }

    std::vector<bool> is_crossing_edges;
    for (int i = 0; i < num_edges; i++) {
        int a = edges[i][0];
        int b = edges[i][1];
        int x = (vertices[a][0] + vertices[b][0]) / 2;
        int y = (vertices[a][1] + vertices[b][1]) / 2;
        if (is_point_in_area[a] && is_point_in_area[b] && point_in_area(std::vector<int>({x, y}))) {
            bool is_crossing = crossing_hole(vertices[a], vertices[b]);
            is_crossing_edges.push_back(is_crossing);
            if (is_crossing) cnt_crossing++;
        } else {
            is_crossing_edges.push_back(true);
            cnt_crossing++;
        }

        if (!validate_distance(i)) cnt_distance++;
    }

    std::cerr << "area, crossing, distance: "<< cnt_in_area << ", " << cnt_crossing << ", " << cnt_distance << std::endl;

    return cnt_in_area + cnt_crossing + cnt_distance;
}

void BrainWallSolver::print() {
    std::cerr << std::endl;
    std::cerr << num_vertices_of_hole << std::endl;
    std::cerr << num_vertices << std::endl;
    std::cerr << "num_vertices: " << num_vertices << std::endl;
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 2; j++) std::cerr << " " << std::fixed << std::setw(10) << std::setprecision(6) << vertices[i][j];
        std::cerr << std::endl;
    }
#if 0
    std::cerr << "num_edges: " << num_edges << std::endl;
    for (int i = 0; i < num_edges; i++) {
        for (int j = 0; j < 2; j++) std::cerr << " " << std::fixed << std::setw(10) << std::setprecision(6) << edges[i][j];
        std::cerr << std::endl;
    }
#endif
    std::cerr << "eps: " << eps << std::endl;
    std::cerr << "Invalid count: " << validation() << std::endl;
}

void BrainWallSolver::output() {
    //for (const auto& vertex : vertices) {
    for (int i = 0; i < num_vertices; i++) {
        std::cout << vertices[i][0] << " " << vertices[i][1] << std::endl;
    }
}

void BrainWallSolver::output_json(std::string fname) {
    std::fstream fs(fname, std::ios_base::out);
    fs << "{\"vertices\": [";
    for (int i = 0; i < num_vertices; i++) {
        fs << "[" << vertices[i][0] << "," << vertices[i][1] << "]" << (i < num_vertices - 1 ? "," : "");
    }
    fs << "]}" << std::endl;
    fs.close();
}

void BrainWallSolver::run() {
    std::cerr << "Energy (-): " << std::fixed << std::setw(15) << std::setprecision(5) << energy() << std::endl;
    first_step();
    std::cerr << "Energy (0a): " << std::fixed << std::setw(15) << std::setprecision(5) << energy() << std::endl;
    second_step();
    std::cerr << "Energy (0b): " << std::fixed << std::setw(15) << std::setprecision(5) << energy() << std::endl;
    for (int i = 0; i < num_iters; i++) {
        calc();
        move();
        std::cerr << "Energy (" << std::setw(7) << i + 1 << "): ";
        //std::cerr << std::fixed << std::setw(15) << std::setprecision(5) << energy() << std::endl;
        std::cerr << std::fixed << std::setw(15) << std::setprecision(5) << energy() << std::endl;
        if (validation() == 0) break;
    }
}

int main(int argc, char** argv) {
#if 0
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [num_iters=1]" << std::endl;
        return 1;
    }
#endif

    int num_iters = argc > 1 ? std::stol(argv[1]) : 10000;
    //std::string fname(argc > 2 ? argv[2]: "");
    //std::string fname("");

    BrainWallSolver solver(num_iters);

    solver.read();
    solver.print();
    solver.run();
    solver.print();
    if (solver.validation() == 0) {
        solver.output();
    }
    solver.output_json("test.json");

    return 0;
}
