/////////////////////////////////////////////////////////////////////
// Brain Wall Solver  by foota, 2021/07/09
// g++ brain_wall_solver.cpp -O2 -Wall -std=c++20

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <random>
#include <chrono>
#include <ranges>
//#include "json.hpp"

std::random_device rd;
std::mt19937 mt(rd());
std::chrono::system_clock::time_point start_time;
constexpr double DEFAULT_TIME_LIMIT = 300.0;
//constexpr double DEFAULT_TIME_LIMIT = 60.0;
//constexpr double DEFAULT_TIME_LIMIT = 20.0;

std::map<std::string, int> BONUS_TYPES({{"GLOBALIST", 0}, {"BREAK_A_LEG", 1}, {"WALLHACK", 2}});

class BrainWallSolver {
private:
    int num_vertices_of_hole;
    int num_vertices;
    int num_edges;
    int num_bonus_types;
    int num_bonus_flags;
    std::vector<std::vector<int>> hole;
    std::vector<int> hole_params; // (min_x, max_x, min_y, max_y, center_x, center_y, centroid_x, centroid_y)
    std::vector<std::vector<int>> vertices;
    std::vector<std::vector<int>> vertices0;
    //std::vector<int> vertices_memo;
    std::vector<int> figure_params; // (min_x, max_x, min_y, max_y, center_x, center_y, centroid_x, centroid_y)
    std::vector<std::vector<double>> forces;
    std::vector<std::vector<int>> edges;
    std::vector<int> figure_distances;
    std::vector<std::pair<double, std::vector<std::vector<int>>>> answers;
    std::vector<std::vector<int>> vertex_candidates;
    std::vector<std::vector<int>> bonuses;  // bonus_type, problem_id, x, y
    std::vector<int> bonus_flags;
    std::vector<std::vector<int>> bonus_candidates;
    std::vector<std::pair<int, int>> point_in_hole;
    double eps;
    double time_limit;

public:
    BrainWallSolver(double time_limit) : time_limit(time_limit * 1000.0) {}
    ~BrainWallSolver(){};

    int dislike();
    double energy(bool use_collision=false);
    bool validate_distance(int e, double d);
    bool validate_distance(int e);
    int validation();
    bool crossing(const std::vector<int>& a, const std::vector<int>& b, const std::vector<int>& c, const std::vector<int>& d);
    bool crossing_hole(const std::vector<int>& a, const std::vector<int>& b);
    bool point_in_area(const std::vector<int>& p);
    bool point_in_area(int i);
    double dist_point_to_line(const std::vector<int>& p1, const std::vector<int>& p2, const std::vector<int> point);
    std::vector<double> nearest_point_on_line(const std::vector<int>& p1, const std::vector<int>& p2, const std::vector<int>& point);
    double nearest_point_on_hole(const std::vector<int>& point, std::vector<int>& ans);
    bool has_crossing_edges();

    void method01();
    void method02();
    void method03();
    void method04();
    void method05();
    void method06();
    void method07();
    void method08();
    void all_methods();
    void move();

    void init();
    void read();
    int run();
    void output();
    void output_json(bool has_last=false);
    void print();

    template <typename U, typename T>
    double dist(const std::vector<U>& p, const std::vector<T>& q) {
        double dx = p[0] - q[0];
        double dy = p[1] - q[1];
        return dx * dx + dy * dy;
    }
    bool in_hole(const std::vector<int>& p) {
        return std::binary_search(point_in_hole.begin(), point_in_hole.end(), std::make_pair(p[0], p[1]));
    }
    bool in_hole(int i) {
        return in_hole(vertices[i]);
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

    std::cin >> num_bonus_types;
    std::vector<int> bonus;
    int problem_id, px, py;
    std::string bonus_type;
    for (int i = 0; i < num_bonus_types; i++) {
        std::cin >> bonus_type >> problem_id >> px >> py;
        bonus = {BONUS_TYPES[bonus_type], problem_id, px, py};
        bonuses.push_back(bonus);
        bonus_candidates.push_back({px, py});
    }

    std::cin >> num_bonus_flags;
    for (int i = 0; i < num_bonus_flags; i++) {
        bonus_flags.push_back(BONUS_TYPES[bonus_type]);
    }
}

void BrainWallSolver::init() {
    read();

    vertices0 = vertices;
    //for (int i = 0; i < num_vertices; i++) vertices_memo.push_back(-1);

    for ([[maybe_unused]]int i = 0; i < num_vertices; i++)
        vertex_candidates.push_back(std::vector<int>());

    std::vector<std::pair<double, std::vector<int>>> fig_edges_length;
    for (int e = 0; e < num_edges; e++) {
        int i = edges[e][0];
        int j = edges[e][1];
        fig_edges_length.push_back(std::make_pair(dist(vertices[i], vertices[j]), std::vector<int>({i, j})));
    }
    std::vector<std::pair<double, std::vector<int>>> hole_edges_length;
    for (int i = 0; i < num_vertices_of_hole - 1; i++) {
        hole_edges_length.push_back(std::make_pair(dist(hole[i], hole[i+1]), std::vector<int>({i, i+1})));
    }
    hole_edges_length.push_back(std::make_pair(dist(hole[num_vertices_of_hole-1], hole[0]), std::vector<int>({0, num_vertices_of_hole-1})));

    for (int i = 0; i < num_edges; i++) {
        double fig_len = fig_edges_length[i].first;
        auto fig_ids = fig_edges_length[i].second;
        for (int j = 0; j < num_vertices_of_hole; j++) {
            double hole_len = hole_edges_length[j].first;
            auto hole_ids = hole_edges_length[j].second;
            //if (std::abs(fig_len - hole_len) < fig_len * eps) {
            if (std::abs(fig_len - hole_len) < 0.0001) {
                for (int k = 0; k < 2; k++) {
                    if (vertex_candidates[fig_ids[k]].empty()) {
                        vertex_candidates[fig_ids[k]].push_back(hole_ids[0]);
                        vertex_candidates[fig_ids[k]].push_back(hole_ids[1]);
                    } else {
                        vertex_candidates[fig_ids[k]].push_back(hole_ids[0]);
                        vertex_candidates[fig_ids[k]].push_back(hole_ids[1]);
                        //auto r = std::set_intersection(vertex_candidates[fig_ids[k]].begin(), vertex_candidates[fig_ids[k]].end(), hole_ids.begin(), hole_ids.end(), vertex_candidates[fig_ids[k]].begin());
                        //vertex_candidates[fig_ids[k]].erase(r,  vertex_candidates[fig_ids[k]].end());
                    }
                }
            }
        }
    }
    for (auto& p : vertex_candidates) {
        std::sort(p.begin(), p.end());
        auto r = std::unique(p.begin(), p.end());
        p.erase(r, p.end());
    }

    /// debug
    std::cerr << "vertex_candidates:" << std::endl;
    int idx = 0;
    for (const auto& p : vertex_candidates) {
        std::cerr << idx++ << " :";
        for (const auto& q : p) std::cerr << " " << q;
        std::cerr << std::endl;
    }

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

    for (int i = hole_params[0]; i <= hole_params[1]; i++) {
        for (int j = hole_params[2]; j <= hole_params[3]; j++) {
            //if (point_in_area(std::vector<int>({i, j}))) point_in_hole.push_back(std::make_pair(i, j));
            //if (point_in_area({i, j})) point_in_hole.push_back(std::make_pair(i, j));
            if (point_in_area({i, j})) point_in_hole.push_back(std::make_pair(i, j));
        }
    }
    std::sort(point_in_hole.begin(), point_in_hole.end());
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
bool BrainWallSolver::point_in_area(int i) {
    if (vertices_memo[i] >= 0) return vertices_memo[i] != 0;
    
    bool ret = point_in_area(vertices[i]);
    vertices_memo[i] = ret ? 1 : 0;

    return ret;
}
#endif

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

int BrainWallSolver::dislike() {
    int sum = 0;
    for (int i = 0; i < num_vertices_of_hole; i++) {
        int min_d = dist(vertices[0], hole[i]);
        for (int j = 0; j < num_vertices; j++) {
            int d = dist(vertices[j], hole[i]);
            if (d < min_d) min_d = d;
        }
        sum += min_d;
    }
    return sum;
}

double BrainWallSolver::energy(bool use_collision) {
    const double rate_distance = 0.5;
    const double rate_area = 10.0;
    const double rate_crossing = 10.0;

    double sum = 0.0;
    //for (int i = 0; i < num_vertices; i++) sum += ene[i];
    sum += dislike();
    for (int e = 0; e < num_edges; e++) {
        int i = edges[e][0];
        int j = edges[e][1];
        double d = dist(vertices[i], vertices[j]);
        double d0 = figure_distances[e];
        double d1 = std::abs(d - d0);
        if (!validate_distance(e, d)) sum += d1 * rate_distance;
    }

    if (use_collision) {
        int cnt_in_area = 0;
        std::vector<bool> is_point_in_area;
        for (int i = 0; i < num_vertices; i++) {
            //bool is_in = point_in_area(vertices[i]);
            //bool is_in = point_in_area(i);
            bool is_in = in_hole(i);
            is_point_in_area.push_back(is_in);
            if (!is_in) cnt_in_area++;
        }
        sum += cnt_in_area * rate_area;

        int cnt_crossing = 0;
        std::vector<bool> is_crossing_edges;
        for (int i = 0; i < num_edges; i++) {
            int a = edges[i][0];
            int b = edges[i][1];
            int x = (vertices[a][0] + vertices[b][0]) / 2;
            int y = (vertices[a][1] + vertices[b][1]) / 2;
            //if (is_point_in_area[a] && is_point_in_area[b] && point_in_area(std::vector<int>({x, y}))) {
            if (is_point_in_area[a] && is_point_in_area[b] && in_hole(std::vector<int>({x, y}))) {
            //if (is_point_in_area[a] && is_point_in_area[b]) {
                bool is_crossing = crossing_hole(vertices[a], vertices[b]);
                is_crossing_edges.push_back(is_crossing);
                if (is_crossing) cnt_crossing++;
            } else {
                is_crossing_edges.push_back(true);
                cnt_crossing++;
            }
        }
        sum += cnt_crossing * rate_crossing;
    }

    return sum;
}

// translation of the figure to the center
void BrainWallSolver::method01() {
    //int offset_x = hole_params[4] - figure_params[4]; // center_x
    //int offset_y = hole_params[5] - figure_params[5]; // center_y
    int offset_x = hole_params[4] - figure_params[6]; // centroid_x
    int offset_y = hole_params[5] - figure_params[7]; // centroid_y
    for (int i = 0; i < num_vertices; i++) {
        vertices[i][0] += offset_x;
        //vertices[i][0] += offset_x + 150; // for problem 105
        vertices[i][1] += offset_y;
        //vertices_memo[i] = -1;
    }
}

// vertices of the figure move to points of candidates (matching edges, bonus areas, and so on)
void BrainWallSolver::method02() {
    std::uniform_int_distribution<> rand(0.0, 1.0);
    for (int i = 0; i < static_cast<int>(vertex_candidates.size()); i++) {
        if (!vertex_candidates[i].empty()) {
            std::uniform_int_distribution<> randint(0, vertex_candidates[i].size() - 1);
            int idx = vertex_candidates[i][randint(mt)];
            vertices[i] = hole[idx];
            //vertices_memo[i] = -1;
        } else {
            if (rand(mt) < 0.5) {
                std::uniform_int_distribution<> randcand(0, bonus_candidates.size() - 1);
                vertices[i] = bonus_candidates[randcand(mt)];
                //vertices_memo[i] = -1;
            }
        }
    }
}

// vertices of the figure move to vertices of the hole
void BrainWallSolver::method03() {
    std::uniform_int_distribution<> randint(0, num_vertices_of_hole - 1);
    for (int i = 0; i < num_vertices_of_hole && i < num_vertices; i++) {
        vertices[i] = hole[randint(mt)];
    }
}

// all vertices of the figure move to around center of the hole
void BrainWallSolver::method04() {
    std::uniform_int_distribution<> randint(-10, 10);
    for (int i = 0; i < num_vertices; i++) {
        vertices[i][0] = hole_params[4] + randint(mt);
        vertices[i][1] = hole_params[5] + randint(mt);
        //vertices_memo[i] = -1;
    }
}

// vertices of the figure move to nearest vertices of the hole
void BrainWallSolver::method05() {
    for (int i = 0; i < num_vertices; i++) {
        //if (!point_in_area(vertices[i])) {
        //if (!point_in_area(i)) {
        if (!in_hole(i)) {
            std::vector<int> ans;
            nearest_point_on_hole(vertices[i], ans);
            vertices[i] = ans;
            //vertices_memo[i] = -1;
        }
    }
}

// vertices of the figure move to nearest edges of the hole
void BrainWallSolver::method06() {
    std::vector<std::vector<int>> nearest_points;
    //std::vector<bool> is_point_in_area;
    for (int i = 0; i < num_vertices; i++) {
        //bool in_area = point_in_area(vertices[i]);
        //bool in_area = point_in_area(i);
        bool in_area = in_hole(i);
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
        //vertices_memo[i] = -1;
    }
}

// same positions of two vertices of an edge is divided to vertices of the hole
void BrainWallSolver::method07() {
    for (int e = 0; e < num_edges; e++) {
        int i = edges[e][0];
        int j = edges[e][1];
        if (vertices[i] == vertices[j]) {
#if 0
            for (const auto& pos : hole) {
                if (vertices[j] != pos) {
                    vertices[j] = pos;
                    //vertices_memo[j] = -1;
                }
            }
#endif
#if 1
            //is_check = false;
            bool is_found = false;
            for (int x = std::max(vertices[i][0] - 5, hole_params[0]); x < std::min(vertices[i][0] + 5, hole_params[1]); x++) {
                for (int y = std::max(vertices[i][1] - 5, hole_params[2]); y < std::min(vertices[i][1] + 5, hole_params[3]); y++) {
                    auto pos = std::vector<int>({x, y});
                    //if (point_in_area(pos)) {
                    if (in_hole(pos)) {
                        vertices[j] = pos;
                        //vertices_memo[j] = -1;
                        is_found = true;
                        break;
                    }
                }
                if (is_found) break;
            }
#endif
        }
    }
}

// search nearest vertices of the hole and crossing edges
void BrainWallSolver::method08() {
    std::vector<std::vector<int>> nearest_points;
    std::vector<bool> is_point_in_area;
    for (int i = 0; i < num_vertices; i++) {
        //is_point_in_area.push_back(point_in_area(vertices[i]));
        //is_point_in_area.push_back(point_in_area(i));
        is_point_in_area.push_back(in_hole(i));
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

    /// processing...
}

void BrainWallSolver::all_methods() {
    //method01();
    method02();
    //method03();
    //method04();
    //method05();
    method06();
    //method07();
    //method08();
}

bool BrainWallSolver::has_crossing_edges() {
    for (int ii = 0; ii < num_edges; ii++) {
        int a = edges[ii][0];
        int b = edges[ii][1];
        if (crossing_hole(vertices[a], vertices[b])) return true;
        /*
        int x = (vertices[a][0] + vertices[b][0]) / 2;
        int y = (vertices[a][1] + vertices[b][1]) / 2;
        if (point_in_area(vertices[a]) && point_in_area(vertices[b]) && point_in_area(std::vector<int>({x, y})) {
            is_crossing = crossing_hole(vertices[a], vertices[b]);
            if (is_crossing) break;
        } else {
            is_crossing = true;
            break;
        }
        */
    }
    return false;
}

void BrainWallSolver::move() {
    std::uniform_int_distribution<> randint(0, num_vertices - 1);
    std::uniform_real_distribution<> rand(0.0, 1.0);
    bool has_opt = false;
    for (int ee = 0; ee < num_edges * 2; ee++) {
        int idx = ee % 2;
        int e = ee / 2;
        int i = edges[e][idx];
        int j = edges[e][1-idx];
        double d0 = figure_distances[e];
        double dist0 = std::sqrt(d0 * (1.0 + eps));
        int range_min_x = static_cast<int>(vertices[i][0] - dist0 - 0.5);
        int range_max_x = static_cast<int>(vertices[i][0] + dist0 + 0.5);
        int range_min_y = static_cast<int>(vertices[i][1] - dist0 - 0.5);
        int range_max_y = static_cast<int>(vertices[i][1] + dist0 + 0.5);
        range_min_x = std::max(range_min_x, hole_params[0]);
        range_max_x = std::min(range_max_x, hole_params[1]);
        range_min_y = std::max(range_min_y, hole_params[2]);
        range_max_y = std::min(range_max_y, hole_params[3]);
#if 0
        if (rand(mt) < 0.1) {
            int k = randint(mt);
            std::swap(vertices[j][0], vertices[k][0]);
            std::swap(vertices[j][1], vertices[k][1]);
            //std::swap(vertices[j], vertices[k]);
            //std::swap(vertices_memo[j], vertices_memo[k]);
        }
#endif
        std::vector<std::pair<double, std::vector<int>>> candidates;
        double d = dist(vertices[i], vertices[j]);

        if (!validate_distance(e, d)) {
            has_opt = true;
            for (int x = hole_params[0]; x <= hole_params[1]; x++) {
                if (x < range_min_x || x > range_max_x) continue;
                for (int y = hole_params[2]; y <= hole_params[3]; y++) {
                    if (y < range_min_y || y > range_max_y) continue;
                    auto pos = std::vector<int>({x, y});
                    auto pos2 = std::vector<int>({(vertices[i][0] + x) / 2, (vertices[i][1] + y) / 2});
                    d = dist(vertices[i], pos);
                    
                    //if (std::abs(d / d0 - 1.0) > eps || !point_in_area(pos) || crossing_hole(vertices[i], pos) || !point_in_area(pos2)) continue;
                    //if (!validate_distance(e, d) || !point_in_area(pos) || crossing_hole(vertices[i], pos) || !point_in_area(pos2)) continue;
                    if (!validate_distance(e, d) || !in_hole(pos) || crossing_hole(vertices[i], pos) || !in_hole(pos2)) continue;
                    //if (std::abs(d / d0 - 1.0) > eps || crossing_hole(vertices[i], pos) || !point_in_area(pos2)) continue;
                    //if (!point_in_area(pos) && rand(mt) < 0.95) continue;
                    
                    auto prev = vertices[j];
                    vertices[j] = pos;
                    candidates.push_back(std::make_pair(energy(), pos));
                    vertices[j] = prev;
                }
            }
            std::uniform_int_distribution<> randcd(0, candidates.size() - 1);
            std::uniform_int_distribution<> randhole(0, num_vertices_of_hole - 1);
            std::uniform_int_distribution<> randbin(0, 1);
            if (!candidates.empty()) {
                std::sort(candidates.begin(), candidates.end());
                if (rand(mt) < 0.1) {
                    vertices[j] = candidates[randcd(mt)].second;
                    //vertices_memo[j] = -1;
                } else if (rand(mt) < 0.05) {
                    //vertices[randbin(mt)] = hole[randhole(mt)];
                    vertices[j] = hole[randhole(mt)];
                    //vertices_memo[j] = -1;
                } else {
                    vertices[j] = candidates.front().second;
                    //vertices_memo[j] = -1;
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
                vertices = vertices0;
                //vertices_memo = std::vector<int>(num_vertices, -1);
                //all_methods();
                if (rand(mt) < 0.2) method01();
                if (rand(mt) < 0.1) method03();
                if (rand(mt) < 0.1) method04();
                if (rand(mt) < 0.1) method05();
                if (rand(mt) < 0.2) method06();
                if (rand(mt) < 0.1) method07();
                if (rand(mt) < 0.3) method02();
#if 1
                for (int i = 0; i < num_vertices / 5 + 1; i++) {
                    int k = randint(mt);
                    int l = randint(mt);
                    std::swap(vertices[k], vertices[l]);
                }
                //int k = randint(mt);
                //int l = randint(mt);
                //std::swap(vertices[j], vertices[k]);
                //std::swap(vertices[l][0], vertices[k][0]);
                //std::swap(vertices[l][1], vertices[k][1]);
                //std::swap(vertices[l], vertices[k]);
                //vertices[j] = hole[randhole(mt)];
                //vertices[j] = hole[randhole(mt)];
#endif            
            }
            //if (vertices[i] == vertices[j]) vertices[j][0] += 1;
        }

    }
    if (validation() == 0) {
        std::cerr << "*** found a pose!" << std::endl;
        //answers.push_back(make_pair(energy(true), vertices));
        answers.push_back(make_pair(dislike(), vertices));
    }

    if (!has_opt) {
        for (int i = 0; i < static_cast<int>(vertex_candidates.size()); i++) {
            if (!vertex_candidates[i].empty()) {
                std::uniform_int_distribution<> randint(0, vertex_candidates[i].size() - 1);
                int idx = vertex_candidates[i][randint(mt)];
                vertices[i] = hole[idx];
            } else {
                std::uniform_int_distribution<> randint(-10, 10);
                vertices[i][0] += randint(mt);
                vertices[i][0] = vertices[i][0] < 0 ? 0 : vertices[i][0];
                vertices[i][1] += randint(mt);
                vertices[i][1] = vertices[i][1] < 0 ? 0 : vertices[i][1];
                //vertices_memo[i] = -1;
            }
        }
#if 0
        //first_step();
        //second_step();
        all_methods();
        int range_vertex = 10;
        std::uniform_int_distribution<> randint(-range_vertex, range_vertex);
        //std::uniform_int_distribution<> randintx(hole_params[0], hole_params[1]);
        //std::uniform_int_distribution<> randinty(hole_params[2], hole_params[3]);
        std::vector<std::pair<double, std::vector<std::vector<int>>>> vertices_candidates;
        auto original_vertices = vertices;
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < num_vertices; j++) {
                vertices[j][0] += randint(mt);
                vertices[j][0] = vertices[j][0] < 0 ? 0 : vertices[j][0];
                vertices[j][1] += randint(mt);
                vertices[j][1] = vertices[j][1] < 0 ? 0 : vertices[j][1];
                //vertices_memo[j] = -1;
            }
            vertices_candidates.push_back(std::make_pair(energy(true), vertices));
            vertices = original_vertices;
        }
        std::sort(vertices_candidates.begin(), vertices_candidates.end());
        vertices = vertices_candidates.front().second;
        //vertices_memo = std::vector<int>(num_vertices, -1);
#endif
    }
}

bool BrainWallSolver::validate_distance(int e, double d) {
    double d0 = figure_distances[e];
    return std::abs(d / d0 - 1.0) <= eps;
}

bool BrainWallSolver::validate_distance(int e) {
    int i = edges[e][0];
    int j = edges[e][1];
    double d = dist(vertices[i], vertices[j]);
    return validate_distance(e, d);
}

int BrainWallSolver::validation() {
    int cnt_in_area = 0;
    int cnt_crossing = 0;
    int cnt_distance = 0;

    std::vector<bool> is_point_in_area;
    for (int i = 0; i < num_vertices; i++) {
        //bool is_in = point_in_area(vertices[i]);
        //bool is_in = point_in_area(i);
        bool is_in = in_hole(i);
        is_point_in_area.push_back(is_in);
        if (!is_in) cnt_in_area++;
    }

    std::vector<bool> is_crossing_edges;
    for (int i = 0; i < num_edges; i++) {
        int a = edges[i][0];
        int b = edges[i][1];
        int x = (vertices[a][0] + vertices[b][0]) / 2;
        int y = (vertices[a][1] + vertices[b][1]) / 2;
        //if (is_point_in_area[a] && is_point_in_area[b] && point_in_area(std::vector<int>({x, y}))) {
        if (is_point_in_area[a] && is_point_in_area[b] && in_hole(std::vector<int>({x, y}))) {
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
    //std::cerr << std::endl;
    std::cerr << "num_vertices_of_hole: " << num_vertices_of_hole << std::endl;
#if 0
    for (int i = 0; i < num_vertices_of_hole; i++) {
        for (int j = 0; j < 2; j++) std::cerr << " " << std::fixed << std::setw(10) << std::setprecision(6) << hole[i][j];
        std::cerr << std::endl;
    }
#endif
    std::cerr << "num_vertices: " << num_vertices << std::endl;
#if 0
    for (int i = 0; i < num_vertices; i++) {
        for (int j = 0; j < 2; j++) std::cerr << " " << std::fixed << std::setw(10) << std::setprecision(6) << vertices[i][j];
        std::cerr << std::endl;
    }
#endif
    std::cerr << "num_edges: " << num_edges << std::endl;
#if 0
    for (int i = 0; i < num_edges; i++) {
        for (int j = 0; j < 2; j++) std::cerr << " " << std::fixed << std::setw(10) << std::setprecision(6) << edges[i][j];
        std::cerr << std::endl;
    }
#endif
    std::cerr << "eps: " << eps << std::endl;
    //validation();
    //std::cerr << "Invalid count: " << validation() << std::endl;
    std::cerr << "num_answers: " << answers.size() << std::endl;
}

void BrainWallSolver::output() {
    if (!answers.empty()) {
        std::sort(answers.begin(), answers.end());
        int dislike = answers.front().first;
        auto ans = answers.front().second;
        for (int i = 0; i < num_vertices; i++) {
            std::cout << ans[i][0] << " " << ans[i][1] << std::endl;
        }
        std::cerr << "dislike: " << dislike << std::endl;
    }
}

void BrainWallSolver::output_json(bool has_last) {
    if (!answers.empty()) {
        std::sort(answers.begin(), answers.end());
        for (int i = 0; i < 5 && i < static_cast<int>(answers.size()); i++) {
            auto ans = answers[i].second;
            std::stringstream ss;
            ss << "best" << i + 1 << ".json";
            std::fstream fs(ss.str(), std::ios_base::out);
            fs << "{\"vertices\": [";
            for (int i = 0; i < num_vertices; i++) {
                fs << "[" << ans[i][0] << "," << ans[i][1] << "]" << (i < num_vertices - 1 ? "," : "");
            }
            fs << "]}" << std::endl;
            fs.close();
        }
    }
    if (has_last) {
        std::fstream fs("last.json", std::ios_base::out);
        fs << "{\"vertices\": [";
        for (int i = 0; i < num_vertices; i++) {
            fs << "[" << vertices[i][0] << "," << vertices[i][1] << "]" << (i < num_vertices - 1 ? "," : "");
        }
        fs << "]}" << std::endl;
        fs.close();
    }
}

int BrainWallSolver::run() {
    std::cerr << "Energy (-): " << std::fixed << std::setw(15) << std::setprecision(5) << energy() << std::endl;
    all_methods();
    std::cerr << "Energy (0): " << std::fixed << std::setw(15) << std::setprecision(5) << energy() << std::endl;

    int num_iters = 0;
    for (int i : std::views::iota(0) | std::views::take_while([this](int i){
            double elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time).count();
            return time_limit - elapsed_time > 0.0;
        })) {
        move();
        if (i % 100 == 0) {
            std::cerr << "Energy (" << std::setw(7) << i + 1 << "): ";
            std::cerr << std::fixed << std::setw(15) << std::setprecision(5) << energy() << std::endl;
        }
        num_iters = i;
    }
    
    return num_iters;
}

int main(int argc, char** argv) {
#if 0
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [num_iters=1]" << std::endl;
        return 1;
    }
#endif

    start_time = std::chrono::system_clock::now();
    //int num_iters = argc > 1 ? std::stol(argv[1]) : 10000;
    double time_limit = argc > 1 ? std::stof(argv[1]) : DEFAULT_TIME_LIMIT;
    //int num_iters = 100; // no use
    //std::string fname(argc > 2 ? argv[2]: "");
    //std::string fname("");

    BrainWallSolver solver(time_limit);

    solver.init();
    solver.print();
    int n = solver.run();
    std::cerr << "num_iters: " << n << std::endl;
    solver.print();
    solver.output();
    solver.output_json(true);

    return 0;
}
