#pragma GCC optimize ("O3")
#pragma GCC optimize ("unroll-loops")
#pragma GCC target ("avx2")

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <queue>
#include <stack>
#include <unordered_set>
#include <random>
#include <cmath>
#include <cassert>

#include <x86intrin.h>

struct xorshift64 {
  unsigned long long int x = 88172645463325252ULL;
  inline unsigned short nextUShort() {
    x = x ^ (x << 7);
    return x = x ^ (x >> 9);
  }
  inline unsigned int nextUShortMod(unsigned long long int mod) {
    x = x ^ (x << 7);
    x = x ^ (x >> 9);
    return ((x & 0x0000ffffffffffff) * mod) >> 48;
  }
  inline unsigned int nextUInt() {
    x = x ^ (x << 7);
    return x = x ^ (x >> 9);
  }
  inline unsigned int nextUIntMod(unsigned long long int mod) {
    x = x ^ (x << 7);
    x = x ^ (x >> 9);
    return ((x & 0x00000000ffffffff) * mod) >> 32;
  }
  inline unsigned long long int nextULL() {
    x = x ^ (x << 7);
    return x = x ^ (x >> 9);
  }
  inline double nextDouble() {
    x = x ^ (x << 7);
    x = x ^ (x >> 9);
    return (double)x * 5.42101086242752217e-20;
  }
};

struct timer {
  double t = 0.0;
  double lastStop = 0.0;
  bool stopped = false;
  timer() {
    restart();
  }
  inline void restart() {
    t = now();
    stopped = false;
  }
  inline void start() {
    if (stopped) {
      t += now() - lastStop;
      stopped = false;
    }
  }
  inline void stop() {
    if (!stopped) {
      lastStop = now();
      stopped = true;
    }
  }
  inline double time() {
    if (stopped) return lastStop - t;
    else return now() - t;
  }
  inline double now() {
    unsigned long long l, h;
    __asm__ ("rdtsc" : "=a"(l), "=d"(h));
    #ifdef LOCAL
    return (double)(l | h << 32) * 2.857142857142857e-10; // 1 / 3.5e9, for local (Ryzen 9 3950X)
    #else
    return (double)(l | h << 32) * 3.5714285714285715e-10; // 1 / 2.8e9, for AWS EC2 C3 (Xeon E5-2680 v2)
    //return (double)(l | h << 32) * 3.4482758620689656e-10; // 1 / 2.9e9, for AWS EC2 C4 (Xeon E5-2666 v3)
    //return (double)(l | h << 32) * 3.333333333333333e-10; // 1 / 3.0e9, for AWS EC2 C5 (Xeon Platinum 8124M / Xeon Platinum 8275CL)
    #endif
  }
};

using namespace std;

typedef long long int ll;
typedef pair<int, int> Pii;

const ll mod = 1000000007;

timer theTimer;
xorshift64 theRandom;
mt19937 theMersenne(1);

// inputs
int h;
vector<pair<double, double> > hole;
int n, m;
vector<pair<double, double> > vertex;
vector<vector<int> > edge;
double eps;

// outputs
vector<pair<double, double> > ans;

// internals
vector<vector<double> > original_distance;

vector<pair<double, double> > last_vertex;

int xmin, ymin, xmax, ymax;

// constants
constexpr double pi = 3.14159265358979323846264338327950;

// parameters
double vertex_out_of_bound_penalty_factor = 1e12;
double edge_out_of_bound_penalty_factor = 1e12;
double edge_distance_penalty_factor = 1e12;
double integer_penalty_factor = 1e12;

// 三角形pabの角度
double point_to_line_angle(pair<double, double> p, pair<double, double> a, pair<double, double> b) {
  double a1 = atan2(a.first - p.first, a.second - p.second);
  double a2 = atan2(b.first - p.first, b.second - p.second);
       if (abs(a2 - a1 + 2.0 * pi) < abs(a2 - a1)) return a2 - a1 + 2.0 * pi;
  else if (abs(a2 - a1 - 2.0 * pi) < abs(a2 - a1)) return a2 - a1 - 2.0 * pi;
  else                                             return a2 - a1;
}

// 線分と点の二乗距離
double point_to_line_square_distance(pair<double, double> p, pair<double, double> a, pair<double, double> b) {
  double x1 = b.first - a.first;
  double y1 = b.second - a.second;
  double x2 = x1 * x1;
  double y2 = y1 * y1;
  double r2 = x2 + y2;
  double tt = -(x1 * (a.first - p.first) + y1 * (a.second - p.second));

       if (tt <  0) return (a.first - p.first) * (a.first - p.first) + (a.second - p.second) * (a.second - p.second);
  else if (tt > r2) return (b.first - p.first) * (b.first - p.first) + (b.second - p.second) * (b.second - p.second);

  double f1 = x1 * (a.second - p.second) - y1 * (a.first - p.first);
  return (f1 * f1) / r2;
}

// 2点の二乗距離
double two_point_square_distance(pair<double, double> u, pair<double, double> v) {
  return (u.first - v.first) * (u.first - v.first) + (u.second - v.second) * (u.second - v.second);
}

// 2線分ab, cdの交差判定
bool is_two_line_intersect(pair<double, double> a, pair<double, double> b, pair<double, double> c, pair<double, double> d) {
  double s, t;
  s = (a.first - b.first) * (c.second - a.second) - (a.second - b.second) * (c.first - a.first);
  t = (a.first - b.first) * (d.second - a.second) - (a.second - b.second) * (d.first - a.first);
  if (s * t >= 0) return false;

  s = (c.first - d.first) * (a.second - c.second) - (c.second - d.second) * (a.first - c.first);
  t = (c.first - d.first) * (b.second - c.second) - (c.second - d.second) * (b.first - c.first);
  if (s * t >= 0) return false;

  return true;
}

// 2線分ab, cdの交点
pair<double, double> intersection_point_of_two_line(pair<double, double> a, pair<double, double> b, pair<double, double> c, pair<double, double> d) {
  double det = (a.first - b.first) * (d.second - c.second) - (d.first - c.first) * (a.second - b.second);
  double t = ((d.second - c.second) * (d.first - b.first) + (c.first - d.first) * (d.second - b.second)) / (det + 1e-18);
  double x = t * a.first + max(0.0, 1.0 - t) * b.first;
  double y = t * a.second + max(0.0, 1.0 - t) * b.second;
  return pair<double, double>(x, y);
}

double calc_score() {
  double score = 0.0;

  for (int v = 0; v < n; v++) {
    // check if vertices are inside the hole
    double angle_sum = 0.0;
    for (int i = 0; i < h; i++) {
      angle_sum += point_to_line_angle(vertex[v], hole[i], hole[(i+1)%h]);
    }
    if (abs(angle_sum) < 1e-4) { // outside
      double min_dist_to_edge = 1e300;
      for (int i = 0; i < h; i++) {
        min_dist_to_edge = min(min_dist_to_edge, point_to_line_square_distance(vertex[v], hole[i], hole[(i+1)%h]));
      }
      score += min_dist_to_edge * vertex_out_of_bound_penalty_factor; // penalty
    }

    // check if edges are inside the hole
    for (auto &u: edge[v]) {
      vector<pair<double, double> > intersections;
      for (int i = 0; i < h; i++) {
        if (is_two_line_intersect(vertex[v], vertex[u], hole[i], hole[(i+1)%h])) {
          intersections.push_back(intersection_point_of_two_line(vertex[v], vertex[u], hole[i], hole[(i+1)%h]));
        }
      }
      intersections.push_back(vertex[v]);
      intersections.push_back(vertex[u]);
      sort(intersections.begin(), intersections.end(), [&](auto &a, auto &b){return   two_point_square_distance(vertex[v], a) - two_point_square_distance(vertex[u], a)
                                                                                    < two_point_square_distance(vertex[v], b) - two_point_square_distance(vertex[u], b);});
      if (intersections[0] != vertex[v]) reverse(intersections.begin(), intersections.end());
      bool is_outside = abs(angle_sum) < 1e-4;
      for (int i = 0; i < (int) intersections.size() - 1; i++) {
        if (is_outside) score += two_point_square_distance(intersections[i], intersections[i+1]) * edge_out_of_bound_penalty_factor; // penalty
        is_outside = !is_outside;
      }
    }

    // check if the edge satisfies distance requirements
    for (auto &u: edge[v]) {
      double dist = two_point_square_distance(vertex[v], vertex[u]);
      if (abs(dist / original_distance[v][u] - 1) > eps * 1e-6) {
        score += (abs(dist / original_distance[v][u] - 1) - eps * 1e-6) * edge_distance_penalty_factor; // penalty
      }
    }

    // check for near-integer
    double fx = abs(vertex[v].first - round(vertex[v].first));
    double fy = abs(vertex[v].second - round(vertex[v].second));
    score += (fx + fy) * integer_penalty_factor;
  }

  // calculate dislikes
  for (int i = 0; i < h; i++) {
    double dislikes = 1e300;
    for (int v = 0; v < n; v++) {
      dislikes = min(dislikes, two_point_square_distance(vertex[v], hole[i]));
    }
    score += dislikes;
  }

  return score;
}

void solve() {
  original_distance = vector<vector<double> >(n, vector<double>(n, -1.0));
  for (int v = 0; v < n; v++) {
    for (auto &u: edge[v]) {
      original_distance[v][u] = two_point_square_distance(vertex[v], vertex[u]);
    }
  }

  xmin = (int) hole[0].first;
  xmax = (int) hole[0].first;
  ymin = (int) hole[0].second;
  ymax = (int) hole[0].second;
  for (int i = 1; i < h; i++) {
    xmin = min(xmin, (int) hole[i].first);
    xmax = max(xmax, (int) hole[i].first);
    ymin = min(ymin, (int) hole[i].second);
    ymax = max(ymax, (int) hole[i].second);
  }

  vector<vector<vector<pair<double, double> > > > possible_delta(n, vector<vector<pair<double, double> > >(n));
  for (int v = 0; v < n; v++) {
    for (auto &u: edge[v]) {
      if (original_distance[v][u] == -1.0) continue;
      for (double dx = -100; dx < 100; dx += 1.0) {
        for (double dy = -100; dy < 100; dy += 1.0) {
          double dist = dx * dx + dy * dy;
          if (abs(dist / original_distance[v][u] - 1.0) <= eps * 1e-6) possible_delta[v][u].emplace_back(dx, dy);
        }
      }
    }
  }
  // for (int v = 0; v < n; v++) {
  //   for (int u = v+1; u < n; u++) {
  //     cerr << "pd " << v << " " << u << endl;
  //     for (auto &x: possible_delta[v][u]) cerr << x.first << " " << x.second << endl;
  //   }
  // }
  // return;

  vector<Pii> check_tree;
  unordered_set<int> checked_vertex;
  checked_vertex.insert(0);
  for (int v = 0; v < n; v++) {
    for (auto &u: edge[v]) {
      if (checked_vertex.find(u) == checked_vertex.end()) {
        check_tree.emplace_back(v, u);
        checked_vertex.insert(u);
      }
    }
  }

  double best_score = 1e300;
  auto backtrack_search = [&](auto self, int e, unordered_set<int> &used_vertex) {
    if (e == n-1) {
      double score = calc_score();
      if (score < best_score) {
        best_score = score;
        ans = vertex;
      }
      return;
    }

    used_vertex.insert(check_tree[e].second);
    for (auto &d: possible_delta[check_tree[e].first][check_tree[e].second]) {
      vertex[check_tree[e].second].first = vertex[check_tree[e].first].first + d.first;
      vertex[check_tree[e].second].second = vertex[check_tree[e].first].second + d.second;
      if (vertex[check_tree[e].second].first < xmin || vertex[check_tree[e].second].first > xmax || vertex[check_tree[e].second].second < ymin || vertex[check_tree[e].second].second > ymax) continue;
      bool good = true;
      for (auto &t: edge[check_tree[e].second]) {
        if (used_vertex.find(t) == used_vertex.end()) continue;
        if (original_distance[check_tree[e].second][t] == -1.0) continue;
        double dist = two_point_square_distance(vertex[check_tree[e].second], vertex[t]);
        if (abs(dist / original_distance[check_tree[e].second][t] - 1.0) > eps * 8.0 * 1e-6) {
          good = false;
          break;
        }
      }
      if (good) {
        self(self, e+1, used_vertex);
      }
    }
    used_vertex.erase(check_tree[e].second);
  };

  unordered_set<int> used_vertex;
  used_vertex.insert(0);
  for (double ox = xmin; ox <= xmax; ox += 1.0) {
    for (double oy = ymin; oy <= ymax; oy += 1.0) {
      cerr << ox << " " << oy << endl;
      vertex[0].first = ox;
      vertex[0].second = oy;
      backtrack_search(backtrack_search, 0, used_vertex);
    }
  }

  cerr << best_score << endl;
}

int main() {
  cin.tie(0);
  ios::sync_with_stdio(false);

  // input
  cin >> h;
  hole = vector<pair<double, double> >(h);
  for (auto &x: hole) cin >> x.first >> x.second;
  cin >> n >> m;
  vertex = vector<pair<double, double> >(n);
  edge = vector<vector<int> >(n);
  for (auto &x: vertex) cin >> x.first >> x.second;
  for (int i = 0; i < m; i++) {
    int a, b;
    cin >> a >> b;
    edge[a].push_back(b);
    edge[b].push_back(a);
  }
  cin >> eps;

  // double angle_sum = 0.0;
  // for (int i = 0; i < h; i++) {
  //   angle_sum += point_to_line_angle(pair<double, double>(24.9, 35.0), hole[i], hole[(i+1)%h]);
  //   cerr << point_to_line_angle(pair<double, double>(24.9, 35.0), hole[i], hole[(i+1)%h]) << endl;
  // }
  // cerr << angle_sum << endl;
  // return 0;

  solve();

  for (auto &x: ans) cout << (int) x.first << " " << (int) x.second << endl;

  return 0;
}
