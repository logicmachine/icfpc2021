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
double vertex_out_of_bound_penalty_factor = 1e6;
double edge_out_of_bound_penalty_factor = 1e6;
double edge_distance_penalty_factor = 1e6;
double integer_penalty_factor = 1e-100;

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
  if (s * t > 0) return false;

  s = (c.first - d.first) * (a.second - c.second) - (c.second - d.second) * (a.first - c.first);
  t = (c.first - d.first) * (b.second - c.second) - (c.second - d.second) * (b.first - c.first);
  if (s * t > 0) return false;

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
      sort(intersections.begin(), intersections.end());
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

double calc_vertex_score(int v) {
  double score = 0.0;

  // check if the vertex is inside the hole
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
    sort(intersections.begin(), intersections.end());
    if (intersections[0] != vertex[v]) reverse(intersections.begin(), intersections.end());
    bool is_outside = abs(angle_sum) < 1e-4;
    for (int i = 0; i < (int) intersections.size() - 1; i++) {
      if (is_outside) score += two_point_square_distance(intersections[i], intersections[i+1]) * edge_out_of_bound_penalty_factor * 2.0; // penalty
      is_outside = !is_outside;
    }
  }

  // check if the edge satisfies distance requirements
  for (auto &u: edge[v]) {
    double dist = two_point_square_distance(vertex[v], vertex[u]);
    if (abs(dist / original_distance[v][u] - 1) > eps * 1e-6) {
      score += (abs(dist / original_distance[v][u] - 1) - eps * 1e-6) * edge_distance_penalty_factor * 2.0; // penalty
    }
  }

  // check for near-integer
  double fx = abs(vertex[v].first - round(vertex[v].first));
  double fy = abs(vertex[v].second - round(vertex[v].second));
  score += (fx + fy) * integer_penalty_factor;

  // calculate dislikes
  for (int i = 0; i < h; i++) {
    double dislikes = 1e300;
    for (int u = 0; u < n; u++) {
      dislikes = min(dislikes, two_point_square_distance(vertex[u], hole[i]));
    }
    score += dislikes;
  }

  return score;
}

void solve() {
  // preprocess
  original_distance = vector<vector<double> >(n, vector<double>(n, -1.0));
  for (int v = 0; v < n; v++) {
    for (auto &u: edge[v]) {
      original_distance[v][u] = two_point_square_distance(vertex[v], vertex[u]);
    }
  }

  // {
  //   double gx = 0.0;
  //   double gy = 0.0;
  //   for (int i = 0; i < h; i++) {
  //     gx += hole[i].first;
  //     gy += hole[i].second;
  //   }
  //   for (int v = 0; v < n; v++) {
  //     vertex[v].first = gx + theRandom.nextDouble() * 1.0 - 0.5;
  //     vertex[v].second = gy + theRandom.nextDouble() * 1.0 - 0.5;
  //   }
  // }

  xmin = (int) hole[0].first;
  xmax = (int) hole[0].first;
  ymin = (int) hole[0].second;
  ymax = (int) hole[0].second;
  for (int i = 1; i < h; i++) {
    xmin = min(xmin, (int) hole[0].first);
    xmax = max(xmax, (int) hole[0].first);
    ymin = min(ymin, (int) hole[0].second);
    ymax = max(ymax, (int) hole[0].second);
  }

  last_vertex = vertex;

  eps *= 0.1;

  theRandom.x = (ll) random_device()() << 32 | (ll) random_device()();

  // optimization
  {
    double score = calc_score();
    double last_score = score;
    double best_score = score;

    double base_temperature = 1.0e6;
    double temperature = base_temperature;
    double decay_rate = 4.0e-5;
    double time_limit = 9.950;
    int iter_count = 0;

    theTimer.restart();

    while (theTimer.time() < time_limit) {
      double roll = theRandom.nextDouble();
      if (roll < 0.20) {
        int v = theRandom.nextUIntMod(n);
        double dx = theRandom.nextDouble() * 4.0 - 2.0;
        double dy = theRandom.nextDouble() * 4.0 - 2.0;

        score -= calc_vertex_score(v);
        vertex[v].first += dx;
        vertex[v].second += dy;
        score += calc_vertex_score(v);

        #ifdef DEBUG
        if (iter_count % 10000 == 0) cerr << iter_count << " " << score << " " << last_score << " " << best_score << " " << temperature << " " << theTimer.time() << endl;
        #endif

        if (score <= last_score) {
          last_score = score;
          last_vertex = vertex;
          if (score < best_score) {
            best_score = score;
          }
        }
        else if (theRandom.nextDouble() < exp(double(last_score - score) / temperature)) { // accept
          last_score = score;
          last_vertex = vertex;
        }
        else { // rollback
          vertex = last_vertex;
          score = last_score;
        }
      }
      else if (roll < 1.00) {
        for (int v = 0; v < n; v++) {
          double dx = theRandom.nextDouble() * 4.0 - 2.0;
          double dy = theRandom.nextDouble() * 4.0 - 2.0;

          vertex[v].first += dx;
          vertex[v].second += dy;
        }

        score = calc_score();

        #ifdef DEBUG
        if (iter_count % 10000 == 0) cerr << iter_count << " " << score << " " << last_score << " " << best_score << " " << temperature << " " << theTimer.time() << endl;
        #endif

        if (score <= last_score) {
          last_score = score;
          last_vertex = vertex;
          if (score < best_score) {
            best_score = score;
          }
        }
        else if (theRandom.nextDouble() < exp(double(last_score - score) / temperature)) { // accept
          last_score = score;
          last_vertex = vertex;
        }
        else { // rollback
          vertex = last_vertex;
          score = last_score;
        }
      }

      iter_count++;
      temperature *= 1.0 - decay_rate;
      // temperature = base_temperature * ((time_limit - theTimer.time()) / time_limit);
    }

    cerr << "iter_count  = " << iter_count << endl;
    cerr << "temperature = " << temperature << endl;
    cerr << "score       = " << score << endl;
    cerr << "best_score  = " << best_score << endl;
  }

  for (int v = 0; v < n; v++) {
    vertex[v].first = round(vertex[v].first);
    vertex[v].second = round(vertex[v].second);
  }
  vertex_out_of_bound_penalty_factor = 1e9;
  edge_out_of_bound_penalty_factor = 1e9;
  edge_distance_penalty_factor = 1e6;

  ans = vertex;

  {
    double score = calc_score();
    double last_score = score;
    double best_score = score;

    double base_temperature = 1.0e3;
    double temperature = base_temperature;
    double decay_rate = 4.0e-6;
    double time_limit = 9.950;
    int iter_count = 0;

    theTimer.restart();

    while (theTimer.time() < time_limit) {
      double roll = theRandom.nextDouble();
      if (roll < 0.10) {
        int v = theRandom.nextUIntMod(n);
        double x = theRandom.nextUIntMod(xmax - xmin + 1) + xmin;
        double y = theRandom.nextUIntMod(ymax - ymin + 1) + ymin;

        auto vertex_v_before = vertex[v];

        vertex[v].first = x;
        vertex[v].second = y;
        score = calc_score();

        #ifdef DEBUG
        if (iter_count % 100000 == 0) cerr << iter_count << " " << score << " " << last_score << " " << best_score << " " << temperature << " " << theTimer.time() << endl;
        #endif

        if (score <= last_score) {
          last_score = score;
          last_vertex = vertex;
          if (score < best_score) {
            ans = vertex;
            best_score = score;
          }
        }
        else if (theRandom.nextDouble() < exp(double(last_score - score) / temperature)) { // accept
          last_score = score;
          last_vertex = vertex;
        }
        else { // rollback
          vertex[v] = vertex_v_before;
          score = last_score;
        }
      }
      else if (roll < 0.50) {
        int v = theRandom.nextUIntMod(n);
        double dx = theRandom.nextUIntMod(11) - 5;
        double dy = theRandom.nextUIntMod(11) - 5;

        vertex[v].first += dx;
        vertex[v].second += dy;
        score = calc_score();

        #ifdef DEBUG
        if (iter_count % 100000 == 0) cerr << iter_count << " " << score << " " << last_score << " " << best_score << " " << temperature << " " << theTimer.time() << endl;
        #endif

        if (score <= last_score) {
          last_score = score;
          last_vertex = vertex;
          if (score < best_score) {
            ans = vertex;
            best_score = score;
          }
        }
        else if (theRandom.nextDouble() < exp(double(last_score - score) / temperature)) { // accept
          last_score = score;
          last_vertex = vertex;
        }
        else { // rollback
          vertex[v].first -= dx;
          vertex[v].second -= dy;
          score = last_score;
        }
      }
      else if (roll < 0.90) {
        for (int v = 0; v < n; v++) {
          double dx = theRandom.nextUIntMod(3) - 1;
          double dy = theRandom.nextUIntMod(3) - 1;

          vertex[v].first += dx;
          vertex[v].second += dy;
        }

        score = calc_score();

        #ifdef DEBUG
        if (iter_count % 100000 == 0) cerr << iter_count << " " << score << " " << last_score << " " << best_score << " " << temperature << " " << theTimer.time() << endl;
        #endif

        if (score <= last_score) {
          last_score = score;
          last_vertex = vertex;
          if (score < best_score) {
            ans = vertex;
            best_score = score;
          }
        }
        else if (theRandom.nextDouble() < exp(double(last_score - score) / temperature)) { // accept
          last_score = score;
          last_vertex = vertex;
        }
        else { // rollback
          vertex = last_vertex;
          score = last_score;
        }
      }
      else if (roll < 1.00) {
        double dx = theRandom.nextUIntMod(11) - 5;
        double dy = theRandom.nextUIntMod(11) - 5;

        for (int v = 0; v < n; v++) {
          vertex[v].first += dx;
          vertex[v].second += dy;
        }

        score = calc_score();

        #ifdef DEBUG
        if (iter_count % 100000 == 0) cerr << iter_count << " " << score << " " << last_score << " " << best_score << " " << temperature << " " << theTimer.time() << endl;
        #endif

        if (score <= last_score) {
          last_score = score;
          last_vertex = vertex;
          if (score < best_score) {
            ans = vertex;
            best_score = score;
          }
        }
        else if (theRandom.nextDouble() < exp(double(last_score - score) / temperature)) { // accept
          last_score = score;
          last_vertex = vertex;
        }
        else { // rollback
          vertex = last_vertex;
          score = last_score;
        }
      }

      if (theRandom.nextDouble() < 1e-4) { // restore best
        vertex = ans;
        score = best_score;
        last_score = best_score;
      }

      iter_count++;
      temperature *= 1.0 - decay_rate;
      // temperature = base_temperature * ((time_limit - theTimer.time()) / time_limit);
    }

    cerr << "iter_count  = " << iter_count << endl;
    cerr << "temperature = " << temperature << endl;
    cerr << "score       = " << score << endl;
    cerr << "best_score  = " << best_score << endl;
  }

  {
    cerr << calc_score() << endl;
    double d = 0.0;
    for (int v = 0; v < n; v++) {
      cerr << calc_vertex_score(v) << endl;
      d += calc_vertex_score(v);
    }
    cerr << d << endl;
  }
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

  for (auto &x: vertex) cout << (int) x.first << " " << (int) x.second << endl;

  return 0;
}
