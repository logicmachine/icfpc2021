#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <iterator>
#include <set>
#include <algorithm>
#include <queue>
#include <chrono>

using namespace std;

const double EPS = 1e-8;
const double INF = 1e12;
using ldouble = long double;
typedef complex<ldouble> P;
using point = P;
using polygon = vector<P>;

namespace std {
  bool operator < (const P& a, const P& b) {
    return real(a) != real(b) ? real(a) < real(b) : imag(a) < imag(b);
  }
}
ldouble cross(const P& a, const P& b) {
  return imag(conj(a)*b);
}
ldouble dot(const P& a, const P& b) {
  return real(conj(a)*b);
}

struct L : public vector<P> {
  L(const P &a, const P &b) {
    push_back(a); push_back(b);
  }
};
using line = L;

typedef vector<P> G;
struct C {
  P p; ldouble r;
  C(const P &p, ldouble r) : p(p), r(r) { }
};

int ccw(P a, P b, P c) {
  b -= a; c -= a;
  if (cross(b, c) > 0)   return +1;       // counter clockwise
  if (cross(b, c) < 0)   return -1;       // clockwise
  if (dot(b, c) < 0)     return +2;       // c--a--b on line
  if (norm(b) < norm(c)) return -2;       // a--b--c on line
  return 0;
}

bool intersectLL(const L &l, const L &m) {
  return abs(cross(l[1]-l[0], m[1]-m[0])) > EPS || // non-parallel
         abs(cross(l[1]-l[0], m[0]-l[0])) < EPS;   // same line
}
bool sameLL(const L &l, const L &m) {
  return abs(cross(l[1]-l[0], m[0]-l[0])) < EPS;
}
bool intersectLS(const L &l, const L &s) {
  return cross(l[1]-l[0], s[0]-l[0])*       // s[0] is left of l
         cross(l[1]-l[0], s[1]-l[0]) < EPS; // s[1] is right of l
}
bool intersectLP(const L &l, const P &p) {
  return abs(cross(l[1]-p, l[0]-p)) < EPS;
}
bool intersectSS(const L &s, const L &t) {
  return ccw(s[0],s[1],t[0])*ccw(s[0],s[1],t[1]) <= 0 &&
         ccw(t[0],t[1],s[0])*ccw(t[0],t[1],s[1]) <= 0;
}
// 交差している場合true (重なっている場合はfalse, 1点のみ接している場合はfalse)
bool intersectStrictlySS(const L &s, const L &t) {
  bool is_same_endpoint = (s[0] == t[0] || s[0] == t[1] || s[1] == t[0] || s[1] == t[1]);
  is_same_endpoint &= !((s[0] == t[0] && s[1] == t[0]) || (s[0] == t[1] && s[1] == t[0]));
  return intersectSS(s, t) && !sameLL(s, t) && !is_same_endpoint;
}
bool intersectSP(const L &s, const P &p) {
  return abs(s[0]-p)+abs(s[1]-p)-abs(s[1]-s[0]) < EPS; // triangle inequality
}

P projection(const L &l, const P &p) {
  ldouble t = dot(p-l[0], l[0]-l[1]) / norm(l[0]-l[1]);
  return l[0] + t*(l[0]-l[1]);
}
P reflection(const L &l, const P &p) {
  return p + (ldouble)2.0 * (projection(l, p) - p);
}
ldouble distanceLP(const L &l, const P &p) {
  return abs(p - projection(l, p));
}
ldouble distanceLL(const L &l, const L &m) {
  return intersectLL(l, m) ? 0 : distanceLP(l, m[0]);
}
ldouble distanceLS(const L &l, const L &s) {
  if (intersectLS(l, s)) return 0;
  return min(distanceLP(l, s[0]), distanceLP(l, s[1]));
}
ldouble distanceSP(const L &s, const P &p) {
  const P r = projection(s, p);
  if (intersectSP(s, r)) return abs(r - p);
  return min(abs(s[0] - p), abs(s[1] - p));
}
ldouble distanceSS(const L &s, const L &t) {
  if (intersectSS(s, t)) return 0;
  return min(min(distanceSP(s, t[0]), distanceSP(s, t[1])),
             min(distanceSP(t, s[0]), distanceSP(t, s[1])));
}
P crosspoint(const L &l, const L &m) {
  ldouble A = cross(l[1] - l[0], m[1] - m[0]);
  ldouble B = cross(l[1] - l[0], l[1] - m[0]);
  if (abs(A) < EPS && abs(B) < EPS) return m[0]; // same line
  if (abs(A) < EPS) assert(false); // !!!PRECONDITION NOT SATISFIED!!!
  return m[0] + B / A * (m[1] - m[0]);
}

#define curr(P, i) P[i]
#define next(P, i) P[(i+1)%P.size()]
enum { OUT, ON, IN };
int contains(const polygon& P, const point& p) {
  bool in = false;
  for (int i = 0; i < P.size(); ++i) {
    point a = curr(P,i) - p, b = next(P,i) - p;
    if (imag(a) > imag(b)) swap(a, b);
    if (imag(a) <= 0 && 0 < imag(b))
      if (cross(a, b) < 0) in = !in;
    if (cross(a, b) == 0 && dot(a, b) <= 0) return ON;
  }
  return in ? IN : OUT;
}
// centerを中心にsrcをtheta分回転
P rotate(const P& center, const P& src, double theta) {
  P p(cos(theta), sin(theta));
  return p * (src - center) + center;
}
ldouble distanceSqurePP(const P &p1, const P &p2) {
  return abs(p1 - p2) * abs(p1 - p2);
}

struct Edge {
  Edge(int to, ldouble distance) : to(to), distance(distance) {}
  int to;
  ldouble distance;
};

void OutputJson(const vector<P>& pose) {
  cout << "{\"vertices\":[";
  for (int i = 0; i < pose.size(); i++) {
    cout << "[" << int(real(pose[i])) << "," << int(imag(pose[i])) << "]";
    if (i != pose.size() - 1) cout << ",";
  }
  cout << "]}" << endl;
}

// 点pとholesの頂点の中で一番近い点との距離
double nearesetDistance(const P &p, const polygon& holes_poly) {
  ldouble min_d = 1e9;
  for (const auto& hole : holes_poly) {
    min_d = min(min_d, distanceSqurePP(p, hole));
  }
  return min_d;
}
// dislikesの計算
double calcDislikes(const vector<P>& pose, const vector<int>& fixed, const polygon& holes_poly) {
  double dislike = 0;
  for (const auto& hole : holes_poly) {
    ldouble min_d = 1e9;
    for (int i = 0; i < pose.size(); i++) {
      const auto& p = pose[i];
      // 決定済みのもののみ計算
      if (fixed[i]) min_d = min(min_d, distanceSqurePP(p, hole));
    }
    dislike += min_d;
  }
  return dislike;
}
// 評価値の計算
double calcEvaluation(const vector<P>& pose, const vector<int>& fixed, const polygon& holes_poly) {
  double dislike = 0;
  for (const auto& hole : holes_poly) {
    ldouble min_d = 1e9;
    for (int i = 0; i < pose.size(); i++) {
      const auto& p = pose[i];
      // 決定済みのもののみ計算
      if (fixed[i]) min_d = min(min_d, distanceSqurePP(p, hole));
    }
    dislike += min_d;
  }
  return dislike;
}

bool isIntersectPoly(const polygon& holes_poly, const P& p1, const P& p2) {
  L line(p1, p2);
  for (int i = 0; i < holes_poly.size(); i++) {
    L hole_line(holes_poly[i], holes_poly[(i + 1) % holes_poly.size()]);
    if (intersectStrictlySS(line, hole_line)) return true;
  }
  return false;
}

// 新たに点を追加する時に制約を満たすか確認
bool isValidAddPointToPose(const P& p, const int id,
                           const vector<P>& pose,
                           const vector<int>& fixed,
                           const vector<vector<Edge>>& edges,
                           const vector<vector<int>>& in_holes_map,
                           const polygon& holes_poly,
                           const int eps) {
  int x = (int)real(p), y = (int)imag(p);
  if (!in_holes_map[y][x]) return false; // 領域内でないと駄目

  // 接続している点について
  for (const auto& edge : edges[id]) {
    if (!fixed[edge.to]) continue;
    const auto& adjp = pose[edge.to]; // 隣接点
    ldouble diff = abs(distanceSqurePP(p, adjp) / edge.distance - 1); // 元々の距離と今の点から隣接点の距離

    // cerr << id << " " << edge.to << " : " << diff << " <= " << eps / 1000000. << " ";
    // cerr << p << " - " << adjp << " = " << distanceSqurePP(p, adjp) << " / " << edge.distance << endl;
    // if (isIntersectPoly(holes_poly, p, adjp)) cerr << "intersect fail." << endl;
    if (diff > eps / 1000000.) return false;
    if (isIntersectPoly(holes_poly, p, adjp)) return false; // 交差している場合駄目
  }
  return true;
}

vector<int> getToporogicalOrder(const vector<int>& init_order,
                                const int n_figure,
                                const vector<P>& vertices,
                                const int m_figure,
                                const vector<vector<Edge>>& edges)
{
  vector<int> order;
  set<int> remaining;
  for (int i = 0; i < n_figure; i++) remaining.insert(i);

  for (int ii = 0; ii < n_figure; ii++) {
    int i = init_order[ii];

    if (remaining.find(i) == remaining.end()) continue;
    queue<int> que;
    que.push(i);
    order.push_back(i); // 0スタート
    remaining.erase(i);
    while (!que.empty()) {
      int curr_id = que.front(); que.pop();
      for (const auto& edge : edges[curr_id]) {
        int next = edge.to;
        if (remaining.find(next) == remaining.end()) continue; // もうないのでいかない
        order.push_back(next); // つながっているものを追加していく
        remaining.erase(next);
        que.push(next);
      }
    }
  }
  assert(remaining.empty() && order.size() == n_figure);
  return order;
}

void validateEdges(const vector<vector<Edge>>& edges, const vector<P>& pose, const int eps) {
  for (int u = 0; u < edges.size(); u++) {
    for (const auto& edge : edges[u]) {
      int v = edge.to;
      if (v < u) continue;
      ldouble diff = abs(distanceSqurePP(pose[u], pose[v]) / edge.distance - 1); // 元々の距離と今の点から隣接点の距離
      cerr << u << "-" << v << " : " << diff << " <= " << eps / 1000000. << endl;
      assert(diff < eps);
    }
  }
}

vector<P> tyrMatchAllPattern(const int n_figure,
                        const vector<P>& vertices,
                        const int m_figure,
                        const vector<vector<Edge>>& edges,
                        const polygon& holes_poly,
                        const vector<vector<int>>& in_holes_map,
                        const int eps,
                        const int maxh, const int maxw)
{
  vector<P> min_pose;
  double min_dislikes = 1e18;
  
  const int match_num = min((int)holes_poly.size(), n_figure);

  vector<P> pose(n_figure);
  vector<int> order(n_figure);
  for (int i = 0; i < n_figure; i++) order[i] = i;

  // 対応付けを order で変えながら試す
  int tryCount = 0;
  int endCount = 1e9;
  do {
    // cerr << "tryCount : " << tryCount++ << endl;
    // order[i]番目のposeの頂点を穴に対応させる
    bool isValid = true;
    vector<int> fixed(n_figure); // 頂点の位置決定済み
    for (int i = 0; i < match_num; i++) {
      int id = order[i];
      P p = holes_poly[i];
      // 制約を満たさない場合終了
      if (!isValidAddPointToPose(p, id, pose, fixed, edges, in_holes_map, holes_poly, eps)) {
        isValid = false;
        break;
      }
      // 決定
      pose[id] = p;
      fixed[id] = 1;
      // cerr << id << " : " << p << endl;
    }
    if (!isValid) continue;
    // cerr << "mid OK" << endl;

    // ========================================
    // poseの数の方が多い場合，残りを適当に配置
    for (int i = match_num; i < n_figure; i++) {
      int id = order[i];
      P valid_p(-1, -1);
      for (int y = 0; y < maxh; y++) {
        for (int x = 0; x < maxw; x++) {
          P curr_p(x, y);
          // 制約を満たして配置可能
          if (isValidAddPointToPose(curr_p, id, pose, fixed, edges, in_holes_map, holes_poly, eps)) {
            valid_p = curr_p;
            // cerr << id << " : " << valid_p << endl;
            goto match_end;
          }
        }
      }
match_end:
      if (real(valid_p) != -1) {
        pose[id] = valid_p;
        fixed[id] = 1;
      }
      else {
        isValid = false;
        break;
      }
    }


    // ========================================
    // 実行可能解であれば評価値を計算
    if (isValid) {
      double dislikes = calcDislikes(pose, fixed, holes_poly);
      // 更新
      if (dislikes < min_dislikes) {
        min_dislikes = dislikes;
        min_pose = pose;
      }
    }
    if (tryCount == endCount) break;
  } while (next_permutation(order.begin(), order.end()));
  
  return min_pose;
}

void solve2(const int n_figure,
           const vector<P>& vertices,
           const int m_figure,
           const vector<vector<Edge>>& edges,
           const polygon& holes_poly,
           const vector<vector<int>>& in_holes_map,
           const int eps,
           const int maxh, const int maxw)
{
  if (n_figure >= holes_poly.size() && n_figure <= 12) {
    vector<P> pose = tyrMatchAllPattern(n_figure, vertices, m_figure, edges, holes_poly, in_holes_map, eps, maxh, maxw);
    if (pose.size() == n_figure) {
      OutputJson(pose);
      cerr << "OK" << endl;
    }
    else {
      cerr << "NG" << endl;
    }
  }
  else {
    cerr << "Constraints" << endl;
  }
}

void solve(const int n_figure,
           const vector<P>& vertices,
           const int m_figure,
           const vector<vector<Edge>>& edges,
           const polygon& holes_poly,
           const vector<vector<int>>& in_holes_map,
           const int eps,
           const int maxh, const int maxw,
           const int timeout_sec=10)
{
  auto start_time = chrono::system_clock::now();

  mt19937 engine(1);
  vector<P> pose = vertices;

  int max_t = 100;
  vector<int> order(n_figure);
  for (int i = 0; i < n_figure; i++) order[i] = i;

  vector<P> maxPose;
  int max_ok_num = -1;
  int ok_num = 0;
  double min_score = 1e18;

  for (int t = 0; t < max_t; t++) {
    shuffle(order.begin(), order.end(), engine);
    // 4回に1回は接続順に並び替える（幅優先順序の方が制約の依存関係順になるので良さそう？）
    if (t%4!=3) order = getToporogicalOrder(order, n_figure, vertices, m_figure, edges);
    
    ok_num = 0;
    vector<int> fixed(n_figure); // 位置が決定済みか
    for (int id : order) {
      double min_dislikes = 1e18;
      P min_p;
      fixed[id] = 1; // 一旦決まったことにする最後に失敗時は戻す
      for (int y = 0; y < maxh; y++) {
        for (int x = 0; x < maxw; x++) {
          P curr_p(x, y);
          pose[id] = curr_p; // 一旦決まったことにする最後に失敗時は戻す

          // 制約を満たして配置可能
          if (isValidAddPointToPose(curr_p, id, pose, fixed, edges, in_holes_map, holes_poly, eps)) {
            double dislikes = nearesetDistance(curr_p, holes_poly)
                              + calcDislikes(pose, fixed, holes_poly) * 0.001;
            // cerr << id << " : " << dislikes << endl;
            if (min_dislikes > dislikes) {
              min_dislikes = dislikes;
              min_p = curr_p;
            }
          }
        }
      }
      if (min_dislikes != 1e18) {
        // cerr << id << " : " << min_p << " / " << min_dislikes << endl;
        fixed[id] = 1;
        pose[id] = min_p;
        ok_num++;
        // cerr << id << ": OK" << endl;
      }
      else {
        fixed[id] = 0;
        pose[id] = P(-1, -1);
        // cerr << id << ": MISS" << endl;
        break;
      }
    } // poseの決定
    if (max_ok_num <= ok_num) {
      max_ok_num = ok_num;
    }
    if (ok_num == n_figure) {
      double score = calcDislikes(pose, fixed, holes_poly);
      if (score < min_score) {
        maxPose = pose;
        min_score = score;
      }
      if (score == 0) break;
    }

    auto end_time = chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end_time-start_time).count();
    if (elapsed > timeout_sec) break; // タイムアウトで終了
  }
  
  if (max_ok_num == n_figure) {
    cerr << "ok" << endl;
    cerr << "score: " << min_score << " | " << max_ok_num << " / " << n_figure << endl;
    // validateEdges(edges, maxPose, eps);
    OutputJson(maxPose);
  }
  else {
    cerr << "ng" << endl;
    cerr << max_ok_num << " / " << n_figure << endl;
    // OutputJson(maxPose);
  }
}

int main(void) {
  int n_hole;
  cin >> n_hole;
  
  vector<P> holes;
  polygon holes_poly;
  int maxx = -1, maxy = -1;
  for (int i = 0; i < n_hole; i++) {
    int x, y; cin >> x >> y;
    holes.emplace_back(x, y);
    holes_poly.emplace_back(x, y);

    maxx = max(maxx, x);
    maxy = max(maxy, y);

  }

  // cerr << "holes" << endl;
  // for (int i = 0; i < n_hole; i++) {
  //   for (int j = i + 1; j < n_hole; j++) {
  //     cerr << i << holes_poly[i] << " - " << j << holes_poly[j] << " : " << distanceSqurePP(holes_poly[i], holes_poly[j]) << endl; 
  //   }
  // }

  int n_figure, m_figure;
  cin >> n_figure >> m_figure;

  vector<P> vertices;
  for (int i = 0; i < n_figure; i++) {
    int x, y; cin >> x >> y;
    vertices.emplace_back(x, y);
  }

  vector<vector<Edge>> edges(n_figure);
  // cerr << "edges" << endl;
  for (int i = 0; i < m_figure; i++) {
    int u, v; cin >> u >> v;
    ldouble distance = distanceSqurePP(vertices[u], vertices[v]);
    edges[u].emplace_back(v, distance);
    edges[v].emplace_back(u, distance);
    // cerr << u << " - " << v << " : " << distance << endl;
  }

  int eps;
  cin >> eps;
  // cerr << "eps " << eps << endl;

  int buf = 10;
  const int maxw = maxx + buf;
  const int maxh = maxy + buf;
  vector<vector<int>> in_holes_map(maxh, vector<int>(maxw));
  for (int y = 0; y < maxh; y++) {
    for (int x = 0; x < maxw; x++) {
      point p(x, y);
      int is_contain = contains(holes_poly, p);
      if (is_contain != OUT) {
        in_holes_map[y][x] = 1;
        // cout << "o";
      }
      else {
        // cout << "_";
      }
    }
    // cout << endl;
  }

  // 制約を満たす中で一番評価値が高くなる場所に置く
  solve(n_figure, vertices, m_figure, edges, holes_poly, in_holes_map,
        eps, maxh, maxw);


  return 0;
}


