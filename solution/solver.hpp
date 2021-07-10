#include <bits/stdc++.h>
#include "geometry.hpp"
using ll = long long int;

const int INF = 1 << 29;
const ll LONGINF = 1LL << 60;

struct Solver {
    int N_hole, N_figure, M_figure;
    vector<Point> holes, figures, answers;
    vector< vector<int> > G;
    int eps_num;
    FILE *fp_in, *fp_out;

    // for statistics
    int max_coord;

    template <typename Tp>
    inline bool chmax(Tp &A, Tp B) {
        if(A < B) {
            A = B;
            return true;
        }
        return false;
    }
    
    template <typename Tp>
    inline bool chmin(Tp &A, Tp B) {
        if(A > B) {
            A = B;
            return true;
        }
        return false;
    }
    
    inline pair<double, double> get_coor(Point pt) {
        double x = pt.real(), y = pt.imag();
        return make_pair(x, y);
    }
    inline pair<int, int> get_coor_int(Point pt) {
        int x = pt.real() + 0.5, y = pt.imag() + 0.5;
        return make_pair(x, y);
    }
    inline bool is_valid_segment(Point s1, Point s2, Point h1, Point h2) {
        bool b1 = isec_ss(s1, s2, h1, h2);
        bool b2 = isec_sp(h1, h2, s1);
        bool b3 = isec_sp(h1, h2, s2);
        if(b1 and !(b2 or b3)) return false;
        return true;
    }
    inline bool keeping_dist(int u, int v,
                             const vector<Point>& ref,
                             const vector<Point>& cur,
                             bool verbose = false) {
        if(u == v) return true;
        ll d1, d2;
        {
            int x1, y1, x2, y2;
            tie(x1, y1) = get_coor_int(ref[u]);
            tie(x2, y2) = get_coor_int(ref[v]);
            d1 = squared_dist(x1, y1, x2, y2);
        }
        {
            int x1, y1, x2, y2;
            tie(x1, y1) = get_coor_int(cur[u]);
            tie(x2, y2) = get_coor_int(cur[v]);
            d2 = squared_dist(x1, y1, x2, y2);
        }
        ll val = 1000000LL * abs(d2 - d1);
        if(val > eps_num * d1) {
            if(verbose) {
                fprintf(stderr, "Error: the length of segment [%d, %d]\n", u, v);
                double d_float;
                {
                    double x1_f, y1_f, x2_f, y2_f;
                    tie(x1_f, y1_f) = get_coor(cur[u]);
                    tie(x2_f, y2_f) = get_coor(cur[v]);
                    int x1_s, y1_s, x2_s, y2_s;
                    tie(x1_s, y1_s) = get_coor_int(cur[u]);
                    tie(x2_s, y2_s) = get_coor_int(cur[v]);
                    fprintf(stderr, "position %d: (%f, %f) [double] (%d, %d) [int]\n", u, x1_f, y1_f, x1_s, y1_s);
                    fprintf(stderr, "position %d: (%f, %f) [double] (%d, %d) [int]\n", v, x2_f, y2_f, x2_s, y2_s);
                    d_float = (x1_f - x2_f) * (x1_f - x2_f) + (y1_f - y2_f) * (y1_f - y2_f);
                }
                double val_float = 1000000.0 * abs(1.0 * d2 / d1 - 1);
                fprintf(stderr, "expected: %lld, found: %lld (float: %f), ratio: %.12f > %d\n", d1, d2, d_float, val_float, eps_num);
            }
            return false;
        }
        return true;
    }
    pair<ll, ll> nearest_hole(ll x, ll y) {
        ll min_distance = LONGINF, hole_index = -1;
        for(int i=0; i<N_hole; i++) {
            ll xh, yh; tie(xh, yh) = get_coor_int(holes[i]);
            if(chmin(min_distance, squared_dist(x, y, xh, yh))) hole_index = i;
        }
        return make_pair(min_distance, hole_index);
    }
    ll squared_dist(ll x1, ll y1, ll x2, ll y2) {
        return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    }    
    Solver() {
        fp_in = stdin;
        init();
    }
    void init() {
        max_coord = 0;
        G.clear();
    }
    void set_file(string in_filename, string out_filename) {
        fprintf(stderr, "open %s\n", in_filename.c_str());
        fp_in = fopen(in_filename.c_str(), "r");
        if(fp_in == nullptr) assert(false);

        fp_out = fopen(out_filename.c_str(), "w");
        if(fp_out == nullptr) assert(false);
        init();
    }
    
    void input() {
        fscanf(fp_in, "%d", &N_hole);
        holes.resize(N_hole);
        for(int i=0; i<N_hole; i++) {
            int x, y; fscanf(fp_in, "%d%d", &x, &y);
            chmax(max_coord, max(abs(x), abs(y)));
            holes[i] = Point(x, y);
        }
        fscanf(fp_in, "%d%d", &N_figure, &M_figure);
        figures.resize(N_figure);
        answers.resize(N_figure);
        G.resize(N_figure);        
        for(int i=0; i<N_figure; i++) {
            int x, y; fscanf(fp_in, "%d%d", &x, &y);
            chmax(max_coord, max(abs(x), abs(y)));
            figures[i] = Point(x, y);
        }
        for(int i=0; i<M_figure; i++) {
            int u, v; fscanf(fp_in, "%d%d", &u, &v);
            G[u].emplace_back(v);
            // G[v].emplace_back(u);
        }
        fscanf(fp_in, "%d", &eps_num);
        fprintf(stderr, "# max_coord = %d\n", max_coord);
        fprintf(stderr, "# eps = %d\n", eps_num);
    }

    void output() {
        fprintf(fp_out, "{\n");
        fprintf(fp_out, "  \"vertices\": [\n");
        for(int i=0; i<N_figure; i++) {
            int x, y; tie(x, y) = get_coor_int(answers[i]);
            fprintf(fp_out, "    [%d, %d]", x, y);
            if(i + 1 < N_figure) fprintf(fp_out, ",");
            fprintf(fp_out, "\n");
        }
        fprintf(fp_out, "  ]\n");
        fprintf(fp_out, "}\n");
    }
    
    ll calc_penalty(const vector<Point> &inputs, bool verbose = false) {
        ll penalty = 0;
        
        // 端点が多角形の内部にあるか
        vector<bool> v_in_poly(N_figure);
        for(int i=0; i<N_figure; i++) {
            int x, y; tie(x, y) = get_coor_int(inputs[i]);
            v_in_poly[i] = inPolygon(Point(x, y), holes);
        }

        // 線分の長さは適切か
        for(int i=0; i<N_figure; i++) {
            for(auto j : G[i]) {
                if(!keeping_dist(i, j, figures, inputs)) {
                    if(verbose) {
                        fprintf(stderr, "segment penalty [%d, %d]\n", i, j);
                    }
                    penalty += 1;
                }
            }
        }
        
        // 端点が内部にない場合はペナルティ (最も近い hole との距離)
        for(int i=0; i<N_figure; i++) {
            if(v_in_poly[i]) continue;
            int x, y; tie(x, y) = get_coor_int(inputs[i]);
            ll min_distance = nearest_hole(x, y).first;
            if(verbose) {
                fprintf(stderr, "distance penalty (%d): %lld\n", i, min_distance);
            }
            penalty += 5LL * min_distance * min_distance;
        }
        
        // 線分の両端点が多角形の内部にある場合
        // 線分がはみ出る場合はペナルティ
        int cnt_intersect = 0;
        for(int i=0; i<N_figure; i++) {
            if(!v_in_poly[i]) continue;
            int x1, y1; tie(x1, y1) = get_coor_int(inputs[i]);
            for(auto j : G[i]) {
                if(!v_in_poly[j]) continue;
                int x2, y2; tie(x2, y2) = get_coor_int(inputs[j]);
                bool ok = true;
                for(int k=0; k<N_hole; k++) {
                    int k1 = k, k2 = (k + 1) % N_hole;
                    if(!is_valid_segment(Point(x1, y1), Point(x2, y2), holes[k1], holes[k2])) {
                        ok = false;
                        break;
                    }
                }
                if(!ok) cnt_intersect++;
            }
        }
        penalty += cnt_intersect;
        if(verbose) {
            fprintf(stderr, "intersection penalty: %d\n", cnt_intersect);
            fprintf(stderr, "# total penalty = %lld\n", penalty);
        }
        return penalty;
    }
    
    // 1000 * log_2(v * e * h) * dislikes
    double evaluate() {
        double score = 0;
        ll dislikes = 0;

        // check distance
        for(int i=0; i<N_figure; i++) {
            for(auto j : G[i]) {
                if(!keeping_dist(i, j, figures, answers, true)) goto END_EVALUATION;
            }
        }

        // check if each segment is inside polygon
        for(int i=0; i<N_figure; i++) {
            // 端点が多角形の内部にある
            int x1, y1; tie(x1, y1) = get_coor_int(answers[i]);
            if(!inPolygon(Point(x1, y1), holes)) {
                fprintf(stderr, "Error: point %d is not in polygon\n", i);
                goto END_EVALUATION;
            }
            for(auto j : G[i]) {
                // figure の線分が holes の線分と交わらない
                int x2, y2; tie(x2, y2) = get_coor_int(answers[j]);
                for(int k=0; k<N_hole; k++) {
                    int k1 = k, k2 = (k + 1) % N_hole;                    
                    if(!is_valid_segment(Point(x1, y1), Point(x2, y2), holes[k1], holes[k2])) {
                        fprintf(stderr, "Error: segment [%d, %d] intersects [%d, %d]\n", i, j, k1, k2);
                        goto END_EVALUATION;
                    }
                }
            }            
        }

        // 1000 * log_2(v * e * h) * dislikes
        for(int i=0; i<N_hole; i++) {
            int x1, y1; tie(x1, y1) = get_coor_int(holes[i]);
            ll min_distance = LONGINF;
            for(int j=0; j<N_figure; j++) {
                int x2, y2; tie(x2, y2) = get_coor_int(answers[j]);
                ll d = squared_dist(x1, y1, x2, y2);
                chmin(min_distance, d);
            }
            dislikes += min_distance;
        }
        score = 1000 * log2(N_figure * M_figure * N_hole) * dislikes;

    END_EVALUATION:
        return score;
    }

    virtual void solve() {
    }
    
    void run() {
        input();
        solve();
        output();
        calc_penalty(answers, true);
        int score = evaluate();
        fprintf(stderr, "score: %d\n\n", score);
    }
};
