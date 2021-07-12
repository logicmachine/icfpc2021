#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#include <bits/stdc++.h>
#include "geometry.hpp"
#include "random.hpp"
using ll = long long int;

const int INF = 1 << 29;
const ll LONGINF = 1LL << 60;

struct Solver {
    int N_hole, N_figure, M_figure;
    vector<Point> holes, figures, answers;
    vector< vector<int> > G;
    int eps_num, use_eps_num;
    FILE *fp_in, *fp_out;
    Rand rnd;
    string format;

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
        int sgx = pt.real() < 0 ? -1 : +1;
        int sgy = pt.imag() < 0 ? -1 : +1;
        int x = abs(pt.real()) + 0.5, y = abs(pt.imag()) + 0.5;
        return make_pair(sgx * x, sgy * y);
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
                             bool verbose = false,
                             bool system_test = false) {
        if(u == v) return true;
        int eps = (system_test ? eps_num : use_eps_num);
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
        if(val > eps * d1) {
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
                fprintf(stderr, "expected: %lld, found: %lld (float: %f), ratio: %.12f > %d\n", d1, d2, d_float, val_float, eps);
            }
            return false;
        }
        return true;
    }
    pair<ll, ll> nearest_joint(ll hole_i, const vector<Point> &inputs) {
        ll min_distance = LONGINF, joint_index = -1;
        ll xh, yh; tie(xh, yh) = get_coor_int(holes[hole_i]);
        for(int i=0; i<N_figure; i++) {
            ll xj, yj; tie(xj, yj) = get_coor_int(inputs[i]);
            if(chmin(min_distance, squared_dist(xh, yh, xj, yj))) joint_index = i;
        }
        return make_pair(min_distance, joint_index);
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
        ll m1 = (x1 - x2) * (x1 - x2);
        ll m2 = (y1 - y2) * (y1 - y2);
        return m1 + m2;
    }    
    Solver(long long int seed = 0) {
        fp_in = stdin;
        fp_out = stdout;
        format = "txt";
        init(seed);
    }
    void init(long long int seed = 0) {
        max_coord = 0;
        rnd = Rand(seed);
        G.clear();
    }
    void set_file(string in_filename, string out_filename, string n_format) {
        fprintf(stderr, "open %s\n", in_filename.c_str());
        fp_in = fopen(in_filename.c_str(), "r");
        if(fp_in == nullptr) assert(false);

        fp_out = fopen(out_filename.c_str(), "w");
        if(fp_out == nullptr) assert(false);

        format = n_format;
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
        use_eps_num = eps_num / 1;
        fprintf(stderr, "# max_coord = %d\n", max_coord);
        fprintf(stderr, "# eps = %d\n", eps_num);
    }

    void output() {
        if(format == "json") {
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
        else {
            for(int i=0; i<N_figure; i++) {
                int x, y; tie(x, y) = get_coor_int(answers[i]);
                fprintf(fp_out, "%d %d\n", x, y);
            }
        }
    }
    
    ll calc_penalty(const vector<Point> &inputs, bool verbose = false) {
        ll penalty = 0;
        
        // 端点が多角形の内部にあるか
        vector<bool> v_in_poly(N_figure);
        for(int i=0; i<N_figure; i++) {
            int x, y; tie(x, y) = get_coor_int(inputs[i]);
            v_in_poly[i] = inPolygon(Point(x, y), holes);
        }

        // [max: 1e7] 線分の長さは適切か
        {
            int cnt = 0;
            for(int i=0; i<N_figure; i++) {
                for(auto j : G[i]) {
                    if(!keeping_dist(i, j, figures, inputs)) {
                        if(verbose) {
                            fprintf(stderr, "segment penalty [%d, %d]\n", i, j);
                        }
                        cnt++;
                    }
                }
            }
            penalty += 1.0 * cnt / N_figure * (ll)1e7;
        }
        
        // [max: 1e7] 端点が内部にない場合はペナルティ (最も近い hole との距離)
        {
            const ll DIST_THRESH = 500;
            int cnt = count(v_in_poly.begin(), v_in_poly.end(), false);
            ll sum = 0, max_val = cnt * DIST_THRESH * DIST_THRESH;
            for(int i=0; i<N_figure; i++) {            
                if(v_in_poly[i]) continue;
                int x, y; tie(x, y) = get_coor_int(inputs[i]);
                ll min_distance = nearest_hole(x, y).first;
                if(verbose) {
                    fprintf(stderr, "distance penalty (%d): %lld\n", i, min_distance);
                }
                min_distance = min(min_distance, DIST_THRESH);
                sum += min_distance * min_distance;
            }
            if(cnt > 0) penalty += 1.0 * sum / max_val * (ll)1e7;
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
                // 中点が多角形の内部にある
                double xm = (x1 + x2) / 2.0, ym = (y1 + y2) / 2.0;
                if(!inPolygon(Point(xm, ym), holes)) {
                    cnt_intersect++;
                }                
            }
        }
        penalty += cnt_intersect;
        if(verbose) {
            fprintf(stderr, "intersection penalty: %d\n", cnt_intersect);
            fprintf(stderr, "# total penalty = %lld\n", penalty);
        }

        // 制約に違反してない
        if(penalty == 0) {
            for(int i=0; i<N_hole; i++) {
                ll d = (ll)1e8 - nearest_joint(i, inputs).first;
                assert(d >= 0);
                penalty -= d * d;
            }
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
                if(!keeping_dist(i, j, figures, answers, true, true)) goto END_EVALUATION;
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
                // 中点が多角形の内部にある
                double xm = (x1 + x2) / 2.0, ym = (y1 + y2) / 2.0;
                if(!inPolygon(Point(xm, ym), holes)) {
                    fprintf(stderr, "Error: middle point [%d, %d] is not in polygon\n", i, j);
                    goto END_EVALUATION;
                }
            }            
        }

        // 1000 * log_2(v * e * h) * dislikes
        for(int i=0; i<N_hole; i++) {
            ll d = nearest_joint(i, answers).first;
            dislikes += d;
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
        double score = evaluate();
        fprintf(stderr, "score: %.12f\n\n", score);
    }
};

#endif // !__SOLVER_HPP__

