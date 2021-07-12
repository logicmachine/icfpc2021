#include <cstdio>
using namespace std;
#include "shr_visualizer.hpp"
#include "solver.hpp"
#include "random.hpp"

Graphics g;
Movie mov;

struct VerletSolver : public Solver {
    vector<Point> adjust_vertices(const vector<Point> &inputs) {
        const vector<Point> &ref = figures;
        vector<Point> new_answers = inputs;
        for(int u=0; u<N_figure; u++) {
            for(auto v : G[u]) {
                vector<Point> &cur = new_answers;
                double d1, d2, dx_n, dy_n;
                {
                    double x1, y1, x2, y2;
                    tie(x1, y1) = get_coor(ref[u]);
                    tie(x2, y2) = get_coor(ref[v]);
                    ll dx = x1 - x2;
                    ll dy = y1 - y2;
                    d1 = sqrt(dx*dx + dy*dy);
                }
                {
                    double x1, y1, x2, y2;
                    tie(x1, y1) = get_coor(cur[u]);
                    tie(x2, y2) = get_coor(cur[v]);
                    dx_n = x1 - x2;
                    dy_n = y1 - y2;
                    d2 = sqrt(dx_n*dx_n + dy_n*dy_n);
                }
                double diff = d2 - d1;
                if(keeping_dist(u, v, ref, cur)) continue;                
                double perc = d2 < 1e-7 ? 0 : diff / d2 / 2;
                double offx = dx_n * perc;
                double offy = dy_n * perc;
                if(diff < 1e-7) {
                    offx += rnd.NextInt(-1, +1);
                    offy += rnd.NextInt(-1, +1);
                }
                double x1, y1; tie(x1, y1) = get_coor(cur[u]);
                double x2, y2; tie(x2, y2) = get_coor(cur[v]);
                new_answers[u] = Point(x1-offx, y1-offy);
                new_answers[v] = Point(x2+offx, y2+offy);
            }
        }
        return new_answers;
    }

    void get_initial_solution(const vector<Point> &src,
                              vector<Point> &dst) {
        // 適当に平行移動させるのを繰り返して、最も良さそうなものを選ぶ
        dst = src;
        ll min_penalty = calc_penalty(dst);
        fprintf(stderr, "max_coord = %d, before initial penalty: %lld\n", max_coord, min_penalty);

        // src をワープさせるか dst をちょっとだけ動かすか
        const int stride = 1;
        for(int dx=-max_coord; dx<=max_coord; dx+=stride) {
            for(int dy=-max_coord; dy<=max_coord; dy+=stride) {
                vector<Point> n_dst = dst;
                for(int i=0; i<N_figure; i++) {
                    int x, y; tie(x, y) = get_coor_int(src[i]);
                    n_dst[i] = Point(x + dx, y + dy);
                }
                if(chmin(min_penalty, calc_penalty(n_dst))) {
                    dst = n_dst;
                }
            }
        }
        fprintf(stderr, "initial penalty: %lld\n", min_penalty);
    }
    
    void solve_001(const vector<Point>& src,
                   vector<Point>& dst,
                   bool& updated) {
        if(updated) return;
        if(rnd.NextDouble() >= 0.5) return;
        // その時点でヤバイ点 (外に出ており、なおかつ hole から遠い) を見つけて
        // 最も近い点の方向に動かすのを繰り返す
        vector< pair<int, int> > operate_cand;
        for(int i=0; i<N_figure; i++) {
            int x, y; tie(x, y) = get_coor_int(src[i]);
            if(inPolygon(Point(x, y), holes)) continue;
            ll min_distance, hole_index;
            tie(min_distance, hole_index) = nearest_hole(x, y);
            operate_cand.emplace_back(i, hole_index);
        }            
        
        int cand_i = -1;
        if(operate_cand.size()) cand_i = rnd.NextInt(0, operate_cand.size() - 1);
        if(cand_i >= 0) {
            int operate_i, hole_i;
            tie(operate_i, hole_i) = operate_cand[cand_i];
            // hole_i = rnd.NextInt(0, N_hole - 1);
            ll x1, y1; tie(x1, y1) = get_coor_int(src[operate_i]);
            ll x2, y2; tie(x2, y2) = get_coor_int(holes[hole_i]);
            ll dx = (x2 - x1), dy = (y2 - y1);                
            dst[operate_i] = Point(x1 + dx, y1 + dy);
            updated = true;
        }        
    }

    void solve_002(const vector<Point>& src,
                   vector<Point>& dst,
                   bool& updated) {
        if(updated) return;
        if(calc_penalty(src) > 0) return;
        if(rnd.NextDouble() >= 0.5) return;
        vector< tuple<ll, ll, ll> > joints;
        for(int i=0; i<N_hole; i++) {
            ll d, k; tie(d, k) = nearest_joint(i, src);
            joints.emplace_back(d, k, i);
        }
        // 降順にソート -> 一番やばい hole が最初にくる        
        sort(joints.rbegin(), joints.rend());

        assert(joints.size() > 0);
        int i = rnd.NextInt(0, joints.size() - 1);
        int joint_i = get<1>(joints[i]);
        int hole_i = get<2>(joints[i]);
        {            
            ll x1, y1; tie(x1, y1) = get_coor_int(src[joint_i]);
            ll x2, y2; tie(x2, y2) = get_coor_int(holes[hole_i]);
            ll dx = clamp(x2 - x1, -5LL, 5LL);
            ll dy = clamp(y2 - y1, -5LL, 5LL);
            dst[joint_i] = Point(x1 + dx, y1 + dy);
            updated = true;
        }
    }

    void solve_003(const vector<Point>& src,
                   vector<Point>& dst,
                   bool& updated) {
        return;
    }

    void solve_004(const vector<Point>& src,
                   vector<Point>& dst,
                   bool& updated) {
        if(updated) return;
        if(rnd.NextDouble() >= 0.5) return;
        int dx = rnd.NextInt(-5, 5), dy = rnd.NextInt(-5, -5);
        for(int i=0; i<N_figure; i++) {
            int x, y; tie(x, y) = get_coor_int(src[i]);
            dst[i] = Point(x + dx, y + dy);
            updated = true;
        }
    }

    void solve_005(const vector<Point>& src,
                   vector<Point>& dst,
                   bool& updated) {
        if(updated) return;
        int operate_i = rnd.NextInt(0, N_figure - 1);
        int x1, y1; tie(x1, y1) = get_coor_int(src[operate_i]);
        ll dx = rnd.NextInt(-5, +5), dy = rnd.NextInt(-5, +5);
        dst[operate_i] = Point(x1 + dx, y1 + dy);        
    }
    
    virtual void solve() {
        const ll START_TEMP = 1'000'000;
        const ll END_TEMP = 10'000;
        const int END_STEP = 100;
        
        const int WINDOW_SIZE = max_coord + 30;
        g.screen(WINDOW_SIZE, WINDOW_SIZE);

        get_initial_solution(figures, answers);
        ll min_score = calc_penalty(answers);
        ll last_score = min_score;

        for(int step=0; step<END_STEP; step++) {
            // 手元のビジュアライザ (一定のステップごと)
            if(step % 1 == 0) {
                g.clear();
                g.noFill();
                g.stroke(1, 0, 0);
                g.poly(holes);
                g.stroke(0, 0, 0);
                g.fill(0, 0, 0);
                for(int i=0; i<N_figure; i++) {
                    double x, y; tie(x, y) = get_coor(answers[i]);
                    answers[i] = Point(x, y);
                    g.text(to_string(i), x+2, y+2, 2);
                }
                string last_score_str = "";
                {
                    int c = 0;
                    ll s = last_score;
                    while(s>0) {
                        last_score_str += (char)('0' + s % 10);
                        s /= 10;
                        c++;
                        if(s > 0 and c % 3 == 0) last_score_str += ",";
                    }
                }
                reverse(last_score_str.begin(), last_score_str.end());
                g.text("last score = " + last_score_str, 100, 5, 4);
                g.noFill();
                for(int i=0; i<N_figure; i++) {
                    for(auto j : G[i]) {
                        int x1, y1; tie(x1, y1) = get_coor_int(answers[i]);
                        int x2, y2; tie(x2, y2) = get_coor_int(answers[j]);
                        g.line(x1, y1, x2, y2);
                    }
                }
                mov.addFrame(g);
            }
            
            // 近傍選択
            // solve_x で解が変化するなら solve_y は実行しない (x < y)
            vector<Point> new_answers = answers;
            bool changed = false;
            // はみ出しているところを選んで動かす
            solve_001(answers, new_answers, changed);
            // hole から近い点の集合を得て、hole からの距離が遠い順にソート
            solve_002(answers, new_answers, changed);
            // 内部にある joint を選んで hole に近い場所に動かす
            solve_003(answers, new_answers, changed);
            // 全体を動かす
            solve_004(answers, new_answers, changed);
            // どこか適当にうごかす？
            solve_005(answers, new_answers, changed);

            if(step % 1 == 0) {
                g.clear();
                g.noFill();
                g.stroke(1, 0, 0);
                g.poly(holes);
                g.stroke(0, 0, 0);
                g.fill(0, 0, 0);
                for(int i=0; i<N_figure; i++) {
                    double x, y; tie(x, y) = get_coor(new_answers[i]);
                    new_answers[i] = Point(x, y);
                    g.text(to_string(i), x+2, y+2, 2);
                }
                string last_score_str = "";
                {
                    int c = 0;
                    ll s = last_score;
                    while(s>0) {
                        last_score_str += (char)('0' + s % 10);
                        s /= 10;
                        c++;
                        if(s > 0 and c % 3 == 0) last_score_str += ",";
                    }
                }
                reverse(last_score_str.begin(), last_score_str.end());
                g.text("last score = " + last_score_str, 100, 5, 4);
                g.noFill();
                for(int i=0; i<N_figure; i++) {
                    for(auto j : G[i]) {
                        int x1, y1; tie(x1, y1) = get_coor_int(new_answers[i]);
                        int x2, y2; tie(x2, y2) = get_coor_int(new_answers[j]);
                        g.line(x1, y1, x2, y2);
                    }
                }
                mov.addFrame(g);
            }
            
            // verlet algorithm で長さをいい感じにする
            for(int iter=0; iter<100; iter++) {
                new_answers = adjust_vertices(new_answers);
            }

            // ペナルティ (小さいほうがいい)
            // 違反があるなら正の数・ないなら負の数
            ll score = calc_penalty(new_answers);
            ll delta = min_score - score;
            double ratio = 1.0 * step / END_STEP;
            double temp = START_TEMP + (END_TEMP - START_TEMP) * ratio;
            double prob = exp(delta / temp);
            if(prob > rnd.NextDouble()) {
                chmin(min_score, score);
                last_score = score;
                answers = new_answers;
            }
        }
        get_initial_solution(answers, answers);
    }
};

int main(int argc, char** argv) {
    VerletSolver solver;
    if(argc > 1) {
        int testcase_num = stoi(argv[1]);
    
        char in_buf[128], out_buf[128];
        sprintf(in_buf, "../problems/txt/input_%03d.in", testcase_num);
        sprintf(out_buf, "output/answer_%03d.out", testcase_num);
        string in_filename = in_buf;
        string out_filename = out_buf;
        solver.set_file(in_filename, out_filename, "json");
    }
    solver.run();
    string html = mov.dumpHtml(1);
    ofstream fout;
    fout.open("movie.html", ios::out);
    fout << html << endl;
    fout.close();
    return 0;
}
