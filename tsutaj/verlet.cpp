#include <cstdio>
using namespace std;
#include "shr_visualizer.hpp"
#include "solver.hpp"

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
                double perc = diff / d2 / 2;
                double offx = dx_n * perc;
                double offy = dy_n * perc;
                if(diff < 1e-7) {
                    offx += (rand() % 3) - 1;
                    offy += (rand() % 3) - 1;
                }
                double x1, y1; tie(x1, y1) = get_coor(cur[u]);
                double x2, y2; tie(x2, y2) = get_coor(cur[v]);
                new_answers[u] = Point(x1-offx, y1-offy);
                new_answers[v] = Point(x2+offx, y2+offy);
            }
        }
        return new_answers;
    }

    virtual void solve() {
        const int WINDOW_SIZE = max_coord + 30;
        g.screen(WINDOW_SIZE, WINDOW_SIZE);
        answers = figures;
        ll min_score = calc_penalty(answers);

        for(int step=0; step<3000; step++) {
            if(step % 10 == 0) {
                g.clear();
                g.noFill();
                g.stroke(1, 0, 0);
                g.poly(holes);
                g.stroke(0, 0, 0);
                g.fill(0, 0, 0);
                for(int i=0; i<N_figure; i++) {
                    double x, y; tie(x, y) = get_coor(answers[i]);
                    /*
                    x = min(1.0 * WINDOW_SIZE, max(0.0, x));
                    y = min(1.0 * WINDOW_SIZE, max(0.0, y));
                    */
                    answers[i] = Point(x, y);
                    g.text(to_string(i), x+2, y+2, 2);
                }
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

            vector<Point> new_answers = answers;
            
            // その時点でヤバイ点 (外に出ており、なおかつ hole から遠い) を見つけて
            // 最も近い点の方向に動かすのを繰り返す
            vector< pair<int, int> > operate_cand;
            for(int i=0; i<N_figure; i++) {
                int x, y; tie(x, y) = get_coor_int(answers[i]);
                if(inPolygon(Point(x, y), holes)) continue;
                ll min_distance, hole_index;
                tie(min_distance, hole_index) = nearest_hole(x, y);
                operate_cand.emplace_back(i, hole_index);
            }            

            int cand_i = -1;
            if(operate_cand.size()) cand_i = rand() % operate_cand.size();
            if(cand_i >= 0) {
                int operate_i, hole_i;
                tie(operate_i, hole_i) = operate_cand[cand_i];
                hole_i = rand() % N_hole;
                ll x1, y1; tie(x1, y1) = get_coor_int(answers[operate_i]);
                ll x2, y2; tie(x2, y2) = get_coor_int(holes[hole_i]);
                int scale = rand() % 3 + 1;
                ll dx = (x2 - x1) * scale, dy = (y2 - y1) * scale;                
                new_answers[operate_i] = Point(x1 + dx, y1 + dy);
            }
            else {
                // どこか適当にうごかす？
                int operate_i = rand() % N_figure;
                int x1, y1; tie(x1, y1) = get_coor_int(answers[operate_i]);
                int hole_i = rand() % N_hole;
                int x2, y2; tie(x2, y2) = get_coor_int(holes[hole_i]);
                ll dx = (x2 - x1), dy = (y2 - y1);
                new_answers[operate_i] = Point(x1 + dx, y1 + dy);
            }
            
            for(int iter=0; iter<100; iter++) {
                new_answers = adjust_vertices(new_answers);
            }
            ll score = calc_penalty(new_answers);
            if(score < min_score) {
                min_score = score;
                answers = new_answers;
            }
        }
    }
};

int main() {
    VerletSolver solver;
    /*
    int testcase_num;
    scanf("%d", &testcase_num);
    
    char in_buf[128], out_buf[128];
    sprintf(in_buf, "../problems/txt/input_%03d.in", testcase_num);
    sprintf(out_buf, "output/answer_%03d.out", testcase_num);
    string in_filename = in_buf;
    string out_filename = out_buf;
    solver.set_file(in_filename, out_filename);
    */
    solver.run();
    string html = mov.dumpHtml(5);
    ofstream fout;
    fout.open("movie.html", ios::out);
    fout << html << endl;
    fout.close();
    return 0;
}
