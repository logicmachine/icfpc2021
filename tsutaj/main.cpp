#include <bits/stdc++.h>
using namespace std;
#include "solver.hpp"
using ll = long long int;

void chmax(int &A, int B) { A = max(A, B); }
void chmin(int &A, int B) { A = min(A, B); }

int main() {
    Solver solver;
    int overall_max_coord = 0, overall_min_eps = INF;
    int overall_N_hole = 0, overall_N_figure = 0;
    for(int i=1; i<=132; i++) {
        char in_buf[128], out_buf[128];
        sprintf(in_buf, "../problems/txt/input_%03d.in", i);
        sprintf(out_buf, "output/answer_%03d.out", i);
        string in_filename = in_buf;
        string out_filename = out_buf;
        solver.set_file(in_filename, out_filename, "json");
        solver.run();
        chmax(overall_N_hole, solver.N_hole);
        chmax(overall_N_figure, solver.N_figure);
        chmax(overall_max_coord, solver.max_coord);
        chmin(overall_min_eps, solver.eps_num);
    }
    fprintf(stderr, "# overall: max_coord = %d, min_eps = %d\n", overall_max_coord, overall_min_eps);
    fprintf(stderr, "# overall: max_N_hole = %d, max_N_figure = %d\n", overall_N_hole, overall_N_figure);
}
