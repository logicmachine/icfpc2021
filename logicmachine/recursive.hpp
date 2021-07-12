#ifndef BRAIN_WALL_RECURSIVE_HPP
#define BRAIN_WALL_RECURSIVE_HPP

#include <vector>
#include <climits>
#include "brain_wall.hpp"
#include "length_corrector.hpp"
#include "xorshift128.hpp"

class RecursiveSolver {

private:
	struct Bounds {
		double min_x, max_x, min_y, max_y;
		Bounds()
			: min_x( std::numeric_limits<double>::infinity())
			, max_x(-std::numeric_limits<double>::infinity())
			, min_y( std::numeric_limits<double>::infinity())
			, max_y(-std::numeric_limits<double>::infinity())
		{ }
	};

	lc::Polygon            hole;
	std::vector<lc::Point> original;
	std::vector<Edge>      edges;
	double                 eps;

	std::vector<double> raw_figure_matrix;
	std::vector<double> raw_hole_matrix;

	double figure_matrix(int i, int j) const { return raw_figure_matrix[i * original.size() + j]; }
	double hole_matrix(int i, int j) const { return raw_hole_matrix[i * hole.size() + j]; }

public:
	RecursiveSolver(const Problem& problem)
		: hole(problem.hole)
		, original(problem.vertices)
		, edges(problem.edges)
		, eps(problem.epsilon)
	{
		raw_figure_matrix = compute_figure_matrix();
		raw_hole_matrix   = compute_hole_matrix();
	}

	std::vector<lc::Point> solve(){
		const int n = original.size();
		std::vector<int> solution;
		std::vector<int> fixed(n);
		const int ret = recur(solution, fixed);
		std::cerr << "Ret: " << ret << std::endl;
		return {};
	}

private:
	int recur(
		std::vector<int>& solution,
		std::vector<int>& fixed)
	{
		const int n = original.size();
		const int m = hole.size();
		if(solution.size() == m){
			const auto bounds = compute_bounds(solution);
			if(test_corrected_lengths(fixed, bounds)){
				for(int i = 0; i < solution.size(); ++i){ std::cerr << solution[i] << " "; }
				std::cerr << std::endl;
				return 1;
			}else{
				return 0;
			}
		}
		int result = 0;
		for(int v = 0; v < n; ++v){
			if(fixed[v]){ continue; }
			if(!test_hole_distance(solution, v)){ continue; }
			solution.push_back(v);  fixed[v] = true;
			const auto bounds = compute_bounds(solution);
			if(!test_bounds(bounds)){
				solution.pop_back();  fixed[v] = false;
				continue;
			}
			result += recur(solution, fixed);
			solution.pop_back();  fixed[v] = false;
		}
		return result;
	}

	bool test_hole_distance(const std::vector<int>& solution, int v) const {
		const int n = solution.size();
		for(int i = 0; i < n; ++i){
			const auto f = figure_matrix(solution[i], v);
			const auto h = hole_matrix(i, n);
			if(f < h){ return false; }
		}
		return true;
	}

	std::vector<Bounds> compute_bounds(const std::vector<int>& solution) const {
		const int n = original.size();
		const int k = solution.size();
		std::vector<Bounds> bounds(n);
		for(int i = 0; i < k; ++i){
			const int  u = solution[i];
			const auto p = hole[i];
			for(int v = 0; v < n; ++v){
				const auto d = figure_matrix(u, v);
				bounds[v].min_x = std::min(bounds[v].min_x, p.x + d);
				bounds[v].max_x = std::max(bounds[v].min_x, p.x - d);
				bounds[v].min_y = std::min(bounds[v].min_y, p.y + d);
				bounds[v].max_y = std::max(bounds[v].min_y, p.y - d);
			}
		}
		return bounds;
	}

	bool test_bounds(const std::vector<Bounds>& bounds) const {
		for(const auto& b : bounds){
			if(b.min_x < b.max_x){ return false; }
			if(b.min_y < b.max_y){ return false; }
		}
		return true;
	}

	bool test_corrected_lengths(
		const std::vector<int>& fixed,
		const std::vector<Bounds>& bounds) const
	{
		static const int    MAX_ITERATIONS = 2000;
		static const double ALPHA          = 0.0001;
		const int n = fixed.size();
		std::vector<lc::Point> current(n);
		for(int i = 0; i < n; ++i){
			const auto& b = bounds[i];
			current[i].x = (b.min_x + b.max_x) * 0.5;
			current[i].y = (b.min_y + b.max_y) * 0.5;
		}
		std::vector<lc::Point> delta(n);
		double min_score = std::numeric_limits<double>::infinity();
		for(int iteration = 0; iteration < MAX_ITERATIONS; ++iteration){
			for(int i = 0; i < n; ++i){ delta[i] = lc::Point(0, 0); }
			double score = 0.0;
			for(const auto& e : edges){
				const lc::Segment s(current[e.a], current[e.b]);
				const auto expect = e.squared_length;
				const auto actual = s.squared_length();
				if(actual > 0.0){
					const auto d = (actual - expect) * ALPHA * s.direction();
					if(!fixed[e.a]){ delta[e.a] += d; }
					if(!fixed[e.b]){ delta[e.b] -= d; }
				}
				score += std::max(std::abs(actual / expect - 1.0) - eps, 0.0);
			}
			if(score == 0.0){ return true; }
			for(int i = 0; i < n; ++i){
				current[i] += delta[i];
			}
		}
		return false;
	}

	std::vector<double> compute_figure_matrix() const {
		const double margin = 1.05;
		const int n = original.size();
		std::vector<double> mat(n * n, std::numeric_limits<double>::infinity());
		for(int i = 0; i < n; ++i){ mat[i * n + i] = 0; }
		for(const auto& e : edges){
			const auto d = std::sqrt(static_cast<double>(e.squared_length));
			mat[e.a * n + e.b] = margin * d;
			mat[e.b * n + e.a] = margin * d;
		}
		for(int k = 0; k < n; ++k){
			for(int i = 0; i < n; ++i){
				for(int j = 0; j < n; ++j){
					mat[i * n + j] = std::min(mat[i * n + j], mat[i * n + k] + mat[k * n + j]);
				}
			}
		}
		return mat;
	}

	std::vector<double> compute_hole_matrix() const {
		const double margin = 1.05;
		const int n = hole.size();
		std::vector<double> mat(n * n, std::numeric_limits<double>::infinity());
		for(int i = 0; i < n; ++i){ mat[i * n + i] = 0; }
		for(int i = 0; i < n; ++i){
			for(int j = i + 1; j < n; ++j){
				const lc::Segment s(hole[i], hole[j]);
				bool accept = true;
				for(int k = 0; k < n; ++k){
					const int k0 = k, k1 = (k + 1) % n;
					if(i == k0 || i == k1 || j == k0 || j == k1){ continue; }
					if(lc::intersect(s, hole.side(k))){
						accept = false;
						break;
					}
				}
				lc::Polygon p(j - i + 1);
				for(int k = i; k <= j; ++k){ p[k - i] = hole[k]; }
				if(p.area() < 0.0){ accept = false; }
				if(accept){
					const auto d = s.length();
					mat[i * n + j] = margin * d;
					mat[j * n + i] = margin * d;
				}
			}
		}
		for(int k = 0; k < n; ++k){
			for(int i = 0; i < n; ++i){
				for(int j = 0; j < n; ++j){
					mat[i * n + j] = std::min(
						mat[i * n + j],
						mat[i * n + k] + mat[k * n + j]);
				}
			}
		}
		return mat;
	}

};

static std::vector<lc::Point> solve_recursive(const Problem& problem){
	RecursiveSolver solver(problem);
	return solver.solve();
}

#endif
