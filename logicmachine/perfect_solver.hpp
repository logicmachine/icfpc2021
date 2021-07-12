#ifndef BRAIN_WALL_PERFECT_SOLVER_HPP
#define BRAIN_WALL_PERFECT_SOLVER_HPP

#include <vector>
#include <climits>
#include "brain_wall.hpp"
#include "length_corrector.hpp"
#include "xorshift128.hpp"

namespace perfect_solver {

struct Bounds {
	double min_x, max_x, min_y, max_y;
	Bounds()
		: min_x( std::numeric_limits<double>::infinity())
		, max_x(-std::numeric_limits<double>::infinity())
		, min_y( std::numeric_limits<double>::infinity())
		, max_y(-std::numeric_limits<double>::infinity())
	{ }
};

static std::vector<lc::Point> correct_edge_lengths(
	const std::vector<Edge>& edges,
	const std::vector<lc::Point>& coords,
	const std::vector<uint8_t>& fixed,
	double eps)
{
	static const int    MAX_ITERATIONS = 2000;
	static const double ALPHA          = 0.0001;

	const int n = coords.size();
	std::vector<lc::Point> current = coords;

	std::vector<lc::Point> delta(n);
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
		if(score == 0.0){ break; }
		for(int i = 0; i < n; ++i){
			current[i] += delta[i];
		}
	}

	return current;
}

static double optimize(
	std::vector<lc::Point>& vertices,
	const lc::Polygon& hole,
	const std::vector<Edge>& edges,
	const std::vector<int>& indices,
	const std::vector<double>& matrix,
	double eps)
{
	const int n = indices.size();
	const int m = hole.size();
	vertices.assign(n, lc::Point());
	std::vector<uint8_t> fixed(n);

	std::vector<Bounds> bounds(n);
	for(int i = 0; i < m; ++i){
		const int ii = indices[i];
		fixed[ii] = 1;
		vertices[ii] = hole[i];
		bounds[ii].min_x = bounds[ii].max_x = vertices[ii].x;
		bounds[ii].min_y = bounds[ii].max_y = vertices[ii].y;
		for(int j = m; j < n; ++j){
			const int jj = indices[j];
			const auto r = matrix[ii * n + jj] * 1.05;
			bounds[jj].min_x = std::min(bounds[jj].min_x, vertices[ii].x + r);
			bounds[jj].max_x = std::max(bounds[jj].max_x, vertices[ii].x - r);
			bounds[jj].min_y = std::min(bounds[jj].min_y, vertices[ii].y + r);
			bounds[jj].max_y = std::max(bounds[jj].max_y, vertices[ii].y - r);
		}
	}
	for(int i = 0; i < n; ++i){
		std::swap(bounds[i].min_x, bounds[i].max_x);
		std::swap(bounds[i].min_y, bounds[i].max_y);
	}

	double score = 0.0;
	for(int i = m; i < n; ++i){
		const int ii = indices[i];
		const auto& b = bounds[ii];
		if(b.min_x <= b.max_x && b.min_y <= b.max_y){
			bool found = false;
			for(int i = 0; i < 100; ++i){
				vertices[ii].x = b.min_x + floating_random() * (b.max_x - b.min_x);
				vertices[ii].y = b.min_y + floating_random() * (b.max_y - b.min_y);
				if(hole.contains(vertices[ii])){
					found = true;
					break;
				}
			}
			if(found){ continue; }
		}
		if(b.min_x <= b.max_x){
			vertices[ii].x = b.min_x + floating_random() * (b.max_x - b.min_x);
		}else{
			// score += (b.min_x - b.max_x) * (b.min_x - b.max_x);
			vertices[ii].x = (b.min_x + b.max_x) * 0.5;
		}
		if(b.min_y <= b.max_y){
			vertices[ii].y = b.min_y + floating_random() * (b.max_y - b.min_y);
		}else{
			// score += (b.min_y - b.max_y) * (b.min_y - b.max_y);
			vertices[ii].y = (b.min_y + b.max_y) * 0.5;
		}
	}
	vertices = correct_edge_lengths(edges, vertices, fixed, eps);

	for(const auto& e : edges){
		const lc::Segment s(vertices[e.a], vertices[e.b]);
		const auto expect = e.squared_length;
		const auto actual = s.squared_length();
		score +=  0.01 * std::max(std::abs(actual / expect - 1.0) - eps, 0.0);
		score += 10.0 * distance(hole, s);
	}

	return score;
}

}

static std::vector<lc::Point> solve_perfect(const Problem& problem){
	static const int    NUM_ITERATIONS  = 1000000;
	static const double MAX_TEMPERATURE = 10.0;
	static const double MIN_TEMPERATURE = 0.1;
	static const double EPSILON_SCALE   = 0.1;

	const auto& hole  = problem.hole;
	const auto& edges = problem.edges;
	const auto  eps   = problem.epsilon;
	const int n = problem.vertices.size();
	const int m = hole.size();

	std::vector<double> distance_matrix(n * n, std::numeric_limits<double>::infinity());
	for(const auto& e : edges){
		const auto d = std::sqrt(static_cast<double>(e.squared_length));
		distance_matrix[e.a * n + e.b] = d;
		distance_matrix[e.b * n + e.a] = d;
	}
	for(int k = 0; k < n; ++k){
		for(int i = 0; i < n; ++i){
			for(int j = 0; j < n; ++j){
				distance_matrix[i * n + j] = std::min(
					distance_matrix[i * n + j],
					distance_matrix[i * n + k] + distance_matrix[k * n + j]);
			}
		}
	}

	std::vector<int> current(n);
	for(int i = 0; i < n; ++i){ current[i] = i; }
	for(int i = n - 1; i > 0; --i){
		const int j = modulus_random(i + 1);
		std::swap(current[i], current[j]);
	}
	std::vector<lc::Point> best;
	double current_score = perfect_solver::optimize(best, hole, edges, current, distance_matrix, eps * EPSILON_SCALE);
	double best_score    = current_score;

#ifndef NO_PROGRESS
	std::cerr << "Optimize: " << best_score << std::endl;
#endif

	for(int iteration = 0; iteration < NUM_ITERATIONS; ++iteration){
		const double progress = static_cast<double>(iteration) / NUM_ITERATIONS;
		const double temperature = MAX_TEMPERATURE - (MAX_TEMPERATURE - MIN_TEMPERATURE) * progress;

		int swap_a, swap_b;
		swap_a = modulus_random(m);
		swap_b = modulus_random(n);
		std::swap(current[swap_a], current[swap_b]);

		std::vector<lc::Point> new_vertices;
		const auto new_score = perfect_solver::optimize(new_vertices, hole, edges, current, distance_matrix, eps * EPSILON_SCALE);

		const double e = std::exp((current_score - new_score) / temperature);
		if(floating_random() <= e){
			current_score = new_score;
			if(new_score < best_score){
				best_score = new_score;
				best       = new_vertices;
#ifndef NO_PROGRESS
				std::cerr << "Optimize: " << best_score << " (" << iteration << " / " << NUM_ITERATIONS << ")" << std::endl;
#endif
			}
		}else{
			std::swap(current[swap_a], current[swap_b]);
		}
	}
	return best;
}

#endif
