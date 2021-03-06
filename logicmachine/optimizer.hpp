#ifndef BRAIN_WALL_OPTIMIZER_HPP
#define BRAIN_WALL_OPTIMIZER_HPP

#include "brain_wall.hpp"
#include "length_corrector.hpp"
#include "xorshift128.hpp"
#include "utility.hpp"

static double optimize_evaluate(
	const lc::Polygon& hole,
	const std::vector<Edge>& edges,
	const std::vector<lc::Point>& vertices,
	double eps)
{
	static const double RCP_HOLE_PENALTY = 1e-4;

	double score = 0.0;
	for(const auto& e : edges){
		const lc::Segment s(vertices[e.a], vertices[e.b]);
		score += distance(hole, s);
	}
	const int n = hole.size();
	for(int i = 0; i < n; ++i){
		const auto& p = hole[i];
		double best = std::numeric_limits<double>::infinity();
		for(const auto& v : vertices){
			best = std::min(best, (p - v).norm());
		}
		score += best * RCP_HOLE_PENALTY;
	}
	return score;
}

static std::vector<lc::Point> optimize(
	const lc::Polygon& hole,
	const std::vector<Edge>& edges,
	const std::vector<lc::Point>& initial,
	double eps)
{
	static const int    NUM_ITERATIONS  = 20000; //100000;
	static const double MAX_TEMPERATURE = 0.1;
	static const double MIN_TEMPERATURE = 0.001;
	static const double ALPHA           = 1e3 * std::sqrt(eps);
	static const double EPSILON_SCALE   = 0.05;

	const int n = initial.size();
	std::vector<lc::Point> current = correct_edge_lengths(edges, initial, eps * EPSILON_SCALE);
	std::vector<lc::Point> best    = current;
	double current_score = optimize_evaluate(hole, edges, current, eps);
	double best_score    = current_score;
#ifndef NO_PROGRESS
	std::cerr << "Optimize: " << best_score << std::endl;
#endif

	for(int iteration = 0; iteration < NUM_ITERATIONS; ++iteration){
		const double progress = static_cast<double>(iteration) / NUM_ITERATIONS;
		const double temperature = MAX_TEMPERATURE - (MAX_TEMPERATURE - MIN_TEMPERATURE) * progress;
		const double alpha = std::max(ALPHA * (1.0 - progress), 0.1);

		const auto target = modulus_random(n);
		const auto dx = (floating_random() * 2.0 - 1.0) * alpha;
		const auto dy = (floating_random() * 2.0 - 1.0) * alpha;
		auto new_coords = current;
		new_coords[target].x += dx;
		new_coords[target].y += dy;
		new_coords = correct_edge_lengths(edges, new_coords, eps);

		const double new_score = optimize_evaluate(hole, edges, new_coords, eps * EPSILON_SCALE);
		const double e = std::exp((current_score - new_score) / temperature);
		if(floating_random() <= e){
			current_score = new_score;
			current       = std::move(new_coords);
			if(current_score < best_score){
				best_score = current_score;
				best       = current;
#ifndef NO_PROGRESS
				std::cerr << "Optimize: " << best_score << " (" << iteration << " / " << NUM_ITERATIONS << ")" << std::endl;
#endif
			}
		}
	}

	return best;
}

#endif
