#ifndef BRAIN_WALL_DISCRETISE_HPP
#define BRAIN_WALL_DISCRETISE_HPP

#include <vector>
#include <cmath>
#include <climits>

#include "brain_wall.hpp"
#include "xorshift128.hpp"
#include "utility.hpp"


static double discretise_evaluate(
	const lc::Polygon& hole,
	const std::vector<Edge>& edges,
	const std::vector<lc::Point>& coords,
	const std::vector<lc::Point>& deltas,
	double eps)
{
	static const double HOLE_PENALTY = 1e3;

	double score = 0.0;
	for(const auto& e : edges){
		// Length penalty
		const auto p0 = coords[e.a] + deltas[e.a];
		const auto p1 = coords[e.b] + deltas[e.b];
		const auto dx = p1.x - p0.x;
		const auto dy = p1.y - p0.y;
		const auto expect = e.squared_length;
		const auto actual = dx * dx + dy * dy;
		const auto r = actual / expect - 1.0;
		score += std::max(std::abs(r) - eps, 0.0);
		// Hole penalty
		const lc::Segment s(p0, p1);
		score += HOLE_PENALTY * distance(hole, s);
	}
	return score;
}

static std::vector<lc::Point> discretise(
	const lc::Polygon& hole,
	const std::vector<Edge>& edges,
	const std::vector<lc::Point>& values,
	double eps)
{
	static const int    NUM_ITERATIONS  = 200000;
	static const double MAX_TEMPERATURE = 0.00001;
	static const double MIN_TEMPERATURE = 0.000001;

	static const double NEIGHBOR_DX[] = { 1, 1, 0, -1, -1, -1, 0, 1 };
	static const double NEIGHBOR_DY[] = { 0, -1, -1, -1, 0, 1, 1, 1 };
	static const double DELTA_LIMIT   = 3.0;

	const int n = values.size();

	std::vector<lc::Point> discretised(n);
	for(int i = 0; i < n; ++i){
		const auto& p = values[i];
		discretised[i] = lc::Point(std::round(p.x), std::round(p.y));
	}

	std::vector<lc::Point> current(n);
	std::vector<lc::Point> best = current;
	double current_score = discretise_evaluate(hole, edges, discretised, current, eps);
	double best_score    = current_score;
	std::cerr << "Discretise: " << best_score << std::endl;

	for(int iteration = 0; iteration < NUM_ITERATIONS; ++iteration){
		const double progress = static_cast<double>(iteration) / NUM_ITERATIONS;
		const double temperature = MAX_TEMPERATURE - (MAX_TEMPERATURE - MIN_TEMPERATURE) * progress;

		const uint32_t random = modulus_random(n * 8);
		const int target    = random >> 3;
		const int direction = random &  7;

		auto new_deltas = current;
		new_deltas[target].x += NEIGHBOR_DX[direction];
		new_deltas[target].y += NEIGHBOR_DY[direction];
		if(new_deltas[target].x < -DELTA_LIMIT || DELTA_LIMIT < new_deltas[target].x){ continue; }
		if(new_deltas[target].y < -DELTA_LIMIT || DELTA_LIMIT < new_deltas[target].y){ continue; }

		const double new_score = discretise_evaluate(hole, edges, discretised, new_deltas, eps);
		const double e = std::exp((current_score - new_score) / temperature);
		if(xorshift128() <= e * std::numeric_limits<uint32_t>::max()){
			current_score = new_score;
			current       = std::move(new_deltas);
			if(current_score < best_score){
				best_score = current_score;
				best       = current;
				std::cerr << "Discretise: " << best_score << std::endl;
				if(best_score == 0.0){ break; }
			}
		}
	}

	for(int i = 0; i < n; ++i){
		discretised[i] += current[i];
	}
	return discretised;
}

#endif
