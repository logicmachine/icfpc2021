#ifndef BRAIN_WALL_LENGTH_CORRECTOR_HPP
#define BRAIN_WALL_LENGTH_CORRECTOR_HPP

#include "brain_wall.hpp"

static std::vector<lc::Point> correct_edge_lengths(
	const std::vector<Edge>& edges,
	const std::vector<lc::Point>& coords,
	double eps)
{
	static const int    MAX_ITERATIONS = 1000;
	static const double ALPHA          = 0.001;

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
			const auto d = (actual - expect) * ALPHA * s.direction();
			delta[e.a] += d;
			delta[e.b] -= d;
			score += std::max(std::abs(actual / expect - 1.0) - eps, 0.0);
		}
		if(score == 0.0){ break; }
		for(int i = 0; i < n; ++i){
			current[i] += delta[i];
		}
	}

	return current;
}

#endif
