#ifndef BRAIN_WALL_INITIALIZER_HPP
#define BRAIN_WALL_INITIALIZER_HPP

#include <random>
#include "brain_wall.hpp"

static std::vector<lc::Point> initialize_passthrough(const Problem& problem){
	return problem.vertices;
}

static std::vector<lc::Point> initialize_random(const Problem& problem){
	const auto& hole = problem.hole;

	const int n = hole.size();
	double min_x = hole[0].x, max_x = hole[0].x;
	double min_y = hole[0].y, max_y = hole[0].y;
	for(int i = 1; i < n; ++i){
		min_x = std::min(min_x, hole[i].x);
		max_x = std::max(max_x, hole[i].x);
		min_y = std::min(min_y, hole[i].y);
		max_y = std::max(max_y, hole[i].y);
	}

	std::default_random_engine engine;
	std::uniform_real_distribution<double> x_dist(min_x, max_x);
	std::uniform_real_distribution<double> y_dist(min_y, max_y);

	const int k = problem.vertices.size();
	std::vector<lc::Point> vertices(k);
	for(int i = 0; i < k; ++i){
		do {
			vertices[i] = lc::Point(x_dist(engine), y_dist(engine));
		} while(hole.contains(vertices[i]) < 0);
	}

	return vertices;
}

#endif
