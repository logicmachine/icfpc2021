#ifndef BRAIN_WALL_BRAIN_WALL_HPP
#define BRAIN_WALL_BRAIN_WALL_HPP

#include <iostream>
#include <vector>
#include "libcomp/geometry/point.hpp"
#include "libcomp/geometry/polygon.hpp"

//---------------------------------------------------------------------------
// Problem definition
//---------------------------------------------------------------------------

struct Edge {
	int    a, b;
	double squared_length;

	Edge() : a(0), b(0), squared_length(0) { }
	Edge(int a, int b, double l) : a(a), b(b), squared_length(l) { }
};


struct Problem {
	lc::Polygon            hole;
	std::vector<lc::Point> vertices;
	std::vector<Edge>      edges;
	double                 epsilon;

	Problem() : hole(), vertices(), edges(), epsilon(0) { }

	static Problem from(std::istream& is){
		Problem ret;

		std::size_t hole_n;
		is >> hole_n;
		ret.hole = lc::Polygon(hole_n);
		for(std::size_t i = 0; i < hole_n; ++i){
			double x, y;
			is >> x >> y;
			ret.hole[i] = lc::Point(x, y);
		}

		std::size_t figure_n, figure_m;
		is >> figure_n >> figure_m;
		ret.vertices = std::vector<lc::Point>(figure_n);
		for(std::size_t i = 0; i < figure_n; ++i){
			double x, y;
			is >> x >> y;
			ret.vertices[i] = lc::Point(x, y);
		}
		ret.edges = std::vector<Edge>(figure_m);
		for(std::size_t i = 0; i < figure_m; ++i){
			int a, b;
			is >> a >> b;
			ret.edges[i] = Edge(a, b, (ret.vertices[a] - ret.vertices[b]).norm());
		}

		int epsilon;
		is >> epsilon;
		ret.epsilon = epsilon / 1e6;

		return ret;
	}
};

#endif
