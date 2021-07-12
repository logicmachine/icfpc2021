#ifndef BRAIN_WALL_UTILITY_HPP
#define BRAIN_WALL_UTILITY_HPP

#include <vector>
#include "libcomp/geometry/intersect.hpp"
#include "libcomp/geometry/distance.hpp"
#include "libcomp/geometry/crossing_points.hpp"
#include "libcomp/geometry/polygon.hpp"

inline lc::Line bisector(const lc::Line& a, const lc::Line& b){
	// performance options: a and b are always crossing
	const auto c = lc::crossing_points(a, b)[0];
	return lc::Line(c, c + a.direction() + b.direction());
}

inline double distance(const lc::Polygon& poly, const lc::Segment& s){
	const size_t n = poly.size();
	std::vector<lc::Point> candidates = { s.a, s.b };
#if 1
	for(size_t i = 0; i < n; ++i){
		const auto& a = poly[i];
		const auto& b = poly[(i + 1) % n];
		const auto& c = poly[(i + 2) % n];
		const auto bs = bisector(lc::Line(b, a), lc::Line(b, c));
		const auto cp = lc::crossing_points(bs, s);
		if(!cp.empty()){ candidates.push_back(cp[0]); }
	}
#else
	for(size_t i = 0; i < n; ++i){
		for(size_t j = i + 1; j < n; ++j){
			const lc::Line l0(poly[i], poly[(i + 1) % n]);
			const lc::Line l1(poly[j], poly[(j + 1) % n]);
			if(lc::intersect(l0, l1)){
				const auto bs = bisector(l0, l1);
				const auto cp0 = lc::crossing_points(bs, s);
				if(!cp0.empty()){ candidates.push_back(cp0[0]); }
				const auto cp1 = lc::crossing_points(
					lc::Line(bs[0], bs[0] + (bs[1] - bs[0]).ortho()), s);
				if(!cp1.empty()){ candidates.push_back(cp1[0]); }
			}
		}
	}
#endif
	double result = 0.0;
	for(const auto& p : candidates){
		if(poly.contains(p) >= 0){ continue; }
		double cur = std::numeric_limits<double>::infinity();
		for(size_t i = 0; i < n; ++i){
			cur = std::min(cur, lc::distance(poly.side(i), p));
		}
		result = std::max(result, cur);
	}
	return result;
}

inline bool strictly_intersect(const lc::Segment& a, const lc::Segment& b){
	// Touch is not intersect
	if(std::abs(lc::ccw(a, b.a)) == 1){ return false; }
	if(std::abs(lc::ccw(a, b.b)) == 1){ return false; }
	if(std::abs(lc::ccw(b, a.a)) == 1){ return false; }
	if(std::abs(lc::ccw(b, a.b)) == 1){ return false; }
	return lc::intersect(a, b);
}

#endif
