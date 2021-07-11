#include <iostream>
#include <algorithm>
#include <random>

#include "brain_wall.hpp"
#include "initializer.hpp"
#include "optimizer.hpp"
#include "discretise.hpp"

using namespace std;

int main(){
	const auto problem = Problem::from(cin);
	const double eps = problem.epsilon;

	auto coords = initialize_random(problem);
	coords = optimize(problem.hole, problem.edges, coords, eps);
	coords = discretise(problem.hole, problem.edges, coords, eps);

	bool first = true;
	std::cout << "{\"vertices\":[";
	for(const auto& p : coords){
		if(!first){ std::cout << ","; }
		first = false;
		std::cout << "[" << p.x << "," << p.y << "]";
	}
	std::cout << "]}" << std::endl;

	return 0;
}

