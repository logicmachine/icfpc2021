#include <iostream>
#include <algorithm>
#include <random>

#include "brain_wall.hpp"
#include "initializer.hpp"
#include "optimizer.hpp"
#include "discretise.hpp"
#include "perfect_solver.hpp"
#include "recursive.hpp"

using namespace std;

int main(){
	const auto problem = Problem::from(cin);
	int n_bonus, n_flags;
	cin >> n_bonus;
	for(int i = 0; i < n_bonus; ++i){
		string _;
		cin >> _ >> _ >> _ >> _;
	}
	cin >> n_flags;
	if(n_flags > 0){ return 0; }

	const double eps = problem.epsilon;

#if 1
	auto coords = initialize_random(problem);
	coords = optimize(problem.hole, problem.edges, coords, eps);
#elif 0
	auto coords = solve_perfect(problem);
#else
	auto coords = solve_recursive(problem);
#endif
	coords = discretise(problem.hole, problem.edges, coords, eps);

#if 0
	bool first = true;
	std::cout << "{\"vertices\":[";
	for(const auto& p : coords){
		if(!first){ std::cout << ","; }
		first = false;
		std::cout << "[" << p.x << "," << p.y << "]";
	}
	std::cout << "]}" << std::endl;
#else
	for(const auto& p : coords){
		std::cout << static_cast<int>(p.x + 0.5) << " " << static_cast<int>(p.y + 0.5) << std::endl;
	}
#endif

	return 0;
}

