#ifndef BRAIN_WALL_XORSHIFT128_HPP
#define BRAIN_WALL_XORSHIFT128_HPP

#include <chrono>
#include <climits>
#include <cstdint>
#include <unistd.h>

inline uint32_t xorshift128(){
	// static uint32_t x = 192479812u;
	// static uint32_t y = 784892731u;
	static uint32_t x = std::chrono::system_clock::now().time_since_epoch().count();
	static uint32_t y = getpid();
	static uint32_t z = 427398108u;
	static uint32_t w = 48382934u; 
	const uint32_t t = x ^ (x << 11);
	x = y; y = z; z = w;
	w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)); 
	return w;
}

inline uint32_t modulus_random(uint32_t mod){
	const auto t = static_cast<uint64_t>(xorshift128()) * mod;
	return static_cast<uint32_t>(t >> 32);
}

inline double floating_random(){
	static const double rcp = 1.0 / std::numeric_limits<uint32_t>::max();
	return xorshift128() * rcp;
}

#endif
