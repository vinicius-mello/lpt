#ifndef GLPT_HASH_UTIL_HPP
#define GLPT_HASH_UTIL_HPP

/*! \file
 * \brief Tiny shared hash-mixing/prime-sizing helpers for every open-
 * addressing table in this project (glpt_tree.hpp's leaf-existence set;
 * glpt_vertex_ids.hpp's edge cache; and, in vgtl2.0's riemann_cp2_glpt
 * example, glpt_crossing.hpp's face-crossing cache) -- factored out
 * 2026-09-23 when glpt_vertex_ids.hpp's std::map-based edge cache and
 * riemann_cp2_glpt's std::map-based face cache were found (by directly
 * measuring peak RSS against riemann_cp2.cpp) to be dominating memory
 * with std::map's own per-node overhead (3 pointers + color, on top of
 * the actual key/value), masking the real per-leaf saving glpt_tree
 * itself already delivers -- rather than duplicate glpt_tree.hpp's own
 * (already tested) hash/prime code a third time, both now share it.
 */

#include <cstdint>
#include <cstddef>

//! splitmix64-style avalanche mix, same as ~/code/lua/lpt/hash.c's
//! lpt_hash.
inline uint64_t glpt_hash64(uint64_t x) {
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
	x = (x ^ (x >> 27)) * 0x94d049b13c66a8edull;
	x = x ^ (x >> 31);
	return x;
}

inline bool glpt_is_prime(size_t n) {
	if(n<2) return false;
	if(n%2==0) return n==2;
	for(size_t d=3; d*d<=n; d+=2) if(n%d==0) return false;
	return true;
}
//! Smallest prime >= n (trial division -- cheap at the table sizes
//! these tables actually reach; ~/code/lua/lpt/hash.c's own
//! precomputed Mersenne-prime table is tuned for much larger ones).
inline size_t glpt_next_prime(size_t n) {
	if(n<3) return 3;
	size_t c = (n%2==0) ? (n+1) : n;
	while(!glpt_is_prime(c)) c += 2;
	return c;
}

#endif // GLPT_HASH_UTIL_HPP
