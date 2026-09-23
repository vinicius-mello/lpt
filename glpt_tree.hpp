#ifndef GLPT_TREE_HPP
#define GLPT_TREE_HPP

/*! \file
 * \brief The actual mesh: an open-addressing set of CURRENTLY-LEAF glpt
 * codes, deliberately NOT a port of lpt.hpp's own lpt_tree<lpt>.
 *
 * ============================================================
 * WHY NOT lpt_tree<glpt>
 * ============================================================
 * lpt_tree<lpt> (lpt.hpp) keeps EVERY cell ever created, leaf or not,
 * in one map forever -- simple_bisect() flips an existing entry's sign
 * to mark it "subdivided" and appends two new entries, but never erases
 * anything, and next_id only grows. That is exactly the append-only
 * memory-blowup pattern (see glpt.hpp's own WHY THIS FILE EXISTS
 * section) this whole effort exists to avoid -- reusing it would
 * silently reintroduce the problem with a smaller per-cell footprint,
 * not actually solve it. It's also not a drop-in template match for
 * glpt anyway (glpt::DIM vs. the expected lowercase dim; glpt::neighbor()
 * already fully resolves crossings and returns a `result` enum, not the
 * raw neighbor_bd() boundary-code lpt_tree's own stitching logic
 * expects to post-process).
 *
 * ============================================================
 * WHAT THIS FILE DOES INSTEAD
 * ============================================================
 * Stores ONLY the mesh's current leaves. Bisecting a leaf is a single
 * atomic operation -- erase(parent), insert(child(0)), insert(child(1))
 * -- so table occupancy is always exactly the current leaf count, never
 * historical operation count (see glpt_vertex_ids.hpp's own header
 * comment for the same growth argument applied to its edge cache).
 *
 * Storage: open addressing, linear probing, one uint64_t per slot --
 * bit 63 (glpt only uses 62 of 64 bits) is a PRESENT flag, so an empty
 * slot is plain 0 and the table starts life calloc-zeroed with no
 * separate "is this slot used" array needed (same trick
 * ~/code/lua/lpt/lpt_inc.c uses for its own PRESENT/LEAF bits, though
 * that file never actually deletes anything -- see above). Growth:
 * doubles (then rounds up to the next prime, via trial division -- the
 * table sizes here are expected to stay in the thousands-to-millions
 * range, where trial division up to sqrt(n) is cheap, unlike
 * ~/code/lua/lpt/hash.c's precomputed Mersenne-prime table tuned for
 * much larger sizes) once load factor exceeds 0.7. Deletion: backward-
 * shift (no tombstones) -- appropriate here specifically because every
 * bisection is one delete plus two inserts, so tombstones would
 * accumulate continuously and never get reclaimed except at a full
 * rehash; backward-shift keeps probe lengths bounded under that churn
 * instead. (Derived independently here rather than translated from a
 * reference implementation -- verified by glpt_tree_test.cpp's own
 * heavy insert/delete churn test, which checks every tracked leaf is
 * still exists()-reachable afterward, not just assumed correct.)
 *
 * ============================================================
 * NEIGHBOR-IN-TREE ("STITCHING") -- SCOPE, MATCHING lpt_tree's OWN
 * ============================================================
 * glpt::neighbor() gives the SAME-simplex_level candidate combinatorially,
 * but that candidate may not itself be a current leaf: the mesh may be
 * coarser there (some ANCESTOR of the candidate is the actual leaf) or
 * finer (the candidate has already been bisected, so some collection of
 * its DESCENDANTS are the actual leaves). locate_leaf() resolves the
 * first case completely (walk up via parent() until a leaf is found, or
 * the root is reached without one) -- this always terminates and is
 * unambiguous, since a cell has exactly one ancestor chain. The second
 * case is NOT resolved into an enumeration of possibly-many finer
 * leaves here: doing that correctly requires deciding which of a
 * bisected cell's children still touch a given (possibly deep) ancestor
 * facet, itself a real question about the specific vertex the bisected
 * edge introduces relative to that facet -- genuinely more machinery,
 * out of scope for this file, and NOT attempted by lpt_tree's own
 * neighbor_bd() either (it only special-cases exactly one level, and
 * only for two of the DIM+1 facet directions). neighbor_leaf() reports
 * this case explicitly (REFINED_FURTHER) rather than guessing or
 * silently returning something wrong.
 */

#include "glpt.hpp"
#include <cstdlib>

static inline uint64_t glpt_tree_hash(uint64_t x) {
	// Same splitmix64-style avalanche as ~/code/lua/lpt/hash.c's lpt_hash.
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
	x = (x ^ (x >> 27)) * 0x94d049b13c66a8edull;
	x = x ^ (x >> 31);
	return x;
}

static inline bool glpt_tree_is_prime(size_t n) {
	if(n<2) return false;
	if(n%2==0) return n==2;
	for(size_t d=3; d*d<=n; d+=2) if(n%d==0) return false;
	return true;
}
static inline size_t glpt_tree_next_prime(size_t n) {
	if(n<3) return 3;
	size_t c = (n%2==0) ? (n+1) : n;
	while(!glpt_tree_is_prime(c)) c += 2;
	return c;
}

class glpt_tree {
	public:
		static const uint64_t PRESENT_BIT = uint64_t(1)<<63;

		explicit glpt_tree(size_t initial_buckets = 257)
			: slots_(0), nbuckets_(0), count_(0)
		{
			alloc_(glpt_tree_next_prime(initial_buckets));
		}
		~glpt_tree() { std::free(slots_); }

		//! Seeds the mesh with all 108 Gaifullin root cells as leaves --
		//! the natural starting point for a whole-CP^2 mesh. Call once,
		//! before any bisect()/exists() (glpt_gaifullin_init() itself is
		//! idempotent and cheap, called here for convenience).
		void seed_all_roots() {
			glpt_gaifullin_init();
			for(int s=0;s<GLPT_NCELLS;++s) insert_(glpt(s).raw());
		}
		//! Seeds a single root (for tests, or a caller managing its own
		//! subset of seeds).
		void seed_root(int seed_index) { insert_(glpt(seed_index).raw()); }

		bool exists(const glpt& c) const { return find_slot_(c.raw()) != size_t(-1); }

		//! Bisects leaf `r` (which MUST currently exists()) into its two
		//! children, atomically: erase(r), insert(child(0)), insert(child(1)).
		void bisect(const glpt& r) {
			assert(exists(r) && "glpt_tree::bisect: r is not a current leaf");
			erase_(r.raw());
			insert_(r.child(0).raw());
			insert_(r.child(1).raw());
		}

		//! Walks up from `c` (checking `c` itself first) until a current
		//! leaf is found. Returns false iff no ancestor of `c`, including
		//! the root, is a leaf -- meaning `c`'s own region has already
		//! been refined AT LEAST as deep as `c` itself (see this file's
		//! own header comment on why that direction isn't resolved here).
		bool locate_leaf(const glpt& c, glpt& result) const {
			glpt cur = c;
			while(true) {
				size_t slot = find_slot_(cur.raw());
				if(slot != size_t(-1)) { result = cur; return true; }
				if(cur.is_root()) return false;
				cur = cur.parent();
			}
		}

		enum neighbor_status { FOUND, REFINED_FURTHER, UNRESOLVED };
		//! Same-level-or-coarser neighbor lookup: glpt::neighbor(r,i)
		//! gives the same-simplex_level candidate, then locate_leaf()
		//! walks up from it if needed. `r` itself need not be a leaf.
		neighbor_status neighbor_leaf(const glpt& r, int i, glpt& result) const {
			glpt n;
			glpt::result gr = r.neighbor(i, n);
			if(gr != glpt::OK) return UNRESOLVED;
			if(locate_leaf(n, result)) return FOUND;
			return REFINED_FURTHER;
		}

		size_t leaf_count() const { return count_; }
		size_t bucket_count() const { return nbuckets_; }
		double load_factor() const { return double(count_)/double(nbuckets_); }

		//! Visits every current leaf. F is called as f(const glpt&).
		template<class F> void for_each_leaf(F f) const {
			for(size_t i=0;i<nbuckets_;++i)
				if(slots_[i]!=0) f(glpt::from_raw(slots_[i] & ~PRESENT_BIT));
		}

	private:
		uint64_t* slots_;
		size_t nbuckets_;
		size_t count_;

		void alloc_(size_t nbuckets) {
			slots_ = (uint64_t*)std::calloc(nbuckets, sizeof(uint64_t));
			assert(slots_!=0 && "glpt_tree: out of memory");
			nbuckets_ = nbuckets;
		}

		size_t find_slot_(uint64_t code) const {
			size_t h = glpt_tree_hash(code) % nbuckets_;
			while(slots_[h]!=0) {
				if((slots_[h] & ~PRESENT_BIT) == code) return h;
				h = (h+1) % nbuckets_;
			}
			return size_t(-1);
		}

		void grow_if_needed_() {
			if(double(count_+1) <= 0.7*double(nbuckets_)) return;
			size_t old_n = nbuckets_;
			uint64_t* old = slots_;
			alloc_(glpt_tree_next_prime(2*old_n));
			size_t old_count = count_;
			count_ = 0;
			for(size_t i=0;i<old_n;++i)
				if(old[i]!=0) insert_(old[i] & ~PRESENT_BIT);
			(void)old_count;
			std::free(old);
		}

		void insert_(uint64_t code) {
			assert(find_slot_(code)==size_t(-1) && "glpt_tree: inserting a code that's already present");
			grow_if_needed_();
			size_t h = glpt_tree_hash(code) % nbuckets_;
			while(slots_[h]!=0) h = (h+1) % nbuckets_;
			slots_[h] = code | PRESENT_BIT;
			++count_;
		}

		// Backward-shift deletion: no tombstones. See this file's own
		// header comment for the derivation and why it's needed here.
		//
		// BUG FIXED 2026-09-23 (caught by glpt_tree_test.cpp's heavy-
		// churn test: leaf_count()/for_each_leaf() matched the ground
		// truth exactly, but exists() still failed on a genuinely-
		// present code -- for_each_leaf scans every slot directly, so
		// that split verdict pointed straight at a broken PROBE CHAIN,
		// not a missing entry). The first version stopped scanning as
		// soon as it found ONE entry that doesn't need to move back into
		// the gap -- wrong: entries further along the SAME probe run can
		// still depend on that gap being filled, even if the entry
		// immediately after it doesn't. Must keep scanning forward with
		// the gap held in place (examining each entry against the SAME
		// gap) until a genuinely empty slot ends the run, moving the gap
		// only when an entry actually needs to shift.
		void erase_(uint64_t code) {
			size_t gap = find_slot_(code);
			assert(gap!=size_t(-1) && "glpt_tree: erasing a code that isn't present");
			size_t scan = gap;
			while(true) {
				scan = (scan+1) % nbuckets_;
				if(slots_[scan]==0) break;
				uint64_t scan_code = slots_[scan] & ~PRESENT_BIT;
				size_t ideal = glpt_tree_hash(scan_code) % nbuckets_;
				if(cyclic_in_range_(gap, ideal, scan)) {
					// this entry's own probe path doesn't reach back
					// through gap -- leave it, keep scanning (gap stays).
					continue;
				}
				slots_[gap] = slots_[scan];
				gap = scan;
			}
			slots_[gap] = 0;
			--count_;
		}

		// Is x in the cyclic interval (a,b] modulo nbuckets_?
		static bool cyclic_in_range_impl_(size_t a, size_t x, size_t b) {
			if(a<=b) return x>a && x<=b;
			return x>a || x<=b;
		}
		bool cyclic_in_range_(size_t a, size_t x, size_t b) const {
			return cyclic_in_range_impl_(a,x,b);
		}

		glpt_tree(const glpt_tree&);
		glpt_tree& operator=(const glpt_tree&);
};

#endif // GLPT_TREE_HPP
