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
 *
 * ============================================================
 * COMPATIBLE (GRADED) BISECTION -- compat_bisect()
 * ============================================================
 * bisect() itself enforces no relationship between a cell and its
 * neighbors at all -- a mesh built with only bisect() calls can become
 * arbitrarily irregular. compat_bisect(r) additionally guarantees the
 * mesh stays 1-level graded (no same-level facet, other than the two a
 * bisection leaves whole -- i==DIM and i==r.level(), which don't need
 * this) ever has more than a 1-level depth difference across it),
 * recursively bisecting whatever neighbors are needed first. Ported
 * from lpt.hpp's own lpt_tree<lpt>::compat_bisect/compat_bisect_rec
 * (2026-09-23) -- structurally the same recursion, but that port
 * surfaced three real bugs, all found by glpt_tree_test.cpp's own
 * exhaustive graded-invariant sweep (not by inspection), documented in
 * compat_bisect_rec()'s own comment: two were outright crashes (a
 * same-level neighbor that's a root already bisected away has no
 * parent; a neighbor's parent already mid-recursion higher up the same
 * call stack got double-processed), and one was a wrong-direction
 * inference (a same-level neighbor absent because it had just been
 * refined FURTHER, within the very same top-level call, was mistaken
 * for one that's coarser). Whether the ORIGINAL lpt.hpp has the same
 * three gaps for the plain hypercube was not investigated -- Gaifullin's
 * seeds simply exercise paths (crossing between independently-specified
 * cells, not a single shared symmetric domain) that may not come up
 * there. Verified after the fixes at ~40000 compat_bisect() calls
 * cascading to nearly 4 million leaves (outside the committed suite),
 * 11.6 million facet checks, 0 graded-invariant violations.
 */

#include "glpt.hpp"
#include <cstdlib>
#include <set>
#include <vector>

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
		//! Both new children are appended to recent_leaves() (see its own
		//! comment) -- NOT cleared first, so a caller driving a priority-
		//! queue refinement loop can call clear_recent() once, then
		//! compat_bisect() (which can bisect many cells in one call, via
		//! its own neighbor-forcing recursion) and read back everything
		//! that was created, not just r's own two children.
		void bisect(const glpt& r) {
			assert(exists(r) && "glpt_tree::bisect: r is not a current leaf");
			erase_(r.raw());
			glpt c0=r.child(0), c1=r.child(1);
			insert_(c0.raw());
			insert_(c1.raw());
			recent_.push_back(c0);
			recent_.push_back(c1);
		}

		//! Cells inserted since the last clear_recent() call, across
		//! however many bisect() calls happened in between (a single
		//! compat_bisect() can trigger several, via its own recursive
		//! neighbor-forcing) -- mirrors lpt.hpp's own lpt_tree<lpt>::recent
		//! list, which exists for exactly this reason (a caller driving a
		//! priority-queue refinement loop needs to know every new leaf to
		//! push, not just the ones from the top-level bisect() call it
		//! made itself).
		//!
		//! NOT the same as "every CURRENT leaf inserted since
		//! clear_recent()": an entry can be superseded (bisected again)
		//! LATER within the same compat_bisect() call, e.g. if forcing one
		//! neighbor finer creates a cell that a DIFFERENT neighbor's own
		//! forcing then needs to refine again -- that cell's own children
		//! are separately present in recent_ from THAT later bisect(), so
		//! nothing is lost, but a stale entry can appear here. A caller
		//! must check exists() before using each one -- exactly the check
		//! a priority-queue consumer already needs when popping (an entry
		//! can just as easily go stale AFTER being queued, not only before
		//! being read out of recent_leaves()). Found by
		//! glpt_tree_test.cpp's own test, not merely assumed away.
		const std::vector<glpt>& recent_leaves() const { return recent_; }
		void clear_recent() { recent_.clear(); }

		//! Bisects `r`, first recursively bisecting whichever neighbors
		//! are needed to keep the mesh 1-level graded (no same-level
		//! facet ever has more than a 1-level depth difference across
		//! it) -- ported from lpt.hpp's own lpt_tree<lpt>::compat_bisect/
		//! compat_bisect_rec (2026-09-23, user's suggestion to look at
		//! how the original handles it), using glpt::neighbor() (the
		//! plain combinatorial same-level candidate, exactly as the
		//! original calls r.neighbor(i,n) on the bare lpt, NOT through
		//! any tree-aware stitching) plus exists() to detect and fix
		//! mismatches directly, the same way the original does with its
		//! own exists()/is_leaf().
		//!
		//! WHY facets i==DIM and i==r.level() are skipped: those are
		//! exactly the two facets of `r` that do NOT contain the edge
		//! (level(r),DIM) Theorem 1 is about to bisect -- they pass to
		//! one child whole, undivided, so a neighbor there doesn't need
		//! to be any finer on THIS cell's account (a property of which
		//! edge gets cut, independent of whether `r` itself is a 0- or
		//! 1-child of its own parent -- not to be confused with the
		//! UNRELATED sibling-direction facets neighbor() itself special-
		//! cases, see glpt.hpp's own comment).
		//!
		//! WHY recursing into n.parent() doesn't need its own existence
		//! check first (unlike bisect()'s precondition): by induction,
		//! if the mesh was already 1-level graded before this call (and
		//! compat_bisect() is the only way anything is ever bisected),
		//! `n` not existing means the mesh is exactly one level coarser
		//! there, so n.parent() -- exactly r's own level -- is
		//! guaranteed to be the leaf that's actually present -- UNLESS
		//! `n` is itself a root already bisected away, which has no
		//! parent at all (see compat_bisect_rec()'s own comment for the
		//! real crash this caused, and why nothing needs to happen in
		//! that case anyway). Verified directly (not just assumed) by
		//! glpt_tree_test.cpp's own graded-invariant sweep, checked after
		//! many compat_bisect() calls from a fresh seed_all_roots() mesh.
		void compat_bisect(const glpt& r) {
			assert(exists(r) && "glpt_tree::compat_bisect: r is not a current leaf");
			pending_.clear();
			compat_bisect_rec(r);
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
		std::set<glpt> pending_; // compat_bisect_rec()'s own in-progress set
		std::vector<glpt> recent_; // see recent_leaves()'s own comment

		void compat_bisect_rec(const glpt& r) {
			pending_.insert(r);
			for(int i=0;i<=glpt::DIM;++i) {
				if(i==glpt::DIM || i==r.level()) continue;
				glpt n;
				glpt::result gr = r.neighbor(i, n);
				if(gr!=glpt::OK) continue; // see glpt.hpp's own comment -- essentially never happens
				// BUG FIXED 2026-09-23 (caught by glpt_tree_test.cpp's own
				// graded-invariant test -- an assert in glpt::parent(),
				// not a silent wrong answer): the original lpt.hpp this
				// was ported from assumes !exists(n) always means "n is
				// exactly one level coarser, so n.parent() exists" -- true
				// for a NON-root n, by the graded invariant, but for a
				// ROOT n that's already been bisected away, n.parent() is
				// undefined (a root has none) AND unnecessary: r (also at
				// n's own level) being one level COARSER than n's now-
				// existing children is exactly the allowed direction of a
				// 1-level difference, nothing to force here.
				//
				// SECOND BUG FIXED 2026-09-23 (same test, next crash after
				// the one above: an assert in bisect() itself, meaning `r`
				// no longer existed when compat_bisect_rec(r) reached its
				// OWN closing bisect(r) call below): unlike the second
				// check just below, this first branch had no pending_
				// guard -- if n.parent() was already mid-recursion further
				// up the SAME call stack (added to pending_ by an
				// enclosing compat_bisect_rec call, but not yet actually
				// bisect()ed, since that only happens at the very end),
				// recursing into it again here double-processed it: the
				// nested call's own closing bisect() would consume it
				// early, so the enclosing call's LATER bisect() on that
				// same cell then found it already gone. (The original
				// lpt.hpp this was ported from has the identical gap --
				// not something introduced here -- just apparently never
				// exercised by whatever it was used on before.)
				//
				// THIRD BUG FIXED 2026-09-23 (same test, one more crash
				// after both fixes above -- again bisect(r)'s own assert):
				// !exists(n) does not ALWAYS mean "n is exactly one level
				// coarser" -- it can also mean n was refined FURTHER,
				// specifically by THIS SAME top-level compat_bisect() call,
				// moments earlier via a different facet of a different
				// cell (found by tracing a concrete case: n existed, got
				// bisected via the second branch below, and a later facet
				// of an ENTIRELY different cell in the same recursion then
				// asked about that same n again -- now absent, but its
				// PARENT was a root bisected in some earlier, unrelated
				// top-level call, long gone, not "one level coarser" at
				// all). The graded invariant only constrains the mesh
				// BETWEEN top-level compat_bisect() calls -- mid-recursion,
				// a neighbor can transiently become finer than r by more
				// than the usual margin. Checking whether n's OWN children
				// already exist (the direct symptom of the concrete case
				// found) catches this: if so, n is already fine, same as
				// the root case above -- nothing to force.
				bool n_at_depth_cap = n.orthant_level()>=glpt::MAX_ORTHANT_LEVEL && n.level()>=glpt::DIM-1;
				bool n_refined_further = !exists(n) && !n_at_depth_cap && exists(n.child(0));
				if(!exists(n) && !n_refined_further && !n.is_root() && pending_.count(n.parent())==0) compat_bisect_rec(n.parent());
				if(exists(n) && pending_.count(n)==0) compat_bisect_rec(n);
			}
			bisect(r);
		}

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
