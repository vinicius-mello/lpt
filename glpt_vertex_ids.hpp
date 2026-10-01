#ifndef GLPT_VERTEX_IDS_HPP
#define GLPT_VERTEX_IDS_HPP

/*! \file
 * \brief Stable, purely-combinatorial global vertex IDs for glpt cells
 * (user's design, 2026-09-22), built on top of glpt.hpp.
 *
 * ============================================================
 * WHY THIS EXISTS
 * ============================================================
 * glpt's own code (level, sigperm, orthant, seed_index) identifies a
 * CELL, and vertex_weights() gives that cell's vertices as barycentric
 * coordinates -- but two DIFFERENT cells that share a vertex (an edge,
 * a face, ...) have NO shared label for it: each only knows its own
 * local vertex index (0..DIM). Any data that needs to be attached to a
 * shared sub-simplex (the motivating case: a root-crossing point cached
 * on a 2-simplex, so that every 4-cell touching that face agrees on
 * exactly where it is, rather than each recomputing it independently
 * with its own Newton seed / rounding and risking a crack in the
 * extracted surface) needs a canonical key for that sub-simplex that
 * does NOT depend on which cell or local order you reached it from.
 *
 * lpt_inc.c's own vertex table (~/code/lua/lpt/) solves an analogous
 * problem by hashing each vertex's REAL coordinates (a Morton code) --
 * workable there, but it means two cells reaching what is geometrically
 * the same point by different arithmetic paths must land on BIT-
 * IDENTICAL floating point to be recognized as the same vertex. glpt
 * has no such risk to begin with: Maubach bisection is exact in binary
 * (see glpt.hpp's vertex_weights_exact()), so instead of hashing
 * coordinates, this file hashes PARENTAGE: a vertex is either one of
 * Gaifullin's 15 original points (fixed IDs 0..14 -- see
 * examples/top/riemann_cp2/riemann_cp2.cpp's gaifullin_points()) or the
 * midpoint of exactly two already-identified vertices, and a small
 * cache (glpt_edge_cache) maps a sorted PAIR of parent IDs to the
 * (single, canonical) ID of their midpoint -- minted once, reused by
 * every cell that later bisects that same edge, however it got there.
 * This is a purely combinatorial identity, never comparing coordinates.
 *
 * ============================================================
 * MEMORY: WHY THIS DOESN'T REINTRODUCE THE nmt<4> BLOWUP
 * ============================================================
 * Every bisection consumes exactly one existing leaf cell and produces
 * exactly two new ones (net +1 leaf) while minting AT MOST one new
 * vertex ID (net +1 edge_cache entry, or 0 if that edge was already
 * bisected from the other side and the cache already had it). So
 * glpt_edge_cache's size stays proportional to the CURRENT mesh's leaf
 * count, not to the total number of bisections ever performed --
 * unlike nmt<4>, which keeps every historical cell around forever.
 *
 * ============================================================
 * HOW A CELL'S VERTEX IDS ARE COMPUTED
 * ============================================================
 * A root's are read directly off glpt_gaifullin_cells[seed] (already
 * IDs 0..14 by definition). A child's are computed from its PARENT's
 * (already known) ids by matching each of the child's DIM+1 vertices,
 * via EXACT-INTEGER vertex_weights_exact() (brought to a common
 * denominator when the child starts a new round, i.e. increases
 * orthant_level), against the parent's DIM+1 vertices: exactly DIM of
 * them match exactly (inherit the parent's id at that position), and
 * exactly one does not -- the new midpoint, whose id comes from
 * glpt_edge_cache keyed on the parent's OWN ids at local positions
 * level() and DIM (the two endpoints Theorem 1 bisects). Matching by
 * exact weight equality, rather than hand-deriving which permutation of
 * parent positions each child inherits, trades a little (utterly
 * negligible: DIM+1=5 vertices, O(DIM^2)=25 integer-vector comparisons
 * per bisection) runtime cost for confidence: this session hit several
 * subtle bugs in hand-derived local-vertex-order/permutation reasoning
 * for this port, so this deliberately reuses the already-validated
 * vertex_weights_exact() as the single source of truth instead of
 * re-deriving the correspondence combinatorially.
 *
 * There is no way to compute a DEEP cell's vertex ids from its bare
 * glpt code alone: id assignment is inherently STATEFUL (which id a
 * given edge's midpoint gets depends on whether some other cell already
 * bisected it), so ids must be threaded through the actual traversal
 * (root -> child -> child -> ...) alongside a single shared
 * glpt_edge_cache, not recomputed independently per cell.
 *
 * ============================================================
 * THE CRITICAL CORRECTNESS QUESTION, AND HOW IT'S CHECKED
 * ============================================================
 * The whole point is that two DIFFERENT cells sharing a facet, having
 * reached their respective vertex ids via generally DIFFERENT
 * bisection histories (different combinatorial paths from possibly
 * different seeds), still agree on the ids of their SHARED vertices.
 * This is not obvious a priori -- it holds only if, whenever two cells
 * are about to bisect what is geometrically the same edge, they always
 * do so with the same (unordered) pair of already-agreed-upon parent
 * ids, so glpt_edge_cache correctly unifies them. glpt_vertex_ids_test.cpp
 * checks this DIRECTLY: build a mesh via many random bisections sharing
 * ONE glpt_edge_cache, then for every neighbor() pair found, verify the
 * shared facet's id SETS match on both sides.
 */

#include "glpt.hpp"
#include "glpt_hash_util.hpp"
#include <cstdlib>

// Gaifullin's 15 original points (gaifullin_points() in
// riemann_cp2.cpp) -- every root cell's vertex ids, read directly off
// glpt_gaifullin_cells[seed], already lie in this range.
static const int GLPT_BASE_VERTEX_COUNT = 15; // Gaifullin; the ACTIVE seed's count is glpt_seed_vertex_count()

// Maps a bisected edge (sorted pair of already-known vertex ids) to the
// single canonical id of its midpoint, minting a fresh one (starting
// just past the 15 base points) the first time any cell bisects that
// edge.
//
// Open addressing (glpt_tree.hpp's own design, reusing its hash/prime
// helpers via glpt_hash_util.hpp), NOT std::map -- changed 2026-09-23
// after directly measuring peak RSS of riemann_cp2_glpt.cpp against
// riemann_cp2.cpp: at ~470000 leaves this table alone had ~470000
// std::map nodes, each paying red-black-tree overhead (3 pointers +
// color, ~32 bytes) on top of the 12 real payload bytes (uint64_t key +
// int value) -- enough on its own to erase glpt_tree's whole per-leaf
// memory saving over riemann_cp2.cpp's append-only nmt<4>. A slot here
// is 16 bytes (key+value+bool, rounded up for alignment) with no tree
// overhead at all.
//
// REVERTED 2026-09-24 from a briefer 20-bit-id-packed single-uint64_t/
// slot design (5cb237d): that packing assumed vertex ids would never
// exceed 2^20-1=1048575, calibrated from a --depth 16 conic measurement
// (71196 vertices) -- a real, user-triggered mesh (elliptic,
// --generic-seed, --depth 20 --threshold 0.02) minted past that and hit
// the packing's own loud assertion rather than silently corrupting
// data. Full 32-bit ids (via a plain uint32_t-range key, not packed
// alongside the value) remove that ceiling entirely -- keys_ below
// packs (a,b) into one uint64_t (32 bits each, exact, no truncation),
// with the minted id and the present flag in their own separate arrays
// instead of squeezed into the same word.
class glpt_edge_cache {
	public:
		explicit glpt_edge_cache(size_t initial_buckets = 257)
			: keys_(0), values_(0), present_(0), nbuckets_(0), count_(0), next_id_(glpt_seed_vertex_count())
		{
			alloc_(glpt_next_prime(initial_buckets));
		}
		~glpt_edge_cache() { std::free(keys_); std::free(values_); std::free(present_); }

		int get_or_create_midpoint(int a, int b) {
			if(a>b) { int t=a; a=b; b=t; }
			uint64_t key = (uint64_t(uint32_t(a))<<32) | uint32_t(b);
			size_t slot = find_slot_(key);
			if(slot!=size_t(-1)) return values_[slot];
			int id = next_id_++;
			insert_(key, id);
			return id;
		}

		size_t size() const { return count_; }
		int next_id() const { return next_id_; }

		// Removal is a plain point-delete, not a sweep: the KEY of the
		// entry that becomes obsolete when a cell bisects is exactly the
		// (sorted parent-id pair) of the edge it just consumed -- known
		// directly at that moment (user's suggestion, 2026-09-22). A
		// no-op if the edge isn't cached (matches std::map::erase's own
		// by-key behavior). Provided for that future use; not called by
		// anything here yet (see this file's own header comment on why
		// growth is already bounded without it). Backward-shift, same
		// derivation and same caveat as glpt_tree.hpp's own erase_() --
		// see that function's comment for the bug an earlier, simpler
		// (stop-at-first-non-mover) version of this exact algorithm had.
		void forget_edge(int a, int b) {
			if(a>b) { int t=a; a=b; b=t; }
			uint64_t key = (uint64_t(uint32_t(a))<<32) | uint32_t(b);
			size_t gap = find_slot_(key);
			if(gap==size_t(-1)) return;
			size_t scan=gap;
			while(true) {
				scan=(scan+1)%nbuckets_;
				if(!present_[scan]) break;
				size_t ideal = glpt_hash64(keys_[scan]) % nbuckets_;
				bool in_range = (gap<=scan) ? (ideal>gap && ideal<=scan) : (ideal>gap || ideal<=scan);
				if(in_range) continue;
				keys_[gap]=keys_[scan]; values_[gap]=values_[scan]; present_[gap]=true;
				gap=scan;
			}
			present_[gap]=false;
			--count_;
		}

	private:
		uint64_t* keys_;
		int* values_;
		bool* present_;
		size_t nbuckets_, count_;
		int next_id_;

		void alloc_(size_t n) {
			keys_ = (uint64_t*)std::calloc(n,sizeof(uint64_t));
			values_ = (int*)std::calloc(n,sizeof(int));
			present_ = (bool*)std::calloc(n,sizeof(bool));
			assert(keys_!=0 && values_!=0 && present_!=0 && "glpt_edge_cache: out of memory");
			nbuckets_ = n;
		}
		size_t find_slot_(uint64_t key) const {
			size_t h = glpt_hash64(key) % nbuckets_;
			while(present_[h]) {
				if(keys_[h]==key) return h;
				h=(h+1)%nbuckets_;
			}
			return size_t(-1);
		}
		void grow_if_needed_() {
			if(double(count_+1) <= 0.7*double(nbuckets_)) return;
			uint64_t* old_k=keys_; int* old_v=values_; bool* old_p=present_;
			size_t old_n=nbuckets_;
			alloc_(glpt_next_prime(2*old_n));
			count_=0;
			for(size_t i=0;i<old_n;++i) if(old_p[i]) insert_(old_k[i], old_v[i]);
			std::free(old_k); std::free(old_v); std::free(old_p);
		}
		void insert_(uint64_t key, int value) {
			grow_if_needed_();
			size_t h = glpt_hash64(key) % nbuckets_;
			while(present_[h]) h=(h+1)%nbuckets_;
			keys_[h]=key; values_[h]=value; present_[h]=true;
			++count_;
		}

		glpt_edge_cache(const glpt_edge_cache&);
		glpt_edge_cache& operator=(const glpt_edge_cache&);
};

// Root cell's vertex ids: directly glpt_gaifullin_cells[seed][k] --
// already global ids (0..14) by definition, no cache lookup needed.
inline void glpt_root_vertex_ids(int seed, int ids[glpt::DIM+1]) {
	for(int k=0;k<=glpt::DIM;++k) ids[k] = glpt_seed_cells[seed][k];
}

// Child (zo=0 or 1) vertex ids, given the PARENT's (already-known) ids.
// See this file's own header comment for the matching-by-exact-weight
// approach and why it was chosen over hand-deriving the permutation.
inline void glpt_child_vertex_ids(const glpt& parent, const int parent_ids[glpt::DIM+1],
		int zo, glpt_edge_cache& cache, int child_ids[glpt::DIM+1]) {
	const int DIM = glpt::DIM;
	glpt child = parent.child(zo);

	int pw[DIM+1][DIM+1], cw[DIM+1][DIM+1];
	int p_shift, c_shift;
	parent.vertex_weights_exact(pw, p_shift);
	child.vertex_weights_exact(cw, c_shift);
	// Bring parent's weights to the child's denominator (child_shift is
	// always p_shift or p_shift+1 -- orthant_level increases by at most
	// 1 per bisection, exactly when a round completes).
	int scale = 1<<(c_shift-p_shift);

	int n_unmatched = 0;
	for(int k=0;k<=DIM;++k) {
		int match_j = -1;
		for(int j=0;j<=DIM;++j) {
			bool eq = true;
			for(int c=0;c<=DIM && eq;++c) if(cw[k][c] != pw[j][c]*scale) eq=false;
			if(eq) { match_j=j; break; }
		}
		if(match_j>=0) {
			child_ids[k] = parent_ids[match_j];
		} else {
			++n_unmatched;
			child_ids[k] = cache.get_or_create_midpoint(parent_ids[parent.level()], parent_ids[DIM]);
		}
	}
	assert(n_unmatched==1 && "child_vertex_ids: expected exactly one new (non-inherited) vertex per bisection");
	(void)n_unmatched;
}

#endif // GLPT_VERTEX_IDS_HPP
