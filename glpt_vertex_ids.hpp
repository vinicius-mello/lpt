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
#include <map>

// Gaifullin's 15 original points (gaifullin_points() in
// riemann_cp2.cpp) -- every root cell's vertex ids, read directly off
// glpt_gaifullin_cells[seed], already lie in this range.
static const int GLPT_BASE_VERTEX_COUNT = 15;

// Maps a bisected edge (sorted pair of already-known vertex ids) to the
// single canonical id of its midpoint, minting a fresh one (starting
// just past the 15 base points) the first time any cell bisects that
// edge. std::map for correctness/simplicity first -- the open-
// addressing design discussed for the cell-existence table applies
// here too if this ever needs to be faster, but a bisection only
// touches this once per cell (not a hot inner loop), and map's O(log n)
// isn't the bottleneck this design is aimed at (the O(1)-per-leaf
// memory footprint, not lookup speed).
class glpt_edge_cache {
	public:
		glpt_edge_cache() : next_id_(GLPT_BASE_VERTEX_COUNT) {}

		int get_or_create_midpoint(int a, int b) {
			if(a>b) { int t=a; a=b; b=t; }
			uint64_t key = (uint64_t(uint32_t(a))<<32) | uint32_t(b);
			std::map<uint64_t,int>::iterator it = table_.find(key);
			if(it!=table_.end()) return it->second;
			int id = next_id_++;
			table_[key] = id;
			return id;
		}

		size_t size() const { return table_.size(); }
		int next_id() const { return next_id_; }

		// Removal is a plain point-delete, not a sweep: the KEY of the
		// entry that becomes obsolete when a cell bisects is exactly the
		// (sorted parent-id pair) of the edge it just consumed -- known
		// directly at that moment (user's suggestion, 2026-09-22).
		// Provided for that future use; not called by anything here yet
		// (see this file's own header comment on why growth is already
		// bounded without it).
		void forget_edge(int a, int b) {
			if(a>b) { int t=a; a=b; b=t; }
			uint64_t key = (uint64_t(uint32_t(a))<<32) | uint32_t(b);
			table_.erase(key);
		}

	private:
		std::map<uint64_t,int> table_;
		int next_id_;
};

// Root cell's vertex ids: directly glpt_gaifullin_cells[seed][k] --
// already global ids (0..14) by definition, no cache lookup needed.
inline void glpt_root_vertex_ids(int seed, int ids[glpt::DIM+1]) {
	for(int k=0;k<=glpt::DIM;++k) ids[k] = glpt_gaifullin_cells[seed][k];
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
