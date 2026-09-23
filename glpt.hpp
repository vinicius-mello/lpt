#ifndef GLPT_HPP
#define GLPT_HPP

/*! \file
 * \brief Gaifullin-seeded LPT code: Atalay & Mount's pointerless
 * hierarchical-simplex location code ("Pointerless Implementation of
 * Hierarchical Simplicial Meshes and Efficient Neighbor Finding in
 * Arbitrary Dimensions", IJCGA 2007), extended from its original
 * domain (the unit hypercube, tiled by Kuhn's d! canonical simplices)
 * to CP^2's real-dimension-4 Gaifullin triangulation (108 four-
 * simplices, from vgtl2.0's examples/top/riemann_cp2/riemann_cp2.cpp).
 *
 * ============================================================
 * WHY THIS FILE EXISTS
 * ============================================================
 * riemann_cp2.cpp triangulates CP^2 adaptively via Maubach bisection
 * on top of vgtl's `nmt<4>` container, which is append-only: a
 * bisected cell is marked not-current but never freed, since other
 * cells' incidence data can reference it. Memory grows monotonically
 * with refinement depth, and this is the binding constraint on how
 * deep/fine a mesh can be built on ordinary hardware (measured
 * directly: ~4-5GB peak RSS around depth 18 on a 7.2GB machine, see
 * examples/top/riemann_cp2/README.md's Memory section).
 *
 * LPT sidesteps this entirely: a cell is identified by a small,
 * self-contained integer CODE (no pointers, no incidence structure
 * stored anywhere), and BOTH the cell's own geometry (Theorem 1 of
 * the paper, not implemented in this file -- see below) AND its
 * same-depth face-neighbors' codes are computable from the code
 * alone, via table lookups, in O(1)-ish time. A leaf cell's storage
 * cost becomes O(1) (one 64-bit integer, optionally with a couple of
 * flag bits for a hash-table-based tree -- see lpt_inc.c in
 * ~/code/lua/lpt/ for that style), and deleting an internal
 * (non-current) cell from an open-addressing hash table is a
 * routine, well-understood operation -- unlike nmt, nothing else
 * ever stores a pointer INTO an LPT-coded cell, so there is no
 * analogue of the incremental-deletion attempt that was tried and
 * abandoned on nmt (two false assumptions about Maubach conformity
 * invariants mid-refinement; see project memory
 * riemann_pc2_bernstein_refinement).
 *
 * ============================================================
 * WHAT CHANGES GOING FROM THE HYPERCUBE TO GAIFULLIN'S CP^2 SEED
 * ============================================================
 * The paper's own d! root simplices (Kuhn's triangulation of
 * [-1,1]^d) are not d! independent geometric objects: they are
 * literally the orbit of Sym(d) acting on ONE reference simplex S_0
 * (Section 3.1-3.2 of the paper). This is why their own neighbor
 * rules (Theorem 2) resolve MOST root-level neighbors (all i with
 * 0<i<d, via Gamma_SWP,i) by a pure permutation-table lookup that
 * never leaves the shared sigperm space: swapping two entries of a
 * Sym(d) permutation always yields another valid Sym(d) permutation,
 * i.e. another one of the SAME d! roots, by closure. Only i=0 and
 * i=d are genuine domain boundaries for their single hypercube (S^(0)
 * and S^(d) are stated to be empty "for all root simplices", with the
 * explicit reason "these faces are faces of the reference hypercube").
 *
 * Gaifullin's 108 cells (examples/top/riemann_cp2/riemann_cp2.cpp's
 * gaifullin_cells_raw) are NOT a group orbit of one reference simplex
 * -- they are 108 independently-specified vertex-order lists, related
 * to each other only by which facets happen to coincide (a fact
 * established by an offline backtracking 5-edge-coloring search, not
 * by any algebraic symmetry). There is no analogue of "swap stays
 * in-family": nothing guarantees that composing a root's permutation
 * code with Gamma_SWP,i lands on another Gaifullin cell's own
 * canonical vertex order. So EVERY one of a Gaifullin root's 5 facets
 * (not just 2 of the hypercube's 5, for d=4) needs external data to
 * resolve -- confirmed computationally (see gaifullin_facet_check.py,
 * run 2026-09-22): all 108*5/2=270 facet-adjacent cell pairs were
 * found, with ZERO exceptions, and critically the shared facet's 4
 * vertices are in IDENTICAL local order on both sides in every case
 * (not just "the same vertex set" -- the same ORDER), which is what
 * makes a permutation-free splice possible at all.
 *
 * This file adds a SEED_INDEX field (0..107) to the LPT code,
 * selecting which Gaifullin cell's descendant tree a given code
 * belongs to. child(), parent(), is_child0(), and the (level,sigperm,
 * orthant) combinatorics are dimension-generic and REUSED VERBATIM
 * from lpt.hpp's lpt_table<4>, since none of it references seed
 * identity.
 *
 * ============================================================
 * NEIGHBOR FINDING: WHY THE PAPER'S OWN Gamma_SWP FORMULA (Theorem 2)
 * IS NOT USED HERE, AND WHAT REPLACES IT
 * ============================================================
 * The paper's Section 6.2 states that same-orthant neighbors (face
 * index 1<=i<=dim, i.e. everything but the "sibling" direction) never
 * leave the current orthant box, for ANY root simplex -- but this is a
 * COROLLARY of Kuhn's d! canonical simplices being the closed orbit of
 * Sym(d) acting on one reference simplex (Section 3.1-3.2): composing
 * a root's signed permutation with Gamma_SWP,i always lands on another
 * member of that SAME orbit, by group closure, so the result can never
 * be outside the orthant that orbit tiles. Gaifullin's 108 cells are
 * NOT such an orbit -- they are 108 independently-specified vertex-
 * order lists related only by which facets happen to coincide (found
 * by an offline backtracking search, not an algebraic symmetry) -- so
 * there is no reason for that closure property to hold here, and
 * empirically it does not: a same-orthant (1<=i<=dim) crossing changes
 * seed_index in roughly a quarter of cases at realistic depths (see
 * glpt_test.cpp test 3's n_seed_changed count), which was independently
 * confirmed by hand for a specific case (seed=99, depth=7, i=3: the
 * cell's facet opposite vertex 3 has, exactly, zero weight at one of
 * the SEED's own 5 vertex positions once vertex 3 is excluded -- a
 * bona fide boundary touch, not a numerical artifact).
 *
 * A first port of Theorem 2's Gamma_SWP/Gamma_RGT/Gamma_LFT/Gamma_NSW
 * formulas (translated field-by-field from lpt.hpp, on the assumption
 * that same-seed neighbors would behave exactly like the hypercube's
 * same-orthant ones) was tried and found to be genuinely buggy: cross-
 * checked against the barycentric-convexity invariant that every cell
 * in a Maubach-bisection tree must satisfy (every vertex is a convex
 * combination of the SEED's own original DIM+1 vertices, since the
 * tree is built entirely by repeated midpointing), it produced a
 * weight outside [0,1] -- a geometric impossibility -- in 610/2651
 * (23.0%) of sampled same-orthant cases (2026-09-22). That rate lines
 * up with the ~24% of same-orthant cases that legitimately DO cross
 * seeds (above): the formula was silently computing garbage for the
 * seed-crossing cases it has no way to detect, since it was ported
 * assuming (per the hypercube's own guarantee) that this could never
 * happen.
 *
 * glpt::neighbor() therefore does not use Theorem 2's formula at all.
 * The SIBLING-direction facet (i==DIM for a 0-child, i==lstar for a
 * 1-child -- the one introduced by a cell's own most recent bisection)
 * is resolved as a plain combinatorial fact, parent().child(other bit),
 * with no floating point at all: by the definition of Maubach bisection
 * this is ALWAYS exactly the correct neighbor, same seed, no ambiguity.
 * (This was not the original design -- routing the sibling facet
 * through the geometric path below broke ~5-15% of sampled cases,
 * because that facet sits at a finer scale than the path's nudge margin
 * is tuned for; found via glpt_test.cpp's test 4 and fixed 2026-09-22.)
 *
 * For every OTHER (non-sibling) i, the O(1)/O(d) Theorem-2 formula IS
 * used after all -- but gated by an EXACT (integer, no floating point)
 * test rather than trusted outright (user's suggestion, 2026-09-22):
 * root_boundary_facet() determines, via the SAME Theorem-1-derived
 * barycentric weights as vertex_weights() but computed as exact dyadic-
 * rational integers (repeated bisection is exact in binary; denominator
 * never exceeds 2^(MAX_ORTHANT_LEVEL+1)=2048), whether facet i lies
 * EXACTLY on one of the seed's own 5 bounding hyperplanes -- the actual
 * geometric fact that determines whether a crossing happens, checked
 * directly rather than inferred from the candidate's own self-
 * consistency (which an earlier attempt, fast_candidate()+traces_to_
 * own_root(), relied on and found unsafe: it can verify a genuine-but-
 * WRONG cell for Gaifullin, unlike Kuhn's hypercube where the group
 * structure rules that out -- see neighbor()'s own comment). Gated this
 * way, fast_candidate() (Theorem 2's formula, called only once
 * root_boundary_facet() has confirmed the facet is interior) matched
 * the independently-verified glpt_neighbor_geometric() in 132799/132801
 * (99.9985%) of sampled cases -- the earlier "23% convexity-violation"
 * finding turns out to have been entirely about the formula's inability
 * to DETECT which facets cross, not an error in the formula itself.
 *
 * glpt_neighbor_geometric() -- reconstruct real (barycentric)
 * coordinates via Theorem 1 + Lemma 6 (vertex_weights()), nudge a query
 * point across the requested facet, and relocate it via the paper's own
 * search() procedure (Fig. 9 / Lemma 7, glpt_locate()), redirecting
 * through the Gaifullin adjacency table when the crossed facet lies on
 * a seed boundary -- is therefore only reached for root cells (every
 * facet is a boundary by definition) and confirmed boundary facets,
 * i.e. only when a seed change is actually needed. Both its ingredients
 * were independently validated (vertex_weights() against a from-scratch
 * Maubach-bisection simulation; glpt_locate() via self-location of a
 * cell's own centroid, 2000/2000). Getting its nudge margin right at
 * extreme depth took a second pass (see its own eps comment): too tight
 * and floating-point rounding silently swallows it past simplex_level
 * ~35 (neighbor() returning the source cell as its own neighbor); too
 * loose and it overshoots into the wrong cell -- a narrow residual
 * failure at the exact maximum representable depth (43) was later found
 * to disappear entirely once the exact gate above routes interior
 * facets away from this function altogether, since that was the only
 * place still asking it to do floating-point work at that depth.
 *
 * Root-level (orthant_level==0, i.e. is_root()) crossings always go
 * through the Gaifullin adjacency table, using the verified identical-
 * facet-order fact (gaifullin_facet_check.py, 2026-09-22: all 108*5/2=
 * 270 facet-adjacent cell pairs found, shared facets in IDENTICAL local
 * order on both sides) -- the result is simply the neighbor's own fresh
 * root (sigperm=0, level=0, orthant_level=0), no permutation correction
 * needed. This is subsumed by glpt_neighbor_geometric() too, since a
 * root cell IS its own seed's whole body and every one of its 5 facets
 * lies on a seed boundary by definition.
 *
 * See glpt_test.cpp for what is verified: child()/parent() round-trip
 * consistency; an involution check (cross via facet i, then back via
 * the SAME i -- valid for any NON-sibling i, since local facet index
 * always matches on both sides whether or not the seed changes, but NOT
 * valid for the sibling direction, whose return crossing is generally a
 * DIFFERENT index -- see test 3's and test 7's own comments) at root
 * level, same-orthant depth, and deep/arbitrary depth up to and
 * including the representable maximum; the sibling relationship against
 * parent()/child(); the convexity invariant; test 7's full sweep of all
 * of the above together, at every facet index, across the complete
 * representable depth range in one pass; and (test 8) the exact
 * boundary gate itself cross-checked directly against the geometric
 * ground truth.
 */

#include <cstdint>
#include <cassert>
#include "lpt.hpp" // reuses lpt_table<4>: dimension-generic (level,
                    // sigperm) combinatorics, independent of seed.

// --- Gaifullin's 108-cell CP^2 seed (copied verbatim from
// examples/top/riemann_cp2/riemann_cp2.cpp's gaifullin_cells_raw,
// already in Maubach-compatible per-cell vertex order -- local
// position i = the color a 5-edge-coloring search assigned to the
// facet opposite that vertex) -- plus the facet-adjacency table
// derived from it below.
static const int GLPT_NCELLS = 108;
static const int glpt_gaifullin_cells[GLPT_NCELLS][5] = {
	{0,3,7,9,13}, {0,3,7,9,14}, {0,3,7,10,12}, {0,3,7,10,14}, {0,3,7,11,12}, {0,3,7,11,13},
	{0,3,8,9,13}, {0,3,8,9,14}, {0,3,8,10,12}, {0,3,8,10,14}, {0,3,8,11,12}, {0,3,8,11,13},
	{0,4,6,9,13}, {0,4,6,9,14}, {0,4,6,10,12}, {0,4,6,10,14}, {0,4,6,11,12}, {0,4,6,11,13},
	{0,4,8,9,13}, {0,4,8,9,14}, {0,4,8,10,12}, {0,4,8,10,14}, {0,4,8,11,12}, {0,4,8,11,13},
	{0,5,6,9,13}, {0,5,6,9,14}, {0,5,6,10,12}, {0,5,6,10,14}, {0,5,6,11,12}, {0,5,6,11,13},
	{0,5,7,9,13}, {0,5,7,9,14}, {0,5,7,10,12}, {0,5,7,10,14}, {0,5,7,11,12}, {0,5,7,11,13},
	{1,3,6,10,13}, {1,3,6,10,14}, {1,3,6,11,13}, {1,3,6,11,14}, {1,3,7,10,12}, {1,3,7,10,14},
	{1,3,7,11,12}, {1,3,7,11,14}, {1,3,8,10,12}, {1,3,8,10,13}, {1,3,8,11,12}, {1,3,8,11,13},
	{1,4,6,9,13}, {1,4,6,9,14}, {1,4,6,11,13}, {1,4,6,11,14}, {1,4,7,9,12}, {1,4,7,9,14},
	{1,4,7,11,12}, {1,4,7,11,14}, {1,4,8,9,12}, {1,4,8,9,13}, {1,4,8,11,12}, {1,4,8,11,13},
	{1,5,6,9,13}, {1,5,6,9,14}, {1,5,6,10,13}, {1,5,6,10,14}, {1,5,7,9,12}, {1,5,7,9,14},
	{1,5,7,10,12}, {1,5,7,10,14}, {1,5,8,9,12}, {1,5,8,9,13}, {1,5,8,10,12}, {1,5,8,10,13},
	{2,3,6,10,13}, {2,3,6,10,14}, {2,3,6,11,13}, {2,3,6,11,14}, {2,3,7,9,13}, {2,3,7,9,14},
	{2,3,7,11,13}, {2,3,7,11,14}, {2,3,8,9,13}, {2,3,8,9,14}, {2,3,8,10,13}, {2,3,8,10,14},
	{2,4,6,10,12}, {2,4,6,10,14}, {2,4,6,11,12}, {2,4,6,11,14}, {2,4,7,9,12}, {2,4,7,9,14},
	{2,4,7,11,12}, {2,4,7,11,14}, {2,4,8,9,12}, {2,4,8,9,14}, {2,4,8,10,12}, {2,4,8,10,14},
	{2,5,6,10,12}, {2,5,6,10,13}, {2,5,6,11,12}, {2,5,6,11,13}, {2,5,7,9,12}, {2,5,7,9,13},
	{2,5,7,11,12}, {2,5,7,11,13}, {2,5,8,9,12}, {2,5,8,9,13}, {2,5,8,10,12}, {2,5,8,10,13},
};

// glpt_gaifullin_neighbor[c][i] = the cell sharing cell c's facet
// opposite its i-th local vertex (or -1 if none was found -- should
// never happen for a closed manifold's own triangulation, and is
// asserted against in glpt_gaifullin_init()). Filled lazily by
// glpt_gaifullin_init(), which the application must call once before
// using any glpt root-level neighbor() call.
static int glpt_gaifullin_neighbor[GLPT_NCELLS][5];
static bool glpt_gaifullin_ready = false;

inline void glpt_gaifullin_init() {
	if(glpt_gaifullin_ready) return;
	for(int c=0;c<GLPT_NCELLS;++c)
		for(int i=0;i<5;++i)
			glpt_gaifullin_neighbor[c][i] = -1;
	for(int c=0;c<GLPT_NCELLS;++c) {
		for(int i=0;i<5;++i) {
			// the facet of cell c opposite local vertex i: its other 4 vertices
			int facet[4]; int k=0;
			for(int j=0;j<5;++j) if(j!=i) facet[k++]=glpt_gaifullin_cells[c][j];
			if(glpt_gaifullin_neighbor[c][i]!=-1) continue; // already found from the other side
			for(int c2=c+1;c2<GLPT_NCELLS;++c2) {
				// count shared vertices between {facet} and cell c2's 5 vertices
				int shared=0, other_i=-1;
				for(int jj=0;jj<5;++jj) {
					bool in_facet=false;
					for(int kk=0;kk<4;++kk) if(glpt_gaifullin_cells[c2][jj]==facet[kk]) { in_facet=true; break; }
					if(in_facet) ++shared; else other_i=jj;
				}
				if(shared==4) {
					glpt_gaifullin_neighbor[c][i]=c2;
					glpt_gaifullin_neighbor[c2][other_i]=c;
					break;
				}
			}
			assert(glpt_gaifullin_neighbor[c][i]>=0
				&& "Gaifullin triangulation has no boundary facets (CP^2 is closed) -- "
				   "a -1 here means glpt_gaifullin_cells doesn't match riemann_cp2.cpp's "
				   "gaifullin_cells_raw, or that data itself changed.");
		}
	}
	glpt_gaifullin_ready = true;
}

// --- The extended LPT code itself: dim=4 fixed (CP^2's real
// dimension), (level, sigperm, orthant_level, orthant[4]) exactly as
// lpt4_48_bits, plus a 7-bit seed_index (0..107) selecting which
// Gaifullin cell. Packed explicitly into one uint64_t (62 of 64 bits
// used) via shift/mask rather than C++ bitfields, whose layout order
// isn't guaranteed by the standard -- lpt.hpp's own lpt4_48_bits uses
// bitfields and works fine on the toolchains it's been built with,
// but this file has no such history yet, so it takes the more
// portable route.
class glpt;
// Forward declaration -- defined after the class, needs its full
// interface (vertex_weights/simplex_level/seed_index/child), but is
// called from within neighbor(), a member function defined inline
// inside the class body itself. See the GEOMETRIC FALLBACK comment
// (above glpt::vertex_weights, inside the class) for what this does.
inline bool glpt_neighbor_geometric(const glpt& c, int local_excluded_i, glpt& result);

class glpt {
	public:
		static const int DIM = 4;
	private:
		static const int SEED_BITS=7, LEVEL_BITS=2, SIGPERM_BITS=9, ORTHLVL_BITS=4, ORTH_BITS=10;
		static const int SEED_OFF=0;
		static const int LEVEL_OFF=SEED_OFF+SEED_BITS;       // 7
		static const int SIGPERM_OFF=LEVEL_OFF+LEVEL_BITS;   // 9
		static const int ORTHLVL_OFF=SIGPERM_OFF+SIGPERM_BITS; // 18
		static const int ORTH_OFF=ORTHLVL_OFF+ORTHLVL_BITS;  // 22 (+10 per coordinate)
		uint64_t code_;

		static uint64_t mask(int bits) { return (bits>=64) ? ~uint64_t(0) : ((uint64_t(1)<<bits)-1); }
		int get(int off,int bits) const { return int((code_>>off) & mask(bits)); }
		void set(int off,int bits,int v) {
			code_ &= ~(mask(bits)<<off);
			code_ |= (uint64_t(v) & mask(bits)) << off;
		}
	public:
		// Each subdivision ROUND appends exactly one bit to every orthant
		// coordinate, so ORTH_BITS (the coordinate width, not 2^ORTH_BITS-1)
		// is the true cap on representable rounds; ORTHLVL_BITS must also be
		// able to STORE that count. min() of the two, matching lpt.hpp's
		// own convention (e.g. lpt4_48_bits: orthant_level_ 4 bits stores
		// up to max_orthant_level=8, the width of its own orthant_[4]).
		static const int MAX_ORTHANT_LEVEL = ORTH_BITS<((1<<ORTHLVL_BITS)-1) ? ORTH_BITS : ((1<<ORTHLVL_BITS)-1); // = 10

		glpt() : code_(0) {}
		//! Constructs the root of Gaifullin cell `seed` (0..107): level=0,
		//! sigperm=0 (identity -- the cell's own vertex list IS already the
		//! canonical Maubach order), orthant_level=0.
		explicit glpt(int seed) : code_(0) { set(SEED_OFF,SEED_BITS,seed); }

		int seed_index() const { return get(SEED_OFF,SEED_BITS); }
		int level() const { return get(LEVEL_OFF,LEVEL_BITS); }
		int sigperm_get() const { return get(SIGPERM_OFF,SIGPERM_BITS); }
		void sigperm_set(int v) { set(SIGPERM_OFF,SIGPERM_BITS,v); }
		int orthant_level() const { return get(ORTHLVL_OFF,ORTHLVL_BITS); }
		void orthant_level_set(int v) { set(ORTHLVL_OFF,ORTHLVL_BITS,v); }
		int orthant_get(int i) const { return get(ORTH_OFF+i*ORTH_BITS,ORTH_BITS); }
		void orthant_set(int i,int v) { set(ORTH_OFF+i*ORTH_BITS,ORTH_BITS,v); }
		int simplex_level() const { return orthant_level()*DIM+level(); }
		bool is_root() const { return level()==0 && orthant_level()==0; }

		uint64_t raw() const { return code_; }
		bool operator==(const glpt& r) const { return code_==r.code_; }
		bool operator<(const glpt& r) const { return code_<r.code_; }

		// --- dimension-generic table lookups, verbatim from lpt.hpp's
		// lpt<4,code>, operating on lpt_table<4> (curve/seed-independent).
		int orth() const { return lpt_table<DIM>::orth[sigperm_get()]; }
		int sigma_inv(int l) const { return lpt_table<DIM>::sigma_inv[sigperm_get()][l]; }
		int sigma(int l) const { return lpt_table<DIM>::sigma[sigperm_get()][l]; }
		int lft(int l) const { return lpt_table<DIM>::lft[sigperm_get()][l]; }
		int rgt(int l) const { return lpt_table<DIM>::rgt[sigperm_get()][l]; }
		int sigperm_entry(int i) const { return lpt_table<DIM>::sigperm[sigperm_get()][i-1]; }
		int inverse() const { return lpt_table<DIM>::inverse[sigperm_get()]; }
		int neg() const { return lpt_table<DIM>::neg[sigperm_get()]; }
		int nsw(int l) const { return lpt_table<DIM>::nsw[sigperm_get()][l-1]; }
		int swp(int i) const { return lpt_table<DIM>::swp[sigperm_get()][i-1]; }

		//! Same as lpt<dim,code>::is_child0 -- dimension-generic, no seed
		//! dependence (see this function's own comment in lpt.hpp; the
		//! "orthant_level==0" branch here is never actually consulted by
		//! neighbor(), since root-level crossings are resolved entirely via
		//! the Gaifullin table before this would matter -- kept for
		//! completeness/parent() correctness on cells just above a root).
		bool is_child0() const {
			if(orthant_level()==0) {
				if(level()!=0) { int s=sigperm_entry(level()); return (s>0); }
				return false;
			} else if(orthant_level()==1) {
				if(level()!=0) {
					int lstar=level();
					int s=sigperm_entry(lstar);
					return (s>0)==((orthant_get((s<0?-s:s)-1)&1)==0);
				} else {
					int s=sigperm_entry(DIM);
					return (s>0);
				}
			} else {
				int lstar=(level()==0)?DIM:level();
				int s=sigperm_entry(lstar);
				bool sbit;
				if(level()!=0) sbit=(orthant_get((s<0?-s:s)-1)&1)!=0;
				else sbit=(orthant_get((s<0?-s:s)-1)&2)!=0;
				return (s>0)==!sbit;
			}
		}

		//! Child 0 or 1 (zo). Identical logic/formulas to lpt.hpp's
		//! lpt<4,code>::child -- dimension-generic, carries seed_index
		//! through unchanged (bisection within a seed's own tree never
		//! changes which seed you're in).
		glpt child(int zo) const {
			assert(zo==0||zo==1);
			assert(orthant_level()<MAX_ORTHANT_LEVEL || level()<DIM-1);
			glpt r; r.set(SEED_OFF,SEED_BITS,seed_index());
			r.set(LEVEL_OFF,LEVEL_BITS,(level()+1)%DIM);
			if(zo==0) r.sigperm_set(sigperm_get());
			else r.sigperm_set(sigma(level()));
			if(r.level()==0) {
				r.orthant_level_set(orthant_level()+1);
				int o=r.orth();
				for(int i=0;i<DIM;++i) {
					r.orthant_set(i, orthant_get(i)<<1);
					if(o&(1<<i)) r.orthant_set(i, r.orthant_get(i)|1);
				}
			} else {
				r.orthant_level_set(orthant_level());
				for(int i=0;i<DIM;++i) r.orthant_set(i, orthant_get(i));
			}
			return r;
		}

		//! Parent, within the SAME seed's tree -- undefined for a root
		//! (assert). Never crosses seeds: only neighbor() does that.
		glpt parent() const {
			assert(!is_root());
			glpt r; r.set(SEED_OFF,SEED_BITS,seed_index());
			if(level()==0) {
				r.set(LEVEL_OFF,LEVEL_BITS,DIM-1);
				r.orthant_level_set(orthant_level()-1);
				for(int i=0;i<DIM;++i) r.orthant_set(i, orthant_get(i)>>1);
			} else {
				r.set(LEVEL_OFF,LEVEL_BITS,level()-1);
				r.orthant_level_set(orthant_level());
				for(int i=0;i<DIM;++i) r.orthant_set(i, orthant_get(i));
			}
			if(is_child0()) r.sigperm_set(sigperm_get());
			else r.sigperm_set(sigma_inv(r.level()));
			return r;
		}

		// ============================================================
		// GEOMETRIC FALLBACK for the deep i=0 seed-crossing case (see
		// this file's own header comment for why the combinatorial/O(1)
		// approach was left unimplemented). Reuses Theorem 1 of the paper
		// (the SAME closed-form reconstruction lpt.hpp's own simplex()
		// already implements for the plain hypercube -- ported here
		// verbatim, since it's 100% dimension-generic and doesn't
		// reference seed identity at all) to get this cell's 5 vertices
		// as REAL barycentric weights relative to the SEED's own 5
		// original vertices, then Lemma 7 (the paper's own point-location
		// descent rule, Fig. 9) to walk a query point down a DIFFERENT
		// seed's tree to the matching depth. This sidesteps the
		// unresolved "does the orthant-axis concept even transfer
		// between seeds" question entirely: everything here is plain
		// real-number barycentric arithmetic, never combinatorial
		// permutation bookkeeping across the seed boundary.
		//
		// Correctness of each piece, and how it's checked (glpt_test.cpp):
		//  - vertex_weights() (Theorem 1 + Lemma 6): cross-checked against
		//    the ALREADY-VERIFIED fast Gamma_SWP same-orthant formula, by
		//    computing a same-seed (i=1..DIM-1) neighbor BOTH ways and
		//    comparing.
		//  - the Gaifullin-table redirect + assumption that a facet's
		//    barycentric weights transfer UNCHANGED between the two
		//    incident seeds: relies on the identical-local-vertex-order
		//    fact verified 2026-09-22 (gaifullin_facet_check.py) -- tested
		//    here via an involution check (cross to the neighbor seed,
		//    then cross back, must return to the start).
		//
		// Computes this cell's DIM+1 vertices, each as a (DIM+1)-vector of
		// barycentric weights relative to the SEED's own DIM+1 original
		// vertices (i.e. relative to glpt_gaifullin_cells[seed_index()]'s
		// own vertex order). w[k][j] = weight of ORIGINAL seed vertex j in
		// the k-th vertex of THIS cell.
		void vertex_weights(double w[DIM+1][DIM+1]) const {
			double pt[DIM+1][DIM];
			double t[DIM];
			for(int j=0;j<DIM;++j) t[j]=0.0;
			int ol = orthant_level();
			for(int i=0;i<ol;++i) {
				for(int j=0;j<DIM;++j) {
					if(orthant_get(j)&(1<<(ol-1-i))) t[j]-=1.0/(2<<i);
					else t[j]+=1.0/(2<<i);
				}
			}
			for(int i=0;i<=DIM;++i) {
				for(int j=0;j<DIM;++j) {
					double v = (i>j) ? 1.0 : ((i<level()) ? 0.0 : -1.0);
					v = v/(1<<ol);
					int p = sigperm_entry(j+1);
					int col = (p<0?-p:p)-1;
					pt[i][col] = ((p<0)?-1.0:1.0)*v;
				}
				for(int j=0;j<DIM;++j) pt[i][j] += t[j];
			}
			// Lemma 6 (paper, Eq. 23): convert each R^DIM point into
			// DIM+1 barycentric weights relative to the base reference
			// simplex S_empty -- derived by hand from S_empty's own
			// staircase pattern (S_empty[i][j] = +1 if i>j else -1):
			// writing p_j = sum_i w_i*S_empty[i][j] and using sum(w)=1
			// telescopes to w_0=(1-p_0)/2, w_k=(p_{k-1}-p_k)/2 for
			// 0<k<DIM, w_DIM=(1+p_{DIM-1})/2 -- matches Eq. 23 exactly
			// (checked against the paper's own printed formula).
			for(int i=0;i<=DIM;++i) {
				w[i][0] = (1.0-pt[i][0])/2.0;
				for(int k=1;k<DIM;++k) w[i][k] = (pt[i][k-1]-pt[i][k])/2.0;
				w[i][DIM] = (1.0+pt[i][DIM-1])/2.0;
			}
		}

		//! EXACT (integer, no floating point at all) test: does facet
		//! local_excluded_i of this cell lie exactly on one of the SEED's
		//! own DIM+1 bounding hyperplanes -- i.e. is this a genuine
		//! seed-crossing facet, as opposed to a facet strictly interior
		//! to the seed's body? If so, returns true and sets orig_excluded
		//! to the seed's own vertex position whose opposite hyperplane it
		//! lies on (see vertex_weights()'s own comment for the numbering
		//! convention -- this is the SAME orig_excluded glpt_neighbor_
		//! geometric computes via floating point, but exact here).
		//!
		//! Every intermediate quantity in vertex_weights() (Theorem 1 +
		//! Lemma 6) is dyadic rational with denominator exactly
		//! 2^(orthant_level()+1) -- repeated bisection only ever halves,
		//! adds, or subtracts -- so it is exactly representable as a
		//! small integer (orthant_level()<=MAX_ORTHANT_LEVEL=10, so the
		//! denominator is at most 2048) with NO precision concerns at any
		//! representable depth, unlike vertex_weights()'s own double
		//! version. A facet lies on seed position j's hyperplane iff
		//! EVERY vertex other than local_excluded_i has EXACTLY zero
		//! weight there (weights are non-negative, so a sum of them is
		//! zero iff every term is) -- checked per vertex directly, no
		//! combination step (and so no need for GLPT_IRRATIONAL_W-style
		//! decorrelation) needed for this determination.
		//!
		//! Verified (2026-09-22, user's suggestion): cross-checked
		//! against glpt_neighbor_geometric() across ~160000 sampled
		//! (cell,i) pairs -- every case this function calls a boundary
		//! facet DID cross seeds (0 false positives), and for every case
		//! it calls interior, fast_candidate() (the Theorem-2 formula,
		//! otherwise found unsafe to trust blindly -- see neighbor()'s
		//! own comment) matched glpt_neighbor_geometric() in 132799/
		//! 132801 (99.9985%) of cases, the 2 exceptions both at the exact
		//! maximum representable depth (43) where glpt_neighbor_geometric
		//! itself has a separate, already-documented floating-point
		//! precision edge case (see its own eps comment) -- i.e. this
		//! function's OWN answer is exact and depth-independent, and the
		//! 2 disagreements are plausibly glpt_neighbor_geometric being
		//! wrong there, not fast_candidate(). This confirms the earlier
		//! "23% convexity-violation" finding was entirely about the
		//! formula's inability to DETECT a crossing, not an error in the
		//! formula itself for facets that don't cross.
		bool root_boundary_facet(int local_excluded_i, int& orig_excluded) const {
			int ol = orthant_level();
			int t_num[DIM];
			for(int j=0;j<DIM;++j) t_num[j]=0;
			for(int i=0;i<ol;++i) {
				int term = 1<<(ol-1-i);
				for(int j=0;j<DIM;++j) {
					if(orthant_get(j)&(1<<(ol-1-i))) t_num[j]-=term;
					else t_num[j]+=term;
				}
			}
			int w_num[DIM+1][DIM+1];
			for(int k=0;k<=DIM;++k) {
				if(k==local_excluded_i) continue;
				int pt_num[DIM];
				for(int j=0;j<DIM;++j) {
					int v_raw = (k>j) ? 1 : ((k<level()) ? 0 : -1);
					int p = sigperm_entry(j+1);
					int col = (p<0?-p:p)-1;
					pt_num[col] = ((p<0)?-v_raw:v_raw) + t_num[col];
				}
				w_num[k][0] = (1<<ol) - pt_num[0];
				for(int c=1;c<DIM;++c) w_num[k][c] = pt_num[c-1]-pt_num[c];
				w_num[k][DIM] = (1<<ol) + pt_num[DIM-1];
			}
			for(int j=0;j<=DIM;++j) {
				bool all_zero=true;
				for(int k=0;k<=DIM && all_zero;++k) {
					if(k==local_excluded_i) continue;
					if(w_num[k][j]!=0) all_zero=false;
				}
				if(all_zero) { orig_excluded=j; return true; }
			}
			return false;
		}

		//! EXACT (integer) version of vertex_weights(), for the SAME
		//! reason as root_boundary_facet(): every quantity is dyadic
		//! rational with denominator exactly 2^(orthant_level()+1) (see
		//! that function's own comment), representable exactly as a
		//! small integer at any depth. w_num[k][j] = w[k][j] *
		//! 2^(orthant_level()+1) exactly; denom_shift is set to
		//! orthant_level()+1 (i.e. the actual denominator is
		//! 1<<denom_shift) so callers combining weights from cells at
		//! DIFFERENT orthant_level (e.g. a parent and its child, when the
		//! child starts a new round) can bring them to a common
		//! denominator before comparing. Used by the vertex-ID layer
		//! (glpt_vertex_ids.hpp) to identify, by EXACT equality rather
		//! than floating-point closeness, which of a child's vertices are
		//! inherited unchanged from the parent and which one is the new
		//! bisection midpoint.
		void vertex_weights_exact(int w_num[DIM+1][DIM+1], int& denom_shift) const {
			int ol = orthant_level();
			denom_shift = ol+1;
			int t_num[DIM];
			for(int j=0;j<DIM;++j) t_num[j]=0;
			for(int i=0;i<ol;++i) {
				int term = 1<<(ol-1-i);
				for(int j=0;j<DIM;++j) {
					if(orthant_get(j)&(1<<(ol-1-i))) t_num[j]-=term;
					else t_num[j]+=term;
				}
			}
			for(int k=0;k<=DIM;++k) {
				int pt_num[DIM];
				for(int j=0;j<DIM;++j) {
					int v_raw = (k>j) ? 1 : ((k<level()) ? 0 : -1);
					int p = sigperm_entry(j+1);
					int col = (p<0?-p:p)-1;
					pt_num[col] = ((p<0)?-v_raw:v_raw) + t_num[col];
				}
				w_num[k][0] = (1<<ol) - pt_num[0];
				for(int c=1;c<DIM;++c) w_num[k][c] = pt_num[c-1]-pt_num[c];
				w_num[k][DIM] = (1<<ol) + pt_num[DIM-1];
			}
		}

		// Sentinel returned by neighbor() when the geometric approach
		// itself can't resolve a result (should be essentially never --
		// see glpt_neighbor_geometric's own comment).
		enum result { OK=0, GLPT_DEEP_CROSSING_UNIMPLEMENTED=1 };

		// ============================================================
		// FAST CANDIDATE (Theorem 2's formula), GATED BY THE EXACT
		// BOUNDARY TEST ABOVE -- second design, 2026-09-22
		// ============================================================
		// First attempt at recovering an O(1)/O(d) path (not used
		// below): reintroduce Theorem 2's Gamma_SWP-family formula as a
		// candidate generator, gated by walking the candidate back up to
		// its own root via parent() and checking sigperm==identity there
		// (traces_to_own_root(), since removed). That check only proves
		// the candidate is A genuine cell of seed_index()'s own tree, not
		// THE correct neighbor -- for Kuhn's hypercube those coincide
		// (Sym(DIM) acts simply transitively on the d! canonical
		// simplices), but Gaifullin's 108 cells have no such group
		// structure, so nothing stops the formula's output from tracing
		// back to identity by pure combinatorial coincidence while
		// denoting some OTHER, unrelated cell. Measured directly: about
		// 1-6% of "verified" candidates were the wrong cell (e.g.
		// seed=80, depth=13, i=0: fast_candidate() gives a self-
		// consistent, root-tracing cell still in seed 80, while the true
		// neighbor is in seed 81) -- unsafe, so abandoned as a gate.
		//
		// The WORKING gate (user's suggestion): root_boundary_facet()
		// above determines, via exact integer arithmetic (not a
		// candidate's after-the-fact self-consistency), whether facet i
		// actually lies on one of the SEED's own 5 bounding hyperplanes
		// -- the precise geometric fact that matters, computed directly
		// rather than inferred. Cross-checked against glpt_neighbor_
		// geometric() across ~160000 sampled (cell,i) pairs (2026-09-22):
		// every facet root_boundary_facet() called a boundary genuinely
		// crossed seeds (0 false positives), and for every facet it
		// called interior, fast_candidate() (this function) matched
		// glpt_neighbor_geometric() in 132799/132801 (99.9985%) of cases
		// -- the 2 exceptions both at the exact maximum representable
		// depth (43), where glpt_neighbor_geometric itself has a
		// separate, already-documented floating-point precision edge
		// case (see its own eps comment), so plausibly ITS error, not
		// this function's. This retroactively explains the earlier
		// "23% convexity-violation" finding: the formula was never wrong
		// for facets that don't cross -- it was simply blind to which
		// facets do, exactly the blindness root_boundary_facet() fixes.
		//
		// *this*'s own sigperm_entry(1) (not the candidate's) is used
		// below for the i==0 orthant-axis test, matching lpt.hpp's own
		// neighbor_bd exactly (that axis is determined by *this*'s own
		// orientation -- which coordinate direction "vertex 0" differs
		// along -- not by whichever sigperm the candidate ends up with).
		// Only called by neighbor() once root_boundary_facet() has
		// already confirmed facet i is interior; empirically (the sweep
		// above) the orthant-extremum branch below (i==0, x==m) never
		// fired incorrectly for those cases, consistent with an interior
		// i==0 facet not being at that extremum, though this was
		// observed rather than separately proved.
		glpt fast_candidate(int i) const {
			glpt r; r.set(SEED_OFF,SEED_BITS,seed_index());
			r.set(LEVEL_OFF,LEVEL_BITS,level());
			r.orthant_level_set(orthant_level());
			if(is_child0()) {
				if(i==0) r.sigperm_set(neg());
				else if(i==DIM) {
					int lminus=(level()==0)?(DIM-1):(level()-1);
					r.sigperm_set(lft(lminus));
				} else r.sigperm_set(swp(i));
			} else {
				int lstar=(level()==0)?DIM:level();
				if(i==0) r.sigperm_set(neg());
				else if(i==lstar) {
					int lminus=(level()==0)?(DIM-1):(level()-1);
					r.sigperm_set(rgt(lminus));
				} else if(i==DIM) r.sigperm_set(nsw(level()));
				else r.sigperm_set(swp(i));
			}
			int p = sigperm_entry(1); // *this*'s own -- see comment above
			if(orthant_level()==0) {
				for(int j=0;j<DIM;++j) r.orthant_set(j,0);
				return r;
			}
			if(i==0) {
				for(int j=0;j<DIM;++j) r.orthant_set(j, orthant_get(j));
				int d = LPT_SGN(p);
				int x = r.orthant_get(LPT_ABS(p)-1);
				int m = (d==1) ? 0 : ((1<<orthant_level())-1);
				if(x!=m) {
					int mm=(1<<orthant_level())-1;
					if(d==1) r.orthant_set(LPT_ABS(p)-1, mm&(x-1));
					else r.orthant_set(LPT_ABS(p)-1, mm&(x+1));
				}
				// x==m: left as the plain copy above -- an orthant
				// extremum, expected to fail root-trace verification.
			} else if(level()==0 && i==DIM) {
				int o = r.orth();
				for(int j=0;j<DIM;++j) {
					r.orthant_set(j, ((orthant_get(j)>>1)<<1));
					if(o&(1<<j)) r.orthant_set(j, r.orthant_get(j)|1);
				}
			} else {
				for(int j=0;j<DIM;++j) r.orthant_set(j, orthant_get(j));
			}
			return r;
		}

		//! Computes the face-i neighbor N^(i)(this) (0<=i<=DIM), at the
		//! SAME simplex_level (this cell need not itself be a "leaf" of any
		//! particular tree -- that bookkeeping is a caller/tree concern).
		//!
		//! Three cases, cheapest first:
		//!  1. SIBLING direction (i==DIM for a 0-child, i==lstar for a
		//!     1-child -- the facet introduced by *this* cell's own most
		//!     recent bisection): a plain combinatorial fact, always
		//!     exactly parent().child(the OTHER bit) -- same seed, no
		//!     ambiguity, no floating point, O(1) given parent()/child().
		//!     Found necessary, not just a nice-to-have: this facet sits
		//!     at a finer scale than glpt_neighbor_geometric's nudge-
		//!     epsilon is tuned for, and routing it through that path
		//!     caused measurable failures (up to ~5-15% of sampled cases).
		//!  2. Non-root, non-sibling, and root_boundary_facet() (exact
		//!     integer test, see its own comment) confirms facet i does
		//!     NOT lie on one of the seed's own bounding hyperplanes:
		//!     fast_candidate() (Theorem 2's Gamma_SWP-family formula) is
		//!     used directly, O(1)/O(d), no floating point -- safe here
		//!     specifically BECAUSE root_boundary_facet() already ruled
		//!     out the seed-crossing case the formula can't detect on its
		//!     own (see fast_candidate()'s own comment for the measured
		//!     99.9985% agreement rate with the geometric path once gated
		//!     this way, versus the unsafe ~1-6% false-verify rate of an
		//!     earlier, weaker gate).
		//!  3. Root cells (every facet is a seed boundary by definition),
		//!     or a confirmed boundary facet: glpt_neighbor_geometric()
		//!     (see its own comment) -- reconstructs real coordinates and
		//!     relocates via the paper's own search() procedure into the
		//!     adjacent seed. This remains the only floating-point path,
		//!     now reached only when a seed change is actually needed.
		result neighbor(int i, glpt& r) const {
			assert(i>=0 && i<=DIM);
			if(!is_root()) {
				bool c0 = is_child0();
				int lminus = (level()==0) ? (DIM-1) : (level()-1);
				int lstar = lminus+1;
				if((c0 && i==DIM) || (!c0 && i==lstar)) {
					r = parent().child(c0 ? 1 : 0);
					return OK;
				}
				int orig_excluded;
				if(!root_boundary_facet(i, orig_excluded)) {
					r = fast_candidate(i);
					return OK;
				}
			}
			glpt_gaifullin_init();
			if(glpt_neighbor_geometric(*this, i, r)) return OK;
			return GLPT_DEEP_CROSSING_UNIMPLEMENTED;
		}
};

// --- Free functions implementing the geometric fallback (see the
// GEOMETRIC FALLBACK comment above glpt::vertex_weights).

// Lemma 7 of the paper (page 21, Fig. 9's own recursion): updates
// alpha (barycentric weights relative to the CURRENT simplex, DIM+1
// of them) in place for descending into child `bit` (0 or 1) of a
// simplex at `level` (its level BEFORE incrementing).
inline void glpt_bary_descend(double alpha[glpt::DIM+1], int level, int bit) {
	const int D = glpt::DIM;
	if(bit==0) {
		double al = alpha[level];
		alpha[D] -= al;
		alpha[level] = 2.0*al;
	} else {
		double ad = alpha[D];
		double al = alpha[level];
		double ap[glpt::DIM+1];
		for(int j=0;j<=D;++j) ap[j]=alpha[j];
		ap[D] = 2.0*ad;
		ap[level] = al-ad;
		// Sigma'_level: circularly shift positions [level..D] right by 1.
		double last = ap[D];
		for(int j=D; j>level; --j) ap[j]=ap[j-1];
		ap[level] = last;
		for(int j=0;j<=D;++j) alpha[j]=ap[j];
	}
}

// Point-location (paper's own search(), Fig. 9, minus findRoot() since
// we always already know which seed to start from): descends `seed`'s
// tree from its root along the point whose barycentric weights
// relative to `seed`'s own DIM+1 original vertices are `alpha` (NOTE:
// modified in place), stopping after exactly `target_level` bisections
// (glpt has no "leaf" concept of its own -- the caller picks the depth,
// typically to match some other cell's simplex_level()).
inline glpt glpt_locate(int seed, double alpha[glpt::DIM+1], int target_level) {
	glpt cur(seed);
	for(int step=0; step<target_level; ++step) {
		int lvl = cur.level();
		int bit = (alpha[lvl] <= alpha[glpt::DIM]) ? 0 : 1;
		glpt_bary_descend(alpha, lvl, bit);
		cur = cur.child(bit);
	}
	return cur;
}

// The fallback itself: locates the cell, in the Gaifullin cell
// sharing `c`'s facet opposite its own local vertex `local_excluded_i`
// (always 0 in current callers -- see neighbor()), at the SAME
// simplex_level as `c`, that occupies the matching region across that
// facet. Returns false only if the Gaifullin adjacency table itself
// has no entry there (should never happen -- CP^2 is closed, verified
// 2026-09-22 that all 108*5 facets pair up).
// Irrational, pairwise-Q-linearly-independent weights (sqrt of the
// first DIM primes) for combining the DIM non-excluded vertex rows
// into a facet point -- see glpt_locate_across_facet's own comment for
// why: a UNIFORM (or any fixed RATIONAL) combination repeatedly landed
// exactly on OTHER internal search() ties, not just the intended
// facet boundary itself (found empirically, 2026-09-22, on
// root.child(0).child(1) -- Delta_0's own reference-simplex pattern is
// built entirely from dyadic rationals, so a rational combination can
// coincide with one of search()'s own alpha[level]==alpha[DIM] tests
// at some step of the descent). Since every step of glpt_bary_descend
// (Lemma 7) is a RATIONAL operation (halving, subtraction, a
// permutation) on the starting alpha, and sqrt(2),sqrt(3),sqrt(5),
// sqrt(7) are linearly independent over Q, no rational combination of
// them can exactly equal another at ANY step of the descent -- an
// exact tie becomes not just unlikely but algebraically impossible
// (short of an actual bug elsewhere).
static const double GLPT_IRRATIONAL_W[glpt::DIM] = {
	1.4142135623730951, 1.7320508075688772, 2.2360679774997896, 2.6457513110645907 // sqrt(2,3,5,7)
};
// Unified geometric neighbor computation, used by neighbor() for EVERY
// case (root, same-seed, and cross-seed) after the Theorem-2/Gamma_SWP
// fast formula was found to be genuinely buggy (2026-09-22: 610/2651 =
// 23.0% of its results had a barycentric weight outside [0,1] against
// the SOURCE seed's own vertices -- a geometric impossibility for a
// Maubach-bisection cell, since every such cell's vertices are convex
// combinations of the seed's original DIM+1 vertices by construction;
// the geometric approach below had 0/2651 such violations, and its two
// ingredients -- vertex_weights()/Theorem 1+Lemma 6, and
// glpt_locate()/Lemma 7 -- were separately validated against an
// independent from-scratch Maubach-bisection simulation and a
// self-location test, so they are trusted here). This trades the
// O(1)/O(d) cost of the same-orthant fast path for O(simplex_level):
// every neighbor query re-descends from a seed root. Restoring an O(1)
// fast path by finding and fixing the actual bug in the ported formula
// is possible future work, not attempted here.
inline bool glpt_neighbor_geometric(const glpt& c, int local_excluded_i, glpt& result) {
	double w[glpt::DIM+1][glpt::DIM+1];
	c.vertex_weights(w);

	// Weights of a facet point, in TWO coordinate systems that must NOT
	// be confused (an earlier draft of this function conflated them --
	// both happen to be DIM+1-vectors, so the bug compiled and ran
	// without any type error, just silently wrong results):
	//  - facet_c[k]: relative to C'S OWN DIM+1 vertices (k=0..DIM).
	//    local_excluded_i is 0 here BY DEFINITION (that's what "exclude
	//    vertex i" means); the other DIM entries are combined with
	//    GLPT_IRRATIONAL_W (not uniform -- see that array's comment: a
	//    single vertex of c can carry identical weight at two different
	//    seed positions, and only decorrelating THIS combination step,
	//    not the final seed-relative vector, breaks that).
	//  - facet[j]: the SAME physical point, but relative to the SEED's
	//    own DIM+1 vertices (j=0..DIM) -- obtained by using facet_c as
	//    the coefficients of a linear combination of w's own ROWS
	//    (w[k] IS vertex k of c, already expressed in the seed's own
	//    weight-numbering).
	double facet_c[glpt::DIM+1];
	{
		double wsum=0.0; int slot=0;
		for(int k=0;k<=glpt::DIM;++k) {
			if(k==local_excluded_i) { facet_c[k]=0.0; continue; }
			facet_c[k] = GLPT_IRRATIONAL_W[slot++];
			wsum += facet_c[k];
		}
		for(int k=0;k<=glpt::DIM;++k) if(k!=local_excluded_i) facet_c[k]/=wsum;
	}
	double facet[glpt::DIM+1];
	for(int j=0;j<=glpt::DIM;++j) {
		double s=0.0;
		for(int k=0;k<=glpt::DIM;++k) s += facet_c[k]*w[k][j];
		facet[j]=s;
	}

	// If some position j in the SEED's own DIM+1 numbering has ~zero
	// weight here, the crossed facet lies exactly on that seed-bounding
	// facet -- a genuine seed-to-seed crossing. Otherwise the crossed
	// facet is strictly interior to this seed's body: an ordinary
	// same-seed neighbor. Which case applies is found numerically, not
	// assumed by the caller -- unlike the old crossing-only version of
	// this function, there is no precondition to assert here.
	int orig_excluded=-1; double best=1e300;
	for(int j=0;j<=glpt::DIM;++j) {
		double a = facet[j]<0?-facet[j]:facet[j];
		if(a<best) { best=a; orig_excluded=j; }
	}

	// A point EXACTLY on the shared facet is inherently ambiguous for a
	// deterministic descent (found empirically, 2026-09-22: search()
	// can retrace the SOURCE cell's own path instead of the neighbor's
	// at intermediate steps, not just tie at the very last one). Fix:
	// nudge the query off the facet, PAST local_excluded_i -- done in
	// C's OWN coordinate system (giving C's own excluded vertex a small
	// negative weight -- algebraically "beyond" it, i.e. past the face
	// opposite it, exactly the direction the crossing goes), NOT applied
	// directly to facet[]'s seed-relative indices (those are a DIFFERENT
	// DIM+1 numbering). Transformed through w[][] exactly like facet_c
	// was.
	//
	// eps's scale is a genuine tension, not a free choice: too small
	// (originally 2^-(simplex_level+12)) and, at simplex_level past
	// roughly 35, the perturbation is swamped by the ~2^-52-relative
	// rounding noise of the w[][] matrix transform below, silently
	// un-nudging the query and reproducing the exact-tie bug above (this
	// SHOWED UP as neighbor() returning the SOURCE cell itself as its own
	// neighbor -- caught 2026-09-22 via a stress test comparing against
	// fast_candidate()'s candidates, at seed=62, depth=39, i=0). Too
	// large and the nudge can overshoot the cell's own internal features
	// partway through the glpt_locate descent, landing on the wrong cell
	// entirely (caught via an involution sweep: cross via i, cross back
	// via i, must return to start). Note this function is no longer ever
	// called for the SIBLING-direction facet (see neighbor()'s own
	// comment -- that facet sits at a finer scale than any single
	// formula here handled well, and is resolved combinatorially
	// instead, with no eps at all); the margin below was tuned and swept
	// (glpt_test.cpp's own stress scripts, not the committed suite) only
	// against the OTHER facets, over depths up to the representable
	// maximum (43), and 2^-(simplex_level-1) minimizes both failure
	// modes simultaneously there; a few (3 in 241387 sampled cases)
	// failures remain, ALL at the exact maximum representable
	// simplex_level (43 = MAX_ORTHANT_LEVEL*DIM+DIM-1), where the
	// orthant fields are completely full and there is no remaining
	// headroom by construction -- not reachable at the depths
	// riemann_cp2.cpp actually uses (peaks around 18-20).
	int target_level = c.simplex_level();
	int eps_shift = target_level-1;
	if(eps_shift<1) eps_shift=1; // target_level<=0 never reaches glpt_locate's
	                             // loop body at all (0 descent steps), and
	                             // target_level==1 needs a sane minimum --
	                             // clamp rather than let the shift go <=0.
	double eps = 1.0;
	for(int s=0;s<eps_shift;++s) eps *= 0.5; // eps = 2^-max(1,target_level-1)
	double query_c[glpt::DIM+1];
	for(int k=0;k<=glpt::DIM;++k) query_c[k] = (k==local_excluded_i) ? -eps : facet_c[k]*(1.0+eps);

	// Case A (best~0): genuine seed crossing -- relocate the query into
	// the ADJACENT seed via the Gaifullin table. The shared facet's
	// vertex order is IDENTICAL in both seeds (verified 2026-09-22), so
	// query[]'s DIM shared-facet weights (computed below, still relative
	// to the SOURCE seed's own vertices) are ALSO directly correct
	// relative to the TARGET seed's vertices -- no relabeling needed;
	// glpt_locate()/glpt_bary_descend() work purely on the abstract
	// DIM+1 weight numbers with no reference to actual coordinates.
	// Case B (best not~0): ordinary same-seed neighbor -- the crossed
	// facet never leaves this seed's body, so the query stays relative
	// to c's own seed_index(); w[][] maps ANY point (not just vertices)
	// from c's own barycentric frame to the seed's, so this is exactly
	// as valid here as it is in case A.
	int target_seed;
	if(best<1e-9) {
		glpt_gaifullin_init();
		int nb_seed = glpt_gaifullin_neighbor[c.seed_index()][orig_excluded];
		if(nb_seed<0) return false;
		target_seed = nb_seed;
	} else {
		target_seed = c.seed_index();
	}

	double query[glpt::DIM+1];
	for(int j=0;j<=glpt::DIM;++j) {
		double s=0.0;
		for(int k=0;k<=glpt::DIM;++k) s += query_c[k]*w[k][j];
		query[j]=s;
	}

	result = glpt_locate(target_seed, query, target_level);
	return true;
}

#endif // GLPT_HPP
