// Self-tests for glpt.hpp. See that file's own header comment for the
// full design rationale; this program verifies everything the design
// claims is SOLID (not the deep-crossing case, which glpt.hpp itself
// declines to implement -- see its GLPT_DEEP_CROSSING_UNIMPLEMENTED).
//
//   1. child()/parent() round-trip: parent(child(x,b)) == x, for many
//      random walks from many random seeds. A basic sanity check that
//      the dimension-generic port from lpt.hpp introduced no bugs.
//   2. Root-level involution: for EVERY (seed, i) pair (108*5=540
//      combinations), crossing seed A -> B via facet i, then B -> A
//      via the SAME facet index i (valid since the local excluded-
//      vertex index always matches on both sides -- verified
//      2026-09-22, see gaifullin_facet_check.py), must return to the
//      exact original root code.
//   3. i=1..DIM-1 crossing involution: for many random deep (non-root)
//      cells, crossing via i and back via i must return to the original
//      code -- regardless of whether the seed changes along the way.
//      (An earlier version of this test also asserted the seed never
//      changes for these i; that assumption, inherited from Kuhn's
//      hypercube, does not hold for Gaifullin's 108-cell complex, which
//      has no group structure guaranteeing it -- see glpt.hpp's
//      neighbor() comment and this file's own comment above test 3.)
//
// Build: g++ -std=c++11 -I. glpt_test.cpp -o glpt_test && ./glpt_test

#include <cstdio>
#include <cstdlib>
#include "glpt.hpp"

static int g_fail = 0;

#define CHECK(cond, msg) do { \
	if(!(cond)) { printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); ++g_fail; } \
} while(0)

void test_child_parent_roundtrip() {
	printf("--- test 1: child()/parent() round-trip ---\n");
	srand(12345);
	int n_checked = 0;
	for(int seed=0; seed<GLPT_NCELLS; seed+=7) { // sample every 7th seed, still >100 seeds
		glpt root(seed);
		glpt cur = root;
		for(int step=0; step<25; ++step) {
			int b = rand()%2;
			glpt child = cur.child(b);
			glpt back = child.parent();
			CHECK(back==cur, "parent(child(x,b)) != x");
			++n_checked;
			cur = child;
		}
	}
	printf("checked %d child/parent round-trips\n", n_checked);
}

void test_root_level_involution() {
	printf("--- test 2: root-level (Gaifullin) crossing involution ---\n");
	glpt_gaifullin_init();
	int n_checked = 0, n_boundary = 0;
	for(int seed=0; seed<GLPT_NCELLS; ++seed) {
		glpt A(seed);
		for(int i=0; i<=glpt::DIM; ++i) {
			glpt B;
			glpt::result rA = A.neighbor(i, B);
			CHECK(rA==glpt::OK, "root-level neighbor() returned DEEP_CROSSING (should never happen at orthant_level==0)");
			if(rA!=glpt::OK) continue;
			CHECK(B.is_root(), "root-level neighbor isn't itself a root");
			CHECK(B.seed_index()!=seed, "root-level neighbor has the SAME seed_index as source (self-loop?)");
			glpt back;
			glpt::result rB = B.neighbor(i, back);
			CHECK(rB==glpt::OK, "return-crossing (B->A via same i) returned DEEP_CROSSING");
			CHECK(back==A, "root-level involution failed: crossing i then i again didn't return to A");
			++n_checked;
		}
	}
	printf("checked %d root-level (seed,i) crossings, %d flagged as unexpected boundary\n", n_checked, n_boundary);
}

// N^(DIM) is NOT universally "the sibling, same index both ways": per
// Theorem 2, a 0-child reaches its sibling via i=DIM (Gamma_RGT/lft),
// but a 1-child reaches ITS sibling via i=lstar=1+((level-1) mod DIM)
// (Gamma_LFT/rgt) -- i=DIM for a 1-child is an ORDINARY same-orthant
// neighbor (Gamma_NSW), unrelated to the sibling, and (found by this
// very test, 2026-09-22) for a 1-child i=lstar can itself fall inside
// 1..DIM-1 (e.g. level=1 => lstar=1), so THAT index is also excluded
// from the plain-swap check for 1-children specifically -- 0-children
// have no such exception (their i=1..DIM-1 branch never references
// lstar at all, only i=DIM is special for them). A same-index round
// trip is therefore only a valid invariant for the remaining i in
// {1,...,DIM-1} (pure Gamma_SWP,i on both sides, symmetric by
// construction). The sibling relationship (i=DIM from a 0-child,
// i=lstar from a 1-child) is checked separately below, against the
// INDEPENDENTLY-computed sibling (parent().child(1-own_child_bit)),
// trusted from test 1's already-passing parent/child round-trip.
int lstar_of(const glpt& g) {
	int lvl = g.level();
	int lminus = (lvl==0) ? (glpt::DIM-1) : (lvl-1);
	return lminus+1;
}

// NOTE (2026-09-22): this test originally also asserted
// nb.seed_index()==cur.seed_index() for i=1..DIM-1 ("same-orthant"
// neighbors never leave the seed). That assumption was inherited
// directly from Kuhn's hypercube, where it holds ONLY because the d!
// canonical simplices form a closed group under Gamma_SWP -- the exact
// same group-theoretic fact behind the Gamma_SWP formula that was
// found to be buggy for THIS port (23% non-convex results, see
// glpt.hpp's neighbor() comment). Gaifullin's 108 cells have no such
// group closure (that is why root-level crossings need a full 108x5
// adjacency table instead of 2 special-cased facets), so there is no
// structural reason a same-orthant (i=1..DIM-1) crossing must stay
// within one seed at arbitrary depth -- and empirically (verified by
// hand, seed=99 depth=7 i=3: excluding that vertex leaves the
// remaining 4 vertices with EXACTLY zero weight at one of the seed's
// own 5 positions, a bona fide boundary touch, not a numerical
// artifact) it sometimes doesn't. So this test now checks only the
// involution property, which holds regardless of whether a seed change
// occurs.
void test_same_orthant_swap_involution() {
	printf("--- test 3: i=1..DIM-1 crossing involution (excl. 1-child's own lstar) ---\n");
	srand(54321);
	int n_checked = 0, n_deep_skipped = 0, n_seed_changed = 0;
	for(int trial=0; trial<2000; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = 3 + rand()%15;
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);
		bool c0 = cur.is_child0();
		int ls = lstar_of(cur);
		for(int i=1; i<glpt::DIM; ++i) { // strictly < DIM: excludes the 0-child sibling-only DIM case
			if(!c0 && i==ls) continue;   // excludes the 1-child sibling-only lstar case
			glpt nb;
			glpt::result r1 = cur.neighbor(i, nb);
			if(r1!=glpt::OK) { ++n_deep_skipped; continue; }
			if(nb.seed_index()!=cur.seed_index()) ++n_seed_changed;
			glpt back;
			glpt::result r2 = nb.neighbor(i, back);
			CHECK(r2==glpt::OK, "return-crossing for i=1..DIM-1 unexpectedly needs deep-crossing");
			if(r2==glpt::OK) CHECK(back==cur, "i=1..DIM-1 crossing involution failed");
			++n_checked;
		}
	}
	printf("checked %d i=1..DIM-1 involutions, %d unexpectedly needed deep-crossing, %d legitimately changed seed\n",
		n_checked, n_deep_skipped, n_seed_changed);
}

void test_sibling_relationship() {
	printf("--- test 4: sibling relationship (N^(DIM) for 0-children, N^(lstar) for 1-children) ---\n");
	srand(24680);
	int n_checked = 0;
	for(int trial=0; trial<2000; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = 1 + rand()%18;
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);
		if(cur.is_root()) continue;

		bool c0 = cur.is_child0();
		glpt expected_sibling = cur.parent().child(c0 ? 1 : 0);

		glpt got;
		glpt::result r;
		if(c0) r = cur.neighbor(glpt::DIM, got);
		else   r = cur.neighbor(lstar_of(cur), got);

		CHECK(r==glpt::OK, "sibling-direction neighbor() unexpectedly needs deep-crossing");
		if(r==glpt::OK) CHECK(got==expected_sibling, "N^(sibling index) != parent().child(other)");
		++n_checked;
	}
	printf("checked %d sibling relationships\n", n_checked);
}

void test_deep_i0_frequency() {
	// Not a correctness test -- just reports how often i=0 neighbor()
	// used to hit the (now-implemented) deep-crossing fallback in
	// practice, to gauge how much of the mesh that path actually
	// serves.
	printf("--- info: i=0 same-orthant vs geometric-fallback frequency ---\n");
	srand(99999);
	int n_same_orthant=0, n_fallback=0, n_fail=0;
	for(int trial=0; trial<5000; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = 1 + rand()%20;
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);
		if(cur.is_root()) continue;
		glpt nb;
		glpt::result r = cur.neighbor(0, nb);
		if(r!=glpt::OK) { ++n_fail; continue; }
		if(nb.seed_index()==cur.seed_index()) ++n_same_orthant; else ++n_fallback;
	}
	CHECK(n_fail==0, "neighbor(0) failed outright somewhere in the frequency sweep");
	printf("i=0 neighbor(): %d resolved same-orthant, %d resolved via geometric fallback, %d failed (out of %d)\n",
		n_same_orthant, n_fallback, n_fail, n_same_orthant+n_fallback+n_fail);
}

// Convexity-invariant stress test. There is no longer a second,
// independent implementation to cross-check neighbor() against (the
// Theorem-2/Gamma_SWP fast formula was removed as buggy -- see
// glpt.hpp's neighbor() comment, and glpt_test.cpp's own history in
// git for the comparison that caught it). The strongest remaining
// correctness signal is the one that originally exposed that bug: every
// vertex of every cell anywhere in a seed's Maubach-bisection tree must
// be a convex combination of that SEED's own DIM+1 original vertices
// (all barycentric weights in [0,1]), since the tree is built entirely
// by repeated midpointing. This test computes many neighbors, across
// every facet index and a range of depths, and checks that invariant on
// the RESULT's own vertex_weights() -- relative to the result's OWN
// seed_index(), which is exactly what would be violated if neighbor()
// ever mixed up which seed a query point belongs to.
void test_neighbor_convexity_invariant() {
	printf("--- test 5: neighbor() result convexity invariant, all i, all depths ---\n");
	srand(13579);
	int n_checked = 0, n_deep_skipped = 0, n_violations = 0;
	for(int trial=0; trial<3000; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = rand()%25;
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);
		for(int i=0; i<=glpt::DIM; ++i) {
			glpt nb;
			glpt::result r = cur.neighbor(i, nb);
			if(r!=glpt::OK) { ++n_deep_skipped; continue; }
			double w[glpt::DIM+1][glpt::DIM+1];
			nb.vertex_weights(w);
			bool ok = true;
			for(int k=0;k<=glpt::DIM && ok;++k)
				for(int j=0;j<=glpt::DIM && ok;++j)
					if(w[k][j] < -1e-9 || w[k][j] > 1.0+1e-9) ok = false;
			if(!ok) ++n_violations;
			CHECK(ok, "neighbor() result has a barycentric weight outside [0,1]");
			++n_checked;
		}
	}
	printf("checked %d neighbor() results, %d unexpectedly needed deep-crossing, %d convexity violations\n",
		n_checked, n_deep_skipped, n_violations);
}

void test_deep_crossing_involution() {
	printf("--- test 6: deep-crossing (geometric fallback) involution ---\n");
	srand(97531);
	int n_checked = 0, n_same_orthant_skipped = 0;
	for(int trial=0; trial<20000; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = 1 + rand()%36; // push well past round 0, deep into orthant territory
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);
		if(cur.is_root()) continue;

		glpt nb;
		glpt::result r1 = cur.neighbor(0, nb);
		CHECK(r1==glpt::OK, "neighbor(0) failed in deep-crossing involution sweep");
		if(r1!=glpt::OK) continue;
		if(nb.seed_index()==cur.seed_index()) { ++n_same_orthant_skipped; continue; } // not a crossing this time

		glpt back;
		glpt::result r2 = nb.neighbor(0, back);
		CHECK(r2==glpt::OK, "return-crossing (neighbor's own i=0) failed");
		if(r2==glpt::OK) CHECK(back==cur, "deep-crossing involution failed: A->B->A != A");
		++n_checked;
	}
	printf("checked %d deep-crossing (seed-changing) involutions, %d same-orthant (not a crossing this trial)\n",
		n_checked, n_same_orthant_skipped);
}

// Full-depth-range stress test, covering EVERY i (0..DIM, not just the
// specific ones tests 3/4/6 each focus on) at depths spanning the
// entire representable range (0..43=MAX_ORTHANT_LEVEL*DIM+DIM-1). This
// is what actually caught two real bugs during development (2026-09-22):
//  - eps too small in glpt_neighbor_geometric's nudge caused neighbor()
//    to return the SOURCE cell as its own neighbor past simplex_level
//    ~35 (floating-point rounding in the w[][] transform swamped the
//    nudge); eps too large caused wrong-cell overshoot at shallow-to-
//    mid depths instead -- the margin now used was swept against both
//    failure modes at once.
//  - routing the SIBLING-direction facet (see neighbor()'s own comment)
//    through the same geometric nudge broke ~5-15% of cases (that facet
//    sits at a finer scale than the margin above is tuned for); fixed
//    by resolving it as a plain parent()/child() combinatorial fact
//    instead, with no floating point at all.
// IMPORTANT test-methodology note, itself found the hard way: "cross
// via i, then cross back via the SAME i" is only a valid involution
// check for the NON-sibling direction. For the sibling direction, the
// return crossing is generally via a DIFFERENT index (i==DIM only
// undoes to i==lstar or vice versa, depending on child parity -- see
// lstar_of()), so the sibling direction here is checked against the
// exact parent()/child() identity instead (same as test 4), and
// excluded from the plain involution check (same as test 3).
void test_full_depth_stress() {
	// Depth covers the FULL representable range, 0..43=MAX_ORTHANT_
	// LEVEL*DIM+DIM-1 inclusive. An earlier version of this test capped
	// one short of 43: at that exact depth, glpt_neighbor_geometric's
	// nudge-epsilon margin had a documented, narrow, very-low-frequency
	// residual involution failure. That is gone as of the exact-integer
	// root_boundary_facet() gate (2026-09-22, see neighbor()'s own
	// comment): non-sibling facets it calls interior now go through
	// fast_candidate() instead, which has no floating point at all, so
	// glpt_neighbor_geometric is only ever reached for genuine seed
	// crossings -- verified by a 750000-check sweep (not this committed
	// suite) at depth 0..43 inclusive with 0 failures of any kind,
	// hence the cap is removed here too.
	printf("--- test 7: full-depth-range stress (all i, depth 0..43) ---\n");
	srand(31415);
	long n_total=0, n_sib_checked=0, n_nonsib_checked=0, n_deep_skipped=0;
	for(int trial=0; trial<80000; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = rand()%(glpt::MAX_ORTHANT_LEVEL*glpt::DIM+glpt::DIM);
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);
		bool c0 = cur.is_root() ? false : cur.is_child0();
		int ls = lstar_of(cur);
		for(int i=0;i<=glpt::DIM;++i) {
			bool sibling_dir = !cur.is_root() && ((c0 && i==glpt::DIM) || (!c0 && i==ls));
			++n_total;
			glpt nb;
			glpt::result r = cur.neighbor(i, nb);
			if(r!=glpt::OK) { ++n_deep_skipped; continue; }
			if(sibling_dir) {
				++n_sib_checked;
				CHECK(nb==cur.parent().child(c0?1:0), "sibling-direction neighbor() != parent().child(other)");
				continue;
			}
			++n_nonsib_checked;
			CHECK(!(nb==cur), "neighbor() returned the source cell as its own neighbor");
			glpt back;
			glpt::result r2 = nb.neighbor(i, back);
			CHECK(r2==glpt::OK, "return-crossing failed in full-depth stress");
			if(r2==glpt::OK) CHECK(back==cur, "full-depth involution failed: A->B->A != A");
		}
	}
	printf("checked %ld total (%ld sibling, %ld non-sibling), %ld unexpectedly needed deep-crossing\n",
		n_total, n_sib_checked, n_nonsib_checked, n_deep_skipped);
}

// Directly validates the exact-integer gate itself (root_boundary_facet
// + fast_candidate), not just its end-to-end effect through neighbor()
// (test 7 exercises that, but would not by itself catch a bug in this
// mechanism that happened to still produce a self-consistent involution
// -- e.g. a wrong-but-still-invertible pairing, in the spirit of the
// original traces_to_own_root() false positive this design replaced).
// Two independent checks against the SAME (cell,i) pair, for every i,
// at every depth (0..43 inclusive), skipping the sibling direction
// (handled separately, see neighbor()'s own comment):
//  1. root_boundary_facet()=false (interior) must mean
//     glpt_neighbor_geometric() stays in the SAME seed -- i.e. no false
//     negatives: the gate never lets a genuine crossing through as
//     "interior".
//  2. Whenever root_boundary_facet()=false, fast_candidate() must equal
//     glpt_neighbor_geometric()'s result exactly -- the actual safety
//     property neighbor() relies on.
void test_exact_boundary_gate() {
	printf("--- test 8: exact boundary-facet gate vs geometric ground truth ---\n");
	// Depth capped one short of the absolute maximum (43): this test
	// calls glpt_neighbor_geometric() DIRECTLY as the ground truth (test
	// 7 goes through neighbor(), which -- as of this same gate -- never
	// reaches glpt_neighbor_geometric for a facet root_boundary_facet()
	// calls interior, so it doesn't hit this). AT depth 43 exactly,
	// glpt_neighbor_geometric() has its own separate, already-documented
	// floating-point edge case (see its eps comment); comparing
	// fast_candidate() against it there would be checking fast_candidate
	// against a wrong answer. Test 7's 0-failure full-range (0..43
	// inclusive) sweep through neighbor() itself is the actual evidence
	// fast_candidate is fine at that depth -- this test's job is
	// establishing that agreement in general, not re-litigating depth 43.
	srand(2222);
	long n_interior=0, n_boundary=0;
	for(int trial=0; trial<40000; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = rand()%(glpt::MAX_ORTHANT_LEVEL*glpt::DIM+glpt::DIM-1);
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);
		bool c0 = cur.is_root() ? false : cur.is_child0();
		int ls = lstar_of(cur);
		for(int i=0;i<=glpt::DIM;++i) {
			if(!cur.is_root() && ((c0 && i==glpt::DIM) || (!c0 && i==ls))) continue; // sibling, not this gate's concern
			int orig_excluded;
			bool is_boundary = cur.is_root() ? true : cur.root_boundary_facet(i, orig_excluded);
			glpt geo;
			if(!glpt_neighbor_geometric(cur, i, geo)) continue;
			if(is_boundary) {
				++n_boundary;
				CHECK(geo.seed_index()!=cur.seed_index(), "root_boundary_facet() said boundary but geometric result stayed in the same seed");
			} else {
				++n_interior;
				glpt fast = cur.fast_candidate(i);
				CHECK(fast==geo, "fast_candidate() disagrees with glpt_neighbor_geometric() on a facet root_boundary_facet() called interior");
			}
		}
	}
	printf("checked %ld interior + %ld boundary facets against the geometric ground truth\n", n_interior, n_boundary);
}

int main() {
	printf("sizeof(glpt) = %zu bytes\n", sizeof(glpt));
	printf("GLPT_NCELLS = %d\n", GLPT_NCELLS);
	printf("glpt::MAX_ORTHANT_LEVEL = %d\n\n", glpt::MAX_ORTHANT_LEVEL);

	test_child_parent_roundtrip();
	test_root_level_involution();
	test_same_orthant_swap_involution();
	test_sibling_relationship();
	test_neighbor_convexity_invariant();
	test_deep_crossing_involution();
	test_full_depth_stress();
	test_exact_boundary_gate();
	test_deep_i0_frequency();

	printf("\n%s (%d failures)\n", g_fail==0 ? "ALL CHECKS PASSED" : "CHECKS FAILED", g_fail);
	return g_fail==0 ? 0 : 1;
}
