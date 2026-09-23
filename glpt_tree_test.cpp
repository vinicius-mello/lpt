// Self-tests for glpt_tree.hpp. See that file's own header comment for
// the design rationale (leaves-only storage, backward-shift deletion,
// and the deliberately limited scope of neighbor_leaf()'s stitching).
//
//   1. Basic seed/bisect sanity.
//   2. Heavy insert/delete churn (many random bisections, crossing
//      several rehashes): leaf_count() and for_each_leaf() must always
//      exactly match an independently-tracked ground-truth set, and
//      every tracked leaf must exists() -- the actual correctness test
//      for the open-addressing table + backward-shift deletion under
//      the specific insert/delete pattern this file's whole design is
//      built around (one delete, two inserts, every bisection).
//   3. neighbor_leaf() transitions: FOUND when the same-level candidate
//      is itself a leaf; REFINED_FURTHER once it's bisected; and the
//      finer side (either child) can find its way back to the coarser
//      leaf via SOME facet index -- checked directly, not assumed.
//   4. compat_bisect(): after many compat_bisect()-only operations on a
//      growing mesh, exhaustively checks the 1-level-graded invariant
//      it's supposed to maintain -- for EVERY current leaf and EVERY
//      non-skipped facet, the same-level neighbor candidate is at most
//      one parent() step from an existing leaf. Ported from lpt.hpp's
//      own lpt_tree<lpt>::compat_bisect/compat_bisect_rec; this test is
//      the actual verification that the port is correct for Gaifullin's
//      seeds too, not just assumed from the hypercube case.
//
// Build: g++ -std=c++11 -I. glpt_tree_test.cpp lpt.o -o glpt_tree_test

#include <cstdio>
#include <cstdlib>
#include <set>
#include <vector>
#include "glpt_tree.hpp"

static int g_fail = 0;
#define CHECK(cond, msg) do { if(!(cond)) { ++g_fail; printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); } } while(0)

void test_basic_seed_bisect() {
	printf("--- test 1: basic seed/bisect sanity ---\n");
	glpt_tree tree;
	tree.seed_all_roots();
	CHECK(tree.leaf_count()==size_t(GLPT_NCELLS), "seed_all_roots() didn't produce 108 leaves");
	for(int s=0;s<GLPT_NCELLS;++s) CHECK(tree.exists(glpt(s)), "a seeded root isn't exists()");

	glpt r(42);
	glpt c0 = r.child(0), c1 = r.child(1);
	CHECK(!tree.exists(c0) && !tree.exists(c1), "children exist before bisect()");
	tree.bisect(r);
	CHECK(tree.leaf_count()==size_t(GLPT_NCELLS)+1, "leaf_count() didn't increase by 1 after one bisect()");
	CHECK(!tree.exists(r), "parent still exists() after bisect()");
	CHECK(tree.exists(c0) && tree.exists(c1), "children don't exists() after bisect()");
	printf("ok\n");
}

void test_heavy_churn() {
	printf("--- test 2: heavy insert/delete churn vs ground truth ---\n");
	glpt_tree tree;
	tree.seed_all_roots();
	std::set<uint64_t> truth;
	for(int s=0;s<GLPT_NCELLS;++s) truth.insert(glpt(s).raw());

	srand(9999);
	size_t max_buckets_seen = tree.bucket_count();
	int n_bisections = 0;
	for(int round=0; round<30000; ++round) {
		// pick a random element of `truth` (linear scan -- fine, truth
		// is small relative to iteration count, this is a test not
		// production code) whose depth is still shallow enough to bisect.
		int skip = rand() % (int)truth.size();
		std::set<uint64_t>::iterator it = truth.begin();
		for(int k=0;k<skip;++k) ++it;
		glpt r = glpt::from_raw(*it);
		if(r.orthant_level()>=glpt::MAX_ORTHANT_LEVEL-1) continue; // stay well clear of the depth cap
		truth.erase(it);
		glpt c0=r.child(0), c1=r.child(1);
		truth.insert(c0.raw());
		truth.insert(c1.raw());
		tree.bisect(r);
		++n_bisections;
		max_buckets_seen = tree.bucket_count()>max_buckets_seen ? tree.bucket_count() : max_buckets_seen;
	}
	CHECK(tree.leaf_count()==truth.size(), "leaf_count() disagrees with the independently-tracked ground truth");
	int n_missing=0;
	for(std::set<uint64_t>::iterator it=truth.begin(); it!=truth.end(); ++it)
		if(!tree.exists(glpt::from_raw(*it))) ++n_missing;
	CHECK(n_missing==0, "a ground-truth leaf isn't exists() in the tree");

	std::set<uint64_t> visited;
	struct Collector {
		std::set<uint64_t>* out;
		void operator()(const glpt& c) const { out->insert(c.raw()); }
	};
	Collector col; col.out=&visited;
	tree.for_each_leaf(col);
	CHECK(visited==truth, "for_each_leaf() doesn't match the ground truth set exactly");

	printf("bisections=%d final_leaves=%zu buckets: start=257 grew_to=%zu (load factor %.2f)\n",
		n_bisections, tree.leaf_count(), tree.bucket_count(), tree.load_factor());
	CHECK(tree.bucket_count()>257, "table never grew despite thousands of insertions -- rehash path untested");
}

int lstar_of(const glpt& g) {
	int lvl = g.level();
	int lminus = (lvl==0) ? (glpt::DIM-1) : (lvl-1);
	return lminus+1;
}

void test_neighbor_leaf_transitions() {
	printf("--- test 3: neighbor_leaf() FOUND / REFINED_FURTHER transitions ---\n");
	srand(24681);
	int n_found_direct=0, n_refined_further=0, n_finer_side_roundtrip_ok=0, n_finer_side_roundtrip_checked=0;
	for(int trial=0; trial<3000; ++trial) {
		glpt_tree tree;
		tree.seed_all_roots();
		int seed = rand()%GLPT_NCELLS;
		glpt r(seed);
		int depth = rand()%8;
		for(int s=0;s<depth;++s) {
			glpt child = r.child(rand()%2);
			tree.bisect(r);
			r = child;
		}
		// r is now a leaf (its ancestors were bisected down to it).
		CHECK(tree.exists(r), "r isn't a leaf after the setup bisection chain");

		int i = rand()%(glpt::DIM+1);
		glpt n;
		glpt::result gr = r.neighbor(i, n);
		if(gr!=glpt::OK) continue;

		glpt result;
		glpt_tree::neighbor_status st = tree.neighbor_leaf(r, i, result);
		if(tree.exists(n)) {
			CHECK(st==glpt_tree::FOUND && result==n, "neighbor_leaf() didn't find an existing same-level leaf directly");
			++n_found_direct;

			// Now bisect n (if that's safe -- stay clear of the depth cap)
			// and confirm the status flips to REFINED_FURTHER.
			if(n.orthant_level()<glpt::MAX_ORTHANT_LEVEL-1) {
				tree.bisect(n);
				glpt result2;
				glpt_tree::neighbor_status st2 = tree.neighbor_leaf(r, i, result2);
				CHECK(st2==glpt_tree::REFINED_FURTHER, "neighbor_leaf() didn't report REFINED_FURTHER after the neighbor was bisected");
				++n_refined_further;

				// From EITHER of n's new children, some facet index should
				// find its way back to r (locate_leaf walking up from a
				// same-level-as-child candidate should reach r).
				++n_finer_side_roundtrip_checked;
				bool roundtrip_ok=false;
				glpt children[2] = { n.child(0), n.child(1) };
				for(int cc=0; cc<2 && !roundtrip_ok; ++cc) {
					for(int j=0; j<=glpt::DIM && !roundtrip_ok; ++j) {
						glpt back;
						glpt_tree::neighbor_status st3 = tree.neighbor_leaf(children[cc], j, back);
						if(st3==glpt_tree::FOUND && back==r) roundtrip_ok=true;
					}
				}
				CHECK(roundtrip_ok, "neither of the bisected neighbor's children can find their way back to r via any facet");
				if(roundtrip_ok) ++n_finer_side_roundtrip_ok;
			}
		} else {
			// n doesn't exist as a leaf -- either coarser (some ancestor
			// is a leaf) or, in this fresh-tree setup, that shouldn't
			// happen since we only ever bisect ALONG r's own chain. Just
			// confirm locate_leaf()'s own contract: st must be FOUND with
			// an ancestor-or-self of n, or REFINED_FURTHER, never UNRESOLVED.
			CHECK(st==glpt_tree::FOUND || st==glpt_tree::REFINED_FURTHER, "neighbor_leaf() returned UNRESOLVED unexpectedly");
			if(st==glpt_tree::FOUND) {
				CHECK(tree.exists(result), "neighbor_leaf() FOUND result doesn't exists()");
				glpt anc = n;
				bool is_ancestor=false;
				while(true) { if(anc==result) { is_ancestor=true; break; } if(anc.is_root()) break; anc=anc.parent(); }
				CHECK(is_ancestor, "neighbor_leaf() FOUND result isn't an ancestor-or-self of the same-level candidate");
			}
		}
	}
	printf("found_direct=%d refined_further=%d finer_side_roundtrip_ok=%d/%d\n",
		n_found_direct, n_refined_further, n_finer_side_roundtrip_ok, n_finer_side_roundtrip_checked);
}

void test_compat_bisect() {
	printf("--- test 4: compat_bisect() maintains the 1-level-graded invariant ---\n");
	glpt_tree tree;
	tree.seed_all_roots();
	srand(135791);

	std::vector<uint64_t> leaves;
	struct Collector { std::vector<uint64_t>* out; void operator()(const glpt& c) const { out->push_back(c.raw()); } };
	for(int s=0;s<GLPT_NCELLS;++s) leaves.push_back(glpt(s).raw());

	int n_ops=0, n_skipped=0;
	for(int round=0; round<4000; ++round) {
		if(round%50==0) {
			leaves.clear();
			Collector col; col.out=&leaves;
			tree.for_each_leaf(col);
		}
		int idx = rand() % (int)leaves.size();
		glpt r = glpt::from_raw(leaves[idx]);
		if(!tree.exists(r)) { ++n_skipped; continue; } // consumed by an earlier round's own recursion
		if(r.orthant_level()>=glpt::MAX_ORTHANT_LEVEL-2) { ++n_skipped; continue; } // stay clear of the depth cap
		tree.compat_bisect(r);
		++n_ops;
	}
	printf("compat_bisect() calls=%d skipped=%d final leaf_count=%zu\n", n_ops, n_skipped, tree.leaf_count());

	// The actual invariant check: for EVERY current leaf and EVERY facet
	// index compat_bisect() actually guarantees gradedness for (i.e. NOT
	// i==DIM or i==L.level() -- see compat_bisect()'s own comment on
	// what it deliberately does not cover), the same-level neighbor
	// candidate must be at most one parent() step from an existing leaf.
	std::vector<uint64_t> final_leaves;
	Collector col; col.out=&final_leaves;
	tree.for_each_leaf(col);
	long n_checked=0, n_violations=0;
	for(size_t li=0; li<final_leaves.size(); ++li) {
		glpt L = glpt::from_raw(final_leaves[li]);
		for(int i=0;i<=glpt::DIM;++i) {
			if(i==glpt::DIM || i==L.level()) continue;
			glpt n;
			glpt::result gr = L.neighbor(i, n);
			if(gr!=glpt::OK) continue;
			++n_checked;
			bool ok = tree.exists(n) || (!n.is_root() && tree.exists(n.parent()));
			if(!ok) {
				++n_violations;
				if(n_violations<=5) printf("VIOLATION: leaf level=%d orthlvl=%d seed=%d i=%d -- neighbor and its parent both missing\n",
					L.level(), L.orthant_level(), L.seed_index(), i);
			}
		}
	}
	CHECK(n_violations==0, "compat_bisect() left the mesh more than 1 level ungraded somewhere");
	printf("graded-invariant check: %ld facet checks, %ld violations\n", n_checked, n_violations);
}

void test_recent_leaves() {
	printf("--- test 5: recent_leaves() reports every leaf compat_bisect() creates ---\n");
	glpt_tree tree;
	tree.seed_all_roots();
	srand(24680);

	std::vector<uint64_t> leaves;
	struct Collector { std::vector<uint64_t>* out; void operator()(const glpt& c) const { out->push_back(c.raw()); } };
	for(int s=0;s<GLPT_NCELLS;++s) leaves.push_back(glpt(s).raw());

	long n_recent_total=0, n_not_existing=0, n_odd_count=0;
	for(int round=0; round<500; ++round) {
		if(round%25==0) { leaves.clear(); Collector col; col.out=&leaves; tree.for_each_leaf(col); }
		int idx = rand() % (int)leaves.size();
		glpt r = glpt::from_raw(leaves[idx]);
		if(!tree.exists(r)) continue;
		if(r.orthant_level()>=glpt::MAX_ORTHANT_LEVEL-2) continue;

		size_t before = tree.leaf_count();
		tree.clear_recent();
		tree.compat_bisect(r);
		size_t after = tree.leaf_count();

		const std::vector<glpt>& recent = tree.recent_leaves();
		n_recent_total += (long)recent.size();
		// NOTE (found by this test, 2026-09-23): an entry here need NOT
		// still exists() by the time the whole compat_bisect() call
		// returns -- a cell created early in the cascade (e.g. while
		// forcing a neighbor finer) can itself get bisected AGAIN later
		// in the SAME call (e.g. because forcing ANOTHER neighbor then
		// required it too), leaving it superseded by its own children
		// (which are separately present in recent_ from THAT bisect()
		// call). This is not a bug: recent_leaves() is "every cell
		// inserted since clear_recent()", not "every current leaf
		// inserted since clear_recent()" -- a priority-queue-driven
		// refinement loop naturally handles this exactly the way it
		// already must handle a POPPED entry going stale (check exists()
		// before using each entry; skip if gone), so this is tracked as
		// an informational count, not a failure.
		for(size_t i=0;i<recent.size();++i) if(!tree.exists(recent[i])) ++n_not_existing;
		// Every bisect() call nets exactly +1 leaf and appends 2 entries to
		// recent_ -- so recent_.size() must be exactly 2*(number of
		// bisections), and leaf_count() grows by exactly recent_.size()/2.
		if(recent.size()%2!=0) ++n_odd_count;
		CHECK(after-before == recent.size()/2, "leaf_count() growth doesn't match recent_leaves().size()/2");
	}
	CHECK(n_odd_count==0, "recent_leaves() had an odd number of entries after some compat_bisect() call");
	printf("total recent entries across all calls=%ld, later-superseded(informational)=%ld, odd-count calls=%ld\n",
		n_recent_total, n_not_existing, n_odd_count);
}

int main() {
	test_basic_seed_bisect();
	test_heavy_churn();
	test_neighbor_leaf_transitions();
	test_compat_bisect();
	test_recent_leaves();

	printf("\n%s (%d failures)\n", g_fail==0 ? "ALL CHECKS PASSED" : "CHECKS FAILED", g_fail);
	return g_fail==0 ? 0 : 1;
}
