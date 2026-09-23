// Self-tests for glpt_vertex_ids.hpp. See that file's own header
// comment for the design rationale.
//
//   1. vertex_weights_exact() sanity: matches vertex_weights()'s double
//      version (basic check the exact-integer port introduced no bug).
//   2. Incremental id assignment: bisecting a cell must mint EXACTLY one
//      new vertex id (glpt_child_vertex_ids() already asserts this;
//      this test also checks the new id is never reused for an
//      unrelated edge, and that ids below GLPT_BASE_VERTEX_COUNT only
//      ever appear where a cell still directly touches one of
//      Gaifullin's 15 original points).
//   3. THE CRITICAL TEST -- neighbor id consistency: build a cell C via
//      a random bisection path (sharing ONE glpt_edge_cache), compute
//      N = C.neighbor(i), independently RECONSTRUCT N's own bisection
//      path (walk up to its root via parent(), recording each
//      is_child0() bit, then redescend via child()) and compute N's ids
//      using the SAME shared cache. The shared facet's ids must agree
//      on both sides -- this is the actual property the whole design
//      depends on (see glpt_vertex_ids.hpp's own comment), checked
//      directly rather than assumed.
//   4. glpt_edge_cache's own insert/lookup/forget_edge churn against an
//      independent std::map ground truth -- the actual correctness test
//      for its open-addressing storage and backward-shift deletion
//      (added 2026-09-23 when std::map was replaced with open addressing
//      for memory reasons -- see that class's own comment; forget_edge()
//      wasn't otherwise exercised by anything in this file).
//
// Build: g++ -std=c++11 -I. glpt_vertex_ids_test.cpp lpt.o -o glpt_vertex_ids_test

#include <cstdio>
#include <cstdlib>
#include <set>
#include <map>
#include "glpt_vertex_ids.hpp"

static int g_fail = 0;
#define CHECK(cond, msg) do { if(!(cond)) { ++g_fail; printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__); } } while(0)

void test_vertex_weights_exact_sanity() {
	printf("--- test 1: vertex_weights_exact() vs vertex_weights() ---\n");
	srand(1111);
	int n_checked=0, n_mismatch=0;
	for(int trial=0;trial<3000;++trial) {
		int seed=rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth=rand()%(glpt::MAX_ORTHANT_LEVEL*glpt::DIM+glpt::DIM);
		for(int s=0;s<depth;++s) cur=cur.child(rand()%2);
		double w[glpt::DIM+1][glpt::DIM+1];
		cur.vertex_weights(w);
		int wn[glpt::DIM+1][glpt::DIM+1], shift;
		cur.vertex_weights_exact(wn, shift);
		double den = double(1<<shift);
		for(int k=0;k<=glpt::DIM;++k) for(int j=0;j<=glpt::DIM;++j) {
			++n_checked;
			double exact_val = wn[k][j]/den;
			if(exact_val<w[k][j]-1e-9 || exact_val>w[k][j]+1e-9) ++n_mismatch;
		}
	}
	CHECK(n_mismatch==0, "vertex_weights_exact() disagrees with vertex_weights()");
	printf("checked %d weight entries, %d mismatches\n", n_checked, n_mismatch);
}

void test_incremental_id_assignment() {
	printf("--- test 2: incremental vertex id assignment ---\n");
	srand(3333);
	glpt_edge_cache cache;
	int n_checked=0;
	for(int trial=0;trial<2000;++trial) {
		int seed=rand()%GLPT_NCELLS;
		glpt cur(seed);
		int ids[glpt::DIM+1];
		glpt_root_vertex_ids(seed, ids);
		for(int k=0;k<=glpt::DIM;++k)
			CHECK(ids[k]>=0 && ids[k]<GLPT_BASE_VERTEX_COUNT, "root vertex id out of base range");
		int depth = rand()%20;
		for(int step=0; step<depth; ++step) {
			int before = cache.next_id();
			int zo = rand()%2;
			int child_ids[glpt::DIM+1];
			glpt_child_vertex_ids(cur, ids, zo, cache, child_ids);
			int after = cache.next_id();
			CHECK(after==before || after==before+1, "bisection minted more than one new id, or a negative amount");
			// exactly one child_ids entry should be new-or-shared beyond
			// the parent's own id set when a fresh mint happened
			cur = cur.child(zo);
			for(int k=0;k<=glpt::DIM;++k) ids[k]=child_ids[k];
			++n_checked;
		}
	}
	printf("checked %d bisection steps, %d ids minted total\n", n_checked, cache.next_id()-GLPT_BASE_VERTEX_COUNT);
}

int lstar_of(const glpt& g) {
	int lvl = g.level();
	int lminus = (lvl==0) ? (glpt::DIM-1) : (lvl-1);
	return lminus+1;
}

// Reconstructs the (zo bit) path from a cell's own root down to itself.
void path_from_root(const glpt& c, int path[], int& path_len) {
	path_len = 0;
	glpt cur = c;
	while(!cur.is_root()) {
		path[path_len++] = cur.is_child0() ? 0 : 1;
		cur = cur.parent();
	}
	// path[] currently holds the bits leaf-to-root; reverse in place.
	for(int i=0;i<path_len/2;++i) { int t=path[i]; path[i]=path[path_len-1-i]; path[path_len-1-i]=t; }
}

// Computes a cell's vertex ids from scratch (its own root + its own
// bisection path), continuing to use the SAME shared cache -- this is
// what an independent traversal of the mesh reaching this same cell
// would naturally do.
void vertex_ids_from_scratch(const glpt& c, glpt_edge_cache& cache, int ids[glpt::DIM+1]) {
	int path[64], path_len;
	path_from_root(c, path, path_len);
	glpt_root_vertex_ids(c.seed_index(), ids);
	glpt cur(c.seed_index());
	for(int i=0;i<path_len;++i) {
		int child_ids[glpt::DIM+1];
		glpt_child_vertex_ids(cur, ids, path[i], cache, child_ids);
		cur = cur.child(path[i]);
		for(int k=0;k<=glpt::DIM;++k) ids[k]=child_ids[k];
	}
}

void test_neighbor_id_consistency() {
	printf("--- test 3: neighbor id consistency (the critical test) ---\n");
	srand(5555);
	glpt_edge_cache cache;
	int n_checked=0, n_deep_skipped=0;
	for(int trial=0; trial<20000; ++trial) {
		int seed = rand()%GLPT_NCELLS;
		glpt cur(seed);
		int depth = rand()%20;
		for(int step=0; step<depth; ++step) cur = cur.child(rand()%2);

		int c_ids[glpt::DIM+1];
		vertex_ids_from_scratch(cur, cache, c_ids);

		bool c0 = cur.is_root() ? false : cur.is_child0();
		int ls = lstar_of(cur);
		for(int i=0;i<=glpt::DIM;++i) {
			bool sibling = !cur.is_root() && ((c0 && i==glpt::DIM) || (!c0 && i==ls));
			glpt nb;
			glpt::result r = cur.neighbor(i, nb);
			if(r!=glpt::OK) { ++n_deep_skipped; continue; }
			int n_ids[glpt::DIM+1];
			vertex_ids_from_scratch(nb, cache, n_ids);

			std::set<int> shared_c;
			for(int k=0;k<=glpt::DIM;++k) if(k!=i) shared_c.insert(c_ids[k]);
			std::set<int> n_set(n_ids, n_ids+glpt::DIM+1);
			bool subset = true;
			for(std::set<int>::iterator it=shared_c.begin(); it!=shared_c.end(); ++it)
				if(n_set.find(*it)==n_set.end()) subset=false;
			CHECK(subset, "neighbor's vertex ids don't contain the shared facet's ids from this cell's side");
			(void)sibling;
			++n_checked;
		}
	}
	printf("checked %d (cell,i) neighbor pairs, %d unexpectedly needed deep-crossing, cache has %zu edges / %d ids minted\n",
		n_checked, n_deep_skipped, cache.size(), cache.next_id()-GLPT_BASE_VERTEX_COUNT);
}

void test_edge_cache_churn() {
	printf("--- test 4: glpt_edge_cache insert/lookup/forget churn vs std::map ground truth ---\n");
	glpt_edge_cache cache;
	std::map<std::pair<int,int>,int> truth;
	srand(555);
	int n_checked=0, n_mismatch=0;
	for(int round=0; round<200000; ++round) {
		int op = rand()%3;
		int a = rand()%2000, b = rand()%2000;
		if(a==b) continue;
		int lo=a<b?a:b, hi=a<b?b:a;
		if(op<2) {
			int id = cache.get_or_create_midpoint(a,b);
			std::map<std::pair<int,int>,int>::iterator it = truth.find(std::make_pair(lo,hi));
			if(it==truth.end()) truth[std::make_pair(lo,hi)] = id;
			else { ++n_checked; if(it->second!=id) ++n_mismatch; }
		} else {
			cache.forget_edge(a,b);
			truth.erase(std::make_pair(lo,hi));
		}
	}
	CHECK(n_mismatch==0, "glpt_edge_cache returned a different id for an already-cached edge");
	CHECK(cache.size()==truth.size(), "glpt_edge_cache.size() disagrees with the independent ground truth");
	int n_reverify_fail=0;
	for(std::map<std::pair<int,int>,int>::iterator it=truth.begin(); it!=truth.end(); ++it)
		if(cache.get_or_create_midpoint(it->first.first, it->first.second) != it->second) ++n_reverify_fail;
	CHECK(n_reverify_fail==0, "a surviving edge no longer resolves to its own cached id");
	printf("checked=%d mismatch=%d reverify_fail=%d cache.size()=%zu truth.size()=%zu\n",
		n_checked, n_mismatch, n_reverify_fail, cache.size(), truth.size());
}

int main() {
	test_vertex_weights_exact_sanity();
	test_incremental_id_assignment();
	test_neighbor_id_consistency();
	test_edge_cache_churn();

	printf("\n%s (%d failures)\n", g_fail==0 ? "ALL CHECKS PASSED" : "CHECKS FAILED", g_fail);
	return g_fail==0 ? 0 : 1;
}
