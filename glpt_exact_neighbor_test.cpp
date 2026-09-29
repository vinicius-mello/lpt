// Build: g++ -O2 -std=c++17 -I. glpt_exact_neighbor_test.cpp lpt.cpp -o glpt_exact_neighbor_test && ./glpt_exact_neighbor_test 1000000
// Exact (integer) verification of glpt::neighbor(): for every sampled
// (cell, facet i), the returned neighbour must be a different cell at the
// same level whose facet i is the SAME set of points, where each point is
// written exactly as integer barycentric weights over the 15 global
// Gaifullin vertices (common denominator 2^(MAX_ORTHANT_LEVEL+1)).
#include <cstdio>
#include <cstdlib>
#include <random>
#include <array>
#include <algorithm>
#include "glpt.hpp"

const int D=glpt::DIM;
typedef std::array<long,15> P;

void points(const glpt& c, P out[D+1]) {
	int w[D+1][D+1], sh; c.vertex_weights_exact(w,sh);
	int scale = (glpt::MAX_ORTHANT_LEVEL+1) - sh;
	for(int k=0;k<=D;++k) {
		out[k].fill(0);
		for(int j=0;j<=D;++j) out[k][glpt_gaifullin_cells[c.seed_index()][j]] += (long)w[k][j] << scale;
	}
}

int main(int argc,char**argv) {
	long N = argc>1 ? atol(argv[1]) : 200000;
	glpt_gaifullin_init();
	std::mt19937_64 rng(2026);
	int maxdepth = glpt::MAX_ORTHANT_LEVEL*D + D-1;
	long cnt[3]={0,0,0}, bad[3]={0,0,0}, badidx[3]={0,0,0};
	const char* name[3]={"sibling","interior (Theorem 2)","seed crossing (same code)"};
	for(long n=0;n<N;++n) {
		glpt c((int)(rng()%GLPT_NCELLS));
		int depth=rng()%(maxdepth+1);
		for(int s=0;s<depth;++s) c=c.child((int)(rng()&1));
		P pc[D+1]; points(c,pc);
		for(int i=0;i<=D;++i) {
			int kind;
			if(c.is_root()) kind=2;
			else {
				bool c0=c.is_child0(); int ls=(c.level()==0)?D:c.level();
				if((c0&&i==D)||(!c0&&i==ls)) kind=0;
				else kind = (c.root_boundary_facet_fast(i)<0) ? 1 : 2;
			}
			glpt r; if(c.neighbor(i,r)!=glpt::OK) { ++bad[kind]; ++cnt[kind]; continue; }
			++cnt[kind];
			P pr[D+1]; points(r,pr);
			std::vector<P> fc, fr;
			for(int k=0;k<=D;++k) if(k!=i) fc.push_back(pc[k]);
			std::sort(fc.begin(),fc.end());
			// does r contain the facet, with the excluded vertex not in it?
			std::vector<P> all(pr,pr+D+1); std::sort(all.begin(),all.end());
			bool contains=std::includes(all.begin(),all.end(),fc.begin(),fc.end());
			bool same_cell=true; { std::vector<P> a(pc,pc+D+1); std::sort(a.begin(),a.end()); same_cell=(a==all); }
			if(!contains || same_cell || r.simplex_level()!=c.simplex_level()) { ++bad[kind]; continue; }
			for(int k=0;k<=D;++k) if(k!=i) fr.push_back(pr[k]);
			std::sort(fr.begin(),fr.end());
			if(fr!=fc) ++badidx[kind]; // facet sits at a different local index in r
		}
	}
	for(int k=0;k<3;++k) printf("%-28s checked=%ld wrong=%ld different-local-index=%ld\n",name[k],cnt[k],bad[k],badidx[k]);
}
