// lf.c -- brute-force leaky forcing checker (standard zero forcing colour change rule, leaks cannot force).
// Build:  gcc -O2 -o lf lf.c
// Leaks: exactly min(ell,n) leaks are placed adversarially (more leaks is never better for the forcer),
// and every leak set of size <= ell is checked.
// Input (stdin): n m ell mode ; then m edges (u v, 0-indexed); if mode==1: s then s vertices (set to test)
// mode 0: compute Z^(ell) (n<=64) by brute force (supersets of low-degree vertices, justified by Lemma 2.2)
// mode 1: test whether given set is ell-leaky forcing (n<=1000, ell<=2); prints OK or FAIL + first failing leak set
// Output mode 0: the number Z^(ell)(G).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
int n,m,ell,mode;
u64 adj[64];
int *nb[1000]; int deg[1000];

static int closure_mask(u64 blue, u64 leak){
  u64 all = (n==64)?~0ULL:((1ULL<<n)-1);
  int changed=1;
  while(changed){
    changed=0;
    u64 cand = blue & ~leak;
    while(cand){
      int v=__builtin_ctzll(cand); cand&=cand-1;
      u64 w = adj[v] & ~blue;
      if(w && !(w&(w-1))){ blue|=w; changed=1; }
    }
  }
  return blue==all;
}
// check all leak sets of size L (recursive)
static int leakok(u64 blue,int L,int start,u64 leak){
  if(L==0) return closure_mask(blue,leak);
  for(int v=start; v<=n-L; v++){
    if(!leakok(blue,L-1,v+1,leak|(1ULL<<v))) return 0;
  }
  return 1;
}
int isleaky(u64 S){
  if(!closure_mask(S,0)) return 0;
  int L = ell<n?ell:n;
  for(int j=1;j<=L;j++) if(!leakok(S,j,0,0)) return 0; // monotone, but cheap pruning order
  return 1;
}
int freev[64],nf; u64 base;
u64 found;
int rec(int k,int start,u64 S){
  if(k==0){ if(isleaky(S)){found=S;return 1;} return 0; }
  for(int i=start;i<=nf-k;i++) if(rec(k-1,i+1,S|(1ULL<<freev[i]))) return 1;
  return 0;
}
// big-graph closure
int *color,*cnt,*isleak,*queue;
int closure_big(int *S,int s,int l1,int l2){
  for(int v=0;v<n;v++){color[v]=0;isleak[v]=0;}
  if(l1>=0) isleak[l1]=1;
  if(l2>=0) isleak[l2]=1;
  for(int i=0;i<s;i++)color[S[i]]=1;
  int nblue=0; for(int v=0;v<n;v++)nblue+=color[v];
  for(int v=0;v<n;v++){cnt[v]=0; for(int j=0;j<deg[v];j++) if(!color[nb[v][j]])cnt[v]++;}
  int qh=0,qt=0;
  for(int v=0;v<n;v++) if(color[v]&&!isleak[v]&&cnt[v]==1) queue[qt++]=v;
  while(qh<qt){
    int v=queue[qh++];
    if(isleak[v]||cnt[v]!=1) continue;
    int w=-1; for(int j=0;j<deg[v];j++) if(!color[nb[v][j]]){w=nb[v][j];break;}
    if(w<0) continue;
    color[w]=1; nblue++;
    for(int j=0;j<deg[w];j++){int u=nb[w][j]; cnt[u]--; if(color[u]&&!isleak[u]&&cnt[u]==1) queue[qt++]=u;}
    if(!isleak[w]&&cnt[w]==1) queue[qt++]=w;
  }
  return nblue==n;
}
int main(){
  if(scanf("%d %d %d %d",&n,&m,&ell,&mode)!=4) return 1;
  int *eu=malloc(sizeof(int)*m),*ev=malloc(sizeof(int)*m);
  for(int i=0;i<m;i++){if(scanf("%d %d",&eu[i],&ev[i])!=2) return 1;}
  if(mode==0){
    memset(adj,0,sizeof adj);
    for(int i=0;i<m;i++){adj[eu[i]]|=1ULL<<ev[i]; adj[ev[i]]|=1ULL<<eu[i];}
    base=0;nf=0;
    for(int v=0;v<n;v++){ if(__builtin_popcountll(adj[v])<=ell) base|=1ULL<<v; else freev[nf++]=v; }
    for(int k=0;k<=nf;k++){ if(rec(k,0,base)){ printf("%d\n",__builtin_popcountll(found)); return 0;} }
    printf("%d\n",n); return 0;
  } else {
    for(int v=0;v<n;v++){deg[v]=0;}
    for(int i=0;i<m;i++){deg[eu[i]]++;deg[ev[i]]++;}
    for(int v=0;v<n;v++){nb[v]=malloc(sizeof(int)*deg[v]);deg[v]=0;}
    for(int i=0;i<m;i++){nb[eu[i]][deg[eu[i]]++]=ev[i]; nb[ev[i]][deg[ev[i]]++]=eu[i];}
    int s; if(scanf("%d",&s)!=1) return 1; int *S=malloc(sizeof(int)*s); for(int i=0;i<s;i++) if(scanf("%d",&S[i])!=1) return 1;
    color=malloc(4*n);cnt=malloc(4*n);isleak=malloc(4*n);queue=malloc(sizeof(int)*(2*n+2*m+16)); /* pushes <= initial n + one per edge-endpoint decrement + one per newly blue vertex */
    if(!closure_big(S,s,-1,-1)){printf("FAIL none\n");return 0;}
    if(ell>=1) for(int a=0;a<n;a++) if(!closure_big(S,s,a,-1)){printf("FAIL %d\n",a);return 0;}
    if(ell>=2) for(int a=0;a<n;a++) for(int b=a+1;b<n;b++) if(!closure_big(S,s,a,b)){printf("FAIL %d %d\n",a,b);return 0;}
    printf("OK\n"); return 0;
  }
}
