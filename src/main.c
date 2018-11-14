/*
 * Copyright (c) 2018 Kai Puolamaki <kai.puolamaki@iki.fi>
 * Copyright (c) 2018 University of Helsinki
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <R.h>
#include "main.h"

#define allequalM(n,x) ((x)[0]==(x)[1] && ((n)==2 || ((x)[0]==(x)[2] && (x)[0]==(x)[3])))

void
samplecycles2C(int *iter, /* # of iterations */
	       int *nbft, /* number of nodes in spanning tree */
	       int *ntriplets, /* number of edges (triplets) */
	       int *nidx0, /* number of cycles */
	       int *bft, /* spanning tree, contains parent nodes for
			    each node (except root for which parent
			    node is undefined) */
	       int *depth, /* spanning tree depth for each node */
	       int *st,  /* subtree indices for each node in [0,nst-1] */
	       int *idx, /* triplet index of edges toward root in
			    spanning tree */
	       int *idx0, /* triplet indices of cycles (edges not in
			     spanning tree excluding lone odd
			     cycles) */
	       int *length, /* lengths of cycles corresponding to
			       idx0 */
	       int *nodd, /* array of the number of odd cycles for
			     each subgraph */
	       int *oddsv, /* flattened list of triplet indices of odd
			      cycles per subgraph */
	       int *triplets1, /* triplets[,1] */
	       int *triplets2, /* triplets[,2] */
	       double *a, /* lower bound for the weight of each triplet */
	       double *b, /* upper bound for the weight of each triplet */
	       int *accepted_moves, /* number of accepted moves */
	       double *randomstate) /* the state for each edge (triplet) */
{
  int nst=1;
  int **odds;
  double deltamax=0.0;
  GSTACK s;

  GetRNGstate();
  initgstack(&s,*ntriplets);

  /* Find number of subtrees nst, it is the largest subtree index
     +1. */
  for(int i=0; i<*nbft; i++) if(st[i]>nst-1) nst = st[i]+1;
#ifdef DEBUG
  printf("samplecycles2C: nst = %d\n",nst);
#endif
  
  /* odds[i] is the vector of triplet indices of length nodd[i] of odd
     cycles in subgraph i. */ 
  odds = smalloc(nst*sizeof(int *));
  for(int i=0; i<nst; i++) {
    odds[i] = oddsv;
    oddsv += nodd[i];
  }

  /* find the global maximum delta. */
  for(int i=0; i<*ntriplets; i++)
    if(b[i]-a[i]>deltamax) deltamax = b[i]-a[i];
  
  while(*iter>0) {
    int n=2, i, c1, c2;
    double dmin=-deltamax, dmax=deltamax, delta;
    int ch[4];
    double wt[4];

    resetgstack(&s);
    (*iter)--;
    
    /* first randomly pick a cycle */
    i = uniform(*nidx0);
    c1 = idx0[i]; /* triplet index of the cycle */
    /* start a chain at the nodes of the edge c1 */    
    ch[0] = triplets1[c1];
    ch[1] = triplets2[c1];
    /* next weight for both chains will be -1. */
    wt[0] = -1.0;
    wt[1] = -1.0;
    /* add edge to stack with weight +1. If the edge happens to form a
       self-loop, i.e., ch[0]==ch[1], then the weight should be +2. */
    addgstack(&s,c1,ch[0]==ch[1] ? 2.0 : 1.0); 

#ifdef DEBUG
    printf("samplecycles2C: i = %d c1 = %d ch[0] = %d ch[1] = %d\n",
	   i,c1,ch[0],ch[1]);
#endif


    if(length[i]%2==1) {
      /* if cycle is odd, we have two more paths... */
      n = 4;
      /* pick another odd cycle from the same subgraph. If it so
	 happens that c1==c2 then the join should result to zero
	 change to randomstate. */
      c2 = odds[st[ch[0]]][uniform(nodd[st[ch[0]]])];
#ifdef DEBUG
      printf("samplecycles2C: c2 = %d\n",c2);
#endif
      if(c1==c2) continue; /* Bail out here because the following
			      would be a null operation anyway. */
      ch[2] = triplets1[c2];
      ch[3] = triplets2[c2];
      if((depth[ch[0]]+depth[ch[2]])%2==0) {
	/* If the traversal between the cycles is even then the edge
	   c2 comes with a negative weight. */
	addgstack(&s,c2,ch[2]==ch[3] ? -2.0 : -1.0);
	wt[2] = 1.0;
	wt[3] = 1.0;
      } else {
	/* If the traversal between the cycles is even then the edge
	   c2 comes with a positive weight. */
	addgstack(&s,c2,ch[2]==ch[3] ? 2.0 : 1.0);
	wt[2] = -1.0;
	wt[3] = -1.0;
      }
#ifdef DEBUG
      printf("samplecycles2C: ch[2] = %d ch[3] = %d\n",ch[2],ch[3]);
#endif
    }

    /* Repeat this loop, traversing towards the root, until all chains
       meet at the same node */
    while(!allequalM(n,ch)) {
      /* complex (but a bit more efficient way) to set i to
         arg max depth[ch[i]] */
      int i = depth[ch[0]]>depth[ch[1]] ? 0 : 1;
      for(int j=2; j<n; j++) if(depth[ch[j]]>depth[ch[i]]) i=j;
      /* add it to the stack */
      addgstack(&s,idx[ch[i]],wt[i]);
      /* Proceed in the spanning tree towards the root. */
      ch[i] = bft[ch[i]];
      /* The next weight will be of opposite sign. */
      wt[i] = -wt[i];
    }

    /* Find the allowed values for random perturbation of the matrix
       elements [dmin,dmax] such that all values stay in [a,b]. */
    for(int i=0; i<s.n; i++) {
      int j=s.stack[i];
      double rs=randomstate[j], val=s.val[j], aa=a[j], bb=b[j];
#ifdef DEBUG
      printf("samplecycles2C: i = %d s.stack[i] = %d val = %g\n",
	     i,s.stack[i],val);
#endif
      if(val>0.0) {
	if(bb<rs+dmax*val) dmax = (bb-rs)/val;
	if(rs+dmin*val<aa) dmin = (aa-rs)/val;
      } else {
	if(rs+dmax*val<aa) dmax = (aa-rs)/val;
	if(bb<rs+dmin*val) dmin = (bb-rs)/val;
      }
    }
    /* Find a random perturbation from [dmin,dmax] such that all
       values stay in [a,b]. */
    delta = dmin+(dmax-dmin)*unif_rand();

    /* Increment the acceptance rate counter if delta is non-zero */

    *accepted_moves += (delta != 0);
    
    for(int i=0; i<s.n; i++) {
	int j=s.stack[i];
	randomstate[j] += delta*s.val[j];   
    }
  }

  freegstack(&s);
  free(odds);
  PutRNGstate();
}
/* end of samplecycles2C */

/*
 * allequal()
 *     Tests if all items in array are equal.
 * Return:
 *     true or false.
 */

bool
allequal(int n,int *x) {
  for(int i=1; i<n; i++) if(x[i]!=x[0]) return(false);
  return(true);
}
/* end of allequal */

/*
 * samplecyclesC()
 *      R interface which finds *iter cycles in random using the spanning
 *      tree and modifies the random state using the cycles.
 * Return (by pointer):
 *      New randomstate after *iter random cycles.
 */

void
samplecyclesC(int *iter,
	      int *nbft,
	      int *nidx0,
	      int *bft,
	      int *depth,
	      int *idx,
	      int *idx0,
	      int *triplets1,
	      int *triplets2,
	      double *randomstate)
{
  int *p, *m;
  
  GetRNGstate();

  p = smalloc(*nbft*sizeof(int));
  m = smalloc(*nbft*sizeof(int));

  while(*iter>0) {
    int i, ip=0, im=0, left, right;
    double dmin=-1.0, dmax=1.0, delta;
    bool fleft=false, fright=false;

    /* First randomly pick a link that is not in the spanning tree. */
    i = idx0[uniform(*nidx0)];
    left = triplets1[i];
    right = triplets2[i];
    /* Add this link to the "positives" list. */
    p[ip++] = i;

    /* Start traversing the spanning tree from this link towards the
       root node via the two branches (left and right) that correspond
       to the nodes of the initial link above and add links
       alternatingly to the "positives" and "negatives" lists until
       the left and right branches meet, i.e., we have a cycle. */
    while(left!=right) {
      if(depth[left]>depth[right]) {
	if(fleft) 
	  p[ip++] = idx[left];
	else 
	  m[im++] = idx[left];
	fleft = !fleft;
	left = bft[left];
      } else {
	if(fright)
	  p[ip++] = idx[right];
	else
	  m[im++] = idx[right];
	fright = !fright;
	right = bft[right];
      }
    }

    /* Find minimum and maximum possible values for the change in the
       weights of the links, [dmin,dmax]. */
    for(int j=0; j<ip; j++) {
      if(randomstate[p[j]]+dmax>1.0) dmax = 1.0-randomstate[p[j]];
      if(randomstate[p[j]]+dmin<0.0) dmin = -randomstate[p[j]];
    }
    for(int j=0; j<im; j++) {
      if(randomstate[m[j]]-dmax<0.0) dmax = randomstate[m[j]];
      if(randomstate[m[j]]-dmin>1.0) dmin = randomstate[m[j]]-1.0;
    }

    /* Randomly pick an allowed change from [dmin,dmax]. */
    delta = dmin+(dmax-dmin)*unif_rand();

    /* Finally update the random state adding the change to "positive"
       links and reducing it from the "negative" links. */
    for(int j=0; j<ip; j++) randomstate[p[j]] += delta;
    for(int j=0; j<im; j++) randomstate[m[j]] -= delta;
    
    (*iter)--;
  }

  free(p);
  free(m);
  
  PutRNGstate();
}
/* end of samplecyclesC */

/*
 * findlambdasC()
 *      R interface which finds the values of lambda that correspond
 *      to expectations given in s. n is the number of values and m
 *      is the number of constraints. nui contains the number of edges
 *      associated with each value (always two if we have row and column
 *      sum constraints only) and uv contains indices of the respective
 *      constraints. l contains lambdas. Tolerances are specified by
 *      tol (tolerance of lambda), tole (tolerance of expectations),
 *      and tolz (tolerance of expectations, z-normalized). The
 *      computation of tole and tole may be time consuming.
 * Return (by pointer):
 *      Values of lambdas l, final tolerance values tol, tole, and tolz,
 *      as well as iterations used (maxiter is reduced by one for each
 *      iteration).
 */

void
findlambdasC(int *n, int *m,int *nui,int *uv,double *s,double *l,
	     double *tol,double *tole,double *tolz,int *maxiter,
	     int *initlambda) {
  unsigned z=0;
  double delta, deltae=0.0, deltaz=0.0;
  unsigned *nu, **u, *uu;
  BIP b;

  GetRNGstate();
  
  for(unsigned i=0; i<*n; i++) z += nui[i];
  nu = smalloc(*n*sizeof(unsigned));
  u = smalloc(*n*sizeof(unsigned *));
  uu = smalloc(z*sizeof(unsigned));
  for(unsigned i=0; i<z; i++) uu[i] = uv[i];
  for(unsigned i=0; i<*n; i++) {
    u[i] = uu;
    nu[i] = nui[i];
    uu += nu[i];
  }
  initBIP(&b,*n,*m,u,nu);
  free(nu);
  free(u[0]); /* uu */
  free(u);

  if(*initlambda)
    for(unsigned i=0; i<b.m; i++)
      l[i] = findlambda(s[i]/b.nv[i]);

  do {
    delta = fabs(optimlambdas(&b,s,l));
    if(*tole>=0.0 || *tolz>=0.0) {
      deltae = 0.0;
      deltaz = 0.0;
      for(unsigned i=0; i<b.m; i++) {
	double e1, d1;
	e1 = compexpi(i,&b,l,s[i]);
	d1 = fabs(e1-s[i]);
	if(d1>deltae) deltae = d1;
	if(*tolz>=0.0) {
	  double e2, d2;
	  e2 = compexp2i(i,&b,l,s[i]);
	  assert(e2>0.0);
	  if(e2>0.0) {
	    d2 = d1/sqrt(e2);
	    if(d2>deltaz) deltaz = d2;
	  }
	}
      }
    }
    (*maxiter)--;
  } while(*maxiter && delta>*tol && deltae>*tole && deltaz>*tolz);

  *tol = delta;
  *tole = deltae;
  *tolz = deltaz;
    
  freeBIP(&b);

  PutRNGstate();
}
/* end of findlambdasC */


/*
 * smalloc()
 *      malloc that bails out if out of memory
 * Return:
 *      pointer to a memory block
 */

void *
smalloc(size_t n) {
  void *p;
  if((p = malloc(n))==NULL) {
    fprintf(stderr,
	    "smalloc: malloc failure (n=%lu).\n",
	    (unsigned long)n);
    exit(EXIT_FAILURE);
  }
  return(p);
}
/* end of smalloc */

/*
 * uniform()
 *      Produces a random integer uniformly from [0,n-1].
 * Return:
 *      Random integer in [0,n-1].
 */

int
uniform(int n) {
  int i;
  /* unif_rand() should never be one, but just to make sure... */
  do { i = (int)floor(n*unif_rand()); } while(i>=n);
  return(i);
}
/* end of uniform */

/*
 * We represent bipartite graph with structure BIP. The bipartite graph
 * consists of left-hand nodes and right-hand nodes. There are undirected
 * links between some of the left and right handed nodes.
 * There are n left-hand nodes and m the right-hand nodes. In our case
 * the left-hand nodes correspond to values (e.g., movie grades)
 * and right hand-nodes to  constraints (in our case one constraint 
 * for each row and each column). z is the number of edges in the graph
 * (in our case m==2*n, because each value is associated with a row
 * constraint and a column constraint). nu is an array of n unsigned 
 * integers and nv is an array of m unsigned integerss that contains the
 * degrees of the left-hand and the right-hand nodes, respectively. E.g., 
 * nu[i] is the degree of the left-hand node indexed by i. u and v are
 * pointers to arrays of unsigned integers that contain the indices of 
 * nodes linked to the respective node. E.g., u[i] is a array of length 
 * nu[i] that contains the indices of the right-hand nodes which are 
 * linked to the left-hand node i. 
 *
 *
 *
 * initBIP()
 *      initializes structure for a bipartite graph ("BIP").
 * Return (via pointer):
 *      Initialized data structure which should later be free with
 *      freeBIP().
 */

void
initBIP(BIP *b,
	unsigned n,
	unsigned m,
	unsigned **u,
	unsigned *nu) {
  unsigned *uu, *vv;
  
  b->n = n;
  b->m = m;
  b->z = 0;
  for(unsigned i=0; i<n; i++) b->z += nu[i];
  b->u = smalloc(n*sizeof(unsigned *));
  b->v = smalloc(m*sizeof(unsigned *));
  b->nu = smalloc(n*sizeof(unsigned));
  b->nv = smalloc(m*sizeof(unsigned));
  uu = smalloc(b->z*sizeof(unsigned));
  vv = smalloc(b->z*sizeof(unsigned));

  for(unsigned j=0; j<m; j++) {
    b->nv[j] = 0;
  }

  for(unsigned i=0; i<n; i++) {
    b->nu[i] = nu[i];
    b->u[i] = uu;
    for(unsigned k=0; k<nu[i]; k++) {
      *(uu++) = u[i][k];
      b->nv[u[i][k]]++;
    }
  }

  for(unsigned j=0; j<m; j++) {
    b->v[j] = vv;
    vv += b->nv[j];
    b->nv[j] = 0;
  }

  for(unsigned i=0; i<n; i++) 
    for(unsigned k=0; k<nu[i]; k++) 
      b->v[u[i][k]][b->nv[u[i][k]]++] = i;
}
/* end of initBIP */

/*
 * freeBIP()
 *      frees structure for a bipartite graph ("BIP"), should be
 *      called after initBIP() when BIP structure is no longer needed.
 * Return:
 *      None.
 */

void
freeBIP(BIP *b) {
  free(b->u[0]); /* uu */
  free(b->v[0]); /* vv */
  free(b->u);
  free(b->v);
  free(b->nu);
  free(b->nv);
}
/* end of freeBIP */

/* findlambda()
 *      Finds the value of l corresponding to a expectation E[x], inverse 
 *      of fexpU01().
 * Return:
 *      Value of l such that fexpU01(l)=x.
 */

double
findlambda(double x) {
  static double xmax=-1.0, tol;
  if(xmax<0.0) {
    xmax = pow(DBL_MAX,0.25);
    tol = pow(DBL_EPSILON,0.25);
  }
  return(findroot(&rootfexpU01,&x,0.0,-xmax,xmax,1.0,tol,10000));
}
/* end of findlambda */

/*
 * findroot()
 *      Intelligent root finder that takes as input a function for
 *      which root is sought, starting point, and the maximal interval.
 *      The root finder first tries to find a smallest interval in which
 *      the function to be tested changes sign and then uses a standard
 *      root finding methods (now Brent's). We assume here that f is 
 *      increasing and continuous function. This triest to be 
 *      numerically stable.
 * Return:
 *      Location of the root x where approximately f(x)=0, or xmin or xmax
 *      if the root is not in the interval [xmin,xmax].
 */

double
findroot(double (*f)(double,const void *),
		const void *p,
		double x0,
		double xmin,
		double xmax,
		double step,
		double tol,
		int maxiter) {
  double a, b, fa, fb, fx0, gr=(1.0+sqrt(5.0))/2.0;
  int status;

  /* Do nothing if we are already at the solution. */
  a = x0-tol;
  b = x0+tol;
  fa = f(a,p);
  fb = f(b,p);
  if(fa<=0.0 && fb>=0.0) {
    return(x0);
  }
  if(fa<0.0) {
    x0 = b;
    fx0 = fb;
  } else {
    x0 = a;
    fx0 = fa;
  }

  /* If we do not cross the zero then return the closest boundary as a
     solution. */
  fa = f(xmin,p);
  fb = f(xmax,p);
  if(fa>=0.0 && fb>=0.0) {
    return(xmin);
  }
  if(fa<=0.0 && fb<=0.0) {
    return(xmax);
  }

  /* Check that f is increasing function. */
  if(fa>0.0) {
    fprintf(stderr,
	    "findroot: f should be increasing function f(%g) = %g.\n",
	    xmin,fa);
    exit(EXIT_FAILURE);
  }

  /* Find interval xmin, xmax that contains the root. */
  a = x0;
  fa = fx0;
  if(fa<0.0) {
    b = a+step;
    fb = f(b,p);
#ifdef DEBUG
    printf("findroot: [ f(%g) = %g, f(%g) = %g ]\n",a,fa,b,fb);
#endif
    while(fb<0.0) {
      a = b;
      fa = fb;
      b = x0+gr*(b-x0);
      if(b>xmax) { b = xmax; }
      fb = f(b,p);
    }
  } else {
    b = a;
    fb = fa;
    a = b-step;
    fa = f(a,p);
#ifdef DEBUG
    printf("findroot: [ f(%g) = %g, f(%g) = %g ]\n",a,fa,b,fb);
#endif
    while(fa>0.0) {
      b = a;
      fb = fa;
      a = x0-gr*(x0-a);
      if(a<xmin) { a = xmin; }
      fa = f(a,p);
    }
  }

#ifdef DEBUG
  printf("findroot: done [ f(%g) = %g, f(%g) = %g ]\n",a,fa,b,fb);
#endif
  
  return(brent(f,p,a,b,fa,fb,tol,maxiter,&status));
}
/* end of findroot */

/*
 * qinterpolate()
 *      Given three values a, b, c and f(a), f(b), f(c), use inverse
 *      quadratic interpolation to guess location x such that f(x)=0.
 *      The function smartly falls back to the secant method if 
 *      inverse quadratic interpolation is not possible.
 * Return:
 *      Location x such that f(x) might be zero.
 */

double
qinterpolate(double a,double b,double c,
	     double fa,double fb,double fc) {
  double fab, fac, fbc;
  
  if(fb==fa) {
    fswap(&b,&c);
    fswap(&fb,&fc);
  }

#ifdef DEBUG
  if(fb==fa) {
    fprintf(stderr,
	    "qinterpolate: a = %g  b = %g  c = %g  fa = fb = fc = %g\n",
	    a,c,b,fa);
    exit(EXIT_FAILURE);
  }
#endif
  assert(fb!=fa);

  fab = fa-fb;
  if(fa==fc || fb==fc) {
    /* secant method */
    return(-a*fb/fab+b*fa/fab);
  }

  fac = fa-fc;
  fbc = fb-fc;
  /* inverse quadratic interpolation */
  return(a*fb*fc/(fab*fac)-b*fc*fa/(fbc*fab)+c*fa*fb/(fac*fbc));
}
/* end of qinterpolate */

/*
 * brent()
 *      Brent's method to find a root of a function, as described in
 *      https://en.wikipedia.org/wiki/Brent's_method
 * Return:
 *      The location of the root x where f(x)=0 is approximately
 *      satisfied. Furthermore, parameter *status is set to the number
 *      of iterations performed.
 */

double
brent(double (*f)(double,const void *),
      const void *p,
      double a,
      double b,
      double fa,
      double fb,
      double tol,
      int maxiter,
      int *status) {
  double c, d, s, fc, fs, tol1, x1, x2, eps;
  bool mflag=true; /* true if previous step used bisection method */
#ifdef DEBUG
  int count=0;
#endif

  *status = 0;
  
  /* compute machine epsilon */
  eps = DBL_EPSILON;
  
  /* check that fa and fb have different signs */
  if((fa<0.0 && fb<0.0) || (fa>0.0 && fb>0.0)) {
    fputs("brent: fa and fb do not have different signs.\n",stderr);
    exit(EXIT_FAILURE);
  }

  if(fabs(fa)<fabs(fb)) {
    fswap(&a,&b);
    fswap(&fa,&fb);
  }
  c = a;
  fc = fa;

  tol1 = 4.0*eps*fabs(b)+tol;
  while(fabs(a-b)>tol1 && fb!=0.0) {
#ifdef DEBUG
    printf("\nbrent: iteration %d, a = %g  fa = %g  b = %g  "
	   "fb = %g  c = %g  fc = %g\n",
	   ++count,a,fa,b,fb,c,fc);
#endif
    
    if(!(maxiter--)) {
      fputs("brent: maxiter reached.\n",stderr);
      break;
    }

    /* use quadratic or linear interpolation to guess the root */
    s = qinterpolate(a,b,c,fa,fb,fc);
    x1 = (3.0*a+b)/4.0;
    x2 = b;
    if(x2<x1) { fswap(&x1,&x2); }
    if((s<=x1 || x2<=s) || /* condition 1: s is outside (x1,x2) */
       (mflag && fabs(s-b)>=fabs(b-c)/2.0) || /* condition 2 */
       (!mflag && fabs(s-b)>=fabs(c-d)/2.0) || /* condition 3 */
       (mflag && fabs(b-c)<tol1) || /* condition 4 */
       (!mflag && fabs(c-d)<tol1) /* condition 5 */ ) {
      /* bisection method is forced */
      s = (a+b)/2.0;
      mflag = true;
#ifdef DEBUG
      printf("brent: bisection forced, ");
#endif
    } else {
      /* interpolation is accepted */
      mflag = false;
#ifdef DEBUG
      printf("brent: quadratic / linear interpolation, ");
#endif
    }
    fs = f(s,p);
#ifdef DEBUG
    printf("s = %g  f(s) = %g\n",s,fs);
#endif
    d = c; /* d is the previous previous best guess. d is assigned for
	      the first time here; it won't be used above on the first
	      iteration because mflag is set. */
    c = b; /* c is the previous best guess */
    fc = fb;
    if((fa<0.0 && fs>0.0) || (fa>0.0 && fs<0.0)) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if(fabs(fa)<fabs(fb)) {
      fswap(&a,&b);
      fswap(&fa,&fb);
    }
    tol1 = 4.0*eps*fabs(b)+tol;
    (*status)++;
  }

#ifdef DEBUG
  printf("brent: done, b = %g\n",b);
#endif
  
  return(b);
}
/* end of brent */

/*
 * fexpU01()
 *      Computes of E[x]-*y0 under distribution \propto exp(l*x)
 *      where x is in [0,1]. Can be used with root finding algorithm
 *      to find inverse function of fexpU01().
 * Return:
 *      The value of E[x]-*y0 for a given l.
 */

double
fexpU01(double l) {
  double t;
  static double eps = -1.0;
  if(eps<0.0) eps = pow(DBL_EPSILON,0.25);
  if(fabs(l)<eps) return(0.5);
  t = exp(l)-1.0;
  return(1.0+1.0/t-1.0/l);
}
/* end of fexpU01 */


double
rootfexpU01(double x,const void *y0) {
  return(fexpU01(x)-*(double *)y0);
}
/* end of rootfexpU01 */

/*
 * optimlambdas()
 *      Optimize all constraints in random order using findlambdai().
 * Return:
 *      The maximal change in parameters. We may have convergence when
 *      this value is sufficiently small. Modifies the values of lamdba
 *      stored in array l.
 */

double
optimlambdas(const BIP *b,const double *s,double *l) {
  double delta=0.0;
  unsigned *idx;

  idx = smalloc(b->m*sizeof(unsigned));
  for(unsigned i=0; i<b->m; i++) idx[i] = i;
  permute(idx,b->m);

  for(unsigned i=0; i<b->m; i++) {
    unsigned j = idx[i];
    double old = l[j];
    l[j] = findlambdai(j,b,l,s[j]);
    if(fabs(l[j]-old)>fabs(delta)) delta = l[j]-old;
  }

  free(idx);

  return(delta);
}
/* end of optimlambdas */

/*
 * permute()
 *      Permutes array of n unsigned integers uniformly in random using
 *      Knuth's shuffle.
 * Return (via pointer):
 *      Array p is permuted uniformly in random.
 */

void
permute(unsigned *p,unsigned n) {
  /* Knuth shuffle */
  for(unsigned i=0; i<n-1; i++) uswap(&p[i],&p[i+uniform(n-i)]);
}
/* end of permute */

/*
 * fswap()
 *      Swaps the values of two doubles
 */

void
fswap(double *a,double *b) {
  double c=*a;
  *a = *b;
  *b = c;
}
/* end of fswap */


/*
 * uswap()
 *      Swaps the values of two unsigned integers.
 */

void
uswap(unsigned *a,unsigned *b) {
  unsigned c=*a;
  *a = *b;
  *b = c;
}
/* end of uswap */


/*
 * findlambdai()
 *      Find lambda related to the i'th constraint such that the 
 *      expectation equals s (typically sum of values related to the
 *      constraint).
 * Return:
 *      Value of lambda such that the expectation of the constraint 
 *      equals s.
 */

double
findlambdai(unsigned i,const BIP *b,const double *l,double s) {
  double x;
  LAMBDAS a;
  static double xmax=-1.0, tol;

  if(xmax<0.0) {
    xmax = pow(DBL_MAX,0.25);
    tol = pow(DBL_EPSILON,0.25);
  }

  a.n = b->nv[i];
  a.selfcount = smalloc(a.n*sizeof(unsigned));
  a.lambdas = smalloc(a.n*sizeof(double));
  a.s = s;
  
  for(unsigned j=0; j<a.n; j++) {
    a.selfcount[j] = 0;
    a.lambdas[j] = 0.0;
    for(unsigned k=0; k<b->nu[b->v[i][j]]; k++) {
      unsigned ii = b->u[b->v[i][j]][k];
      if(ii==i)
	a.selfcount[j]++; 
      else
	a.lambdas[j] += l[ii];
    }
  }

  x = findroot(&rootcompexp,&a,l[i],-xmax,xmax,1.0,tol,10000);
  
  free(a.selfcount);
  free(a.lambdas);

  return(x);
}
/* end of findlambdai */

/*
 * rootcompexp()
 *      Root version of compexp(), can be used to construct
 *      inverse of compexp().
 * Return:
 *      Expectation minus a given value a->s.
 */

double
rootcompexp(double lambda,const void *p) {
  const LAMBDAS *a=p;
  return(compexp(lambda,p)-a->s);
}
/* end of rootcompexp */

/*
 * compexpaux()
 *      Auxiliary function to compute the expectations for constraints.
 * Return:
 *      Expectation of f.
 */

double
compexpaux(double lambda,const LAMBDAS *a,double (*f)(double)) {
  double y = 0.0;
  for(unsigned i=0; i<a->n; i++)
    y += f(a->selfcount[i]*lambda+a->lambdas[i]);
  return(y);
}
/* end of compexpaux */

/*
 * compexp()
 *      Function to compute the expectations for constraints.
 * Return:
 *      Expectation of fexpU01.
 */

double
compexp(double lambda,const LAMBDAS *a) {
  return(compexpaux(lambda,a,&fexpU01));
}
/* end of compexp */

/*
 * compexpiaux()
 *      Function to compute expectation for the constraint. i is the 
 *      index of the constraint in the bipartite tree b (where the
 *      right-hand nodes correspond to the constraints), l is the vector
 *      of lambdas corresponding to the constraints, s is the target value
 *      for the constraint (e.g., sum computed from the actual data),
 *      and f is the function whose root gives the constraint. We use
 *      structure LAMBDAS as a parameter to f. The idea is that the 
 *      constraint expectation function for which we want to find a root
 *      is f(x)=sum(Ej[x])-a.s, where the expectation j is over a 
 *      probability distribution of form 
 *      \propto exp(a.selfcount[j]*lambda+a.lambdas[j]).
 * Return:
 *      Returns the expectation of the constraint i.
 */

double
compexpiaux(unsigned i,const BIP *b,const double *l,double s,
	    double (*f)(double,const LAMBDAS *)) {
  LAMBDAS a;

  a.n = b->nv[i];
  a.selfcount = smalloc(a.n*sizeof(unsigned));
  a.lambdas = smalloc(a.n*sizeof(double));
  a.s = s;
  
  for(unsigned j=0; j<a.n; j++) {
    a.selfcount[j] = 0;
    a.lambdas[j] = 0.0;
    for(unsigned k=0; k<b->nu[b->v[i][j]]; k++) {
      unsigned ii = b->u[b->v[i][j]][k];
      if(ii==i)
	a.selfcount[j]++; 
      else
	a.lambdas[j] += l[ii];
    }
  }

  s = f(l[i],&a);

  free(a.selfcount);
  free(a.lambdas);

  return(s);
}
/* end of compexpiaux */


/*
 * compexpi()
 *      Interface to compexpiaux() (see above) using compexp().
 * Return:
 *      Expectation related to constraint i.
 */

double
compexpi(unsigned i,const BIP *b,const double *l,double s) {
  return(compexpiaux(i,b,l,s,&compexp));
}
/* end of compexpi */


/*
 * compexp2()
 *      Function to compute the expectations for constraints.
 * Return:
 *      Expectation of fexp2U01.
 */

double
compexp2(double lambda,const LAMBDAS *a) {
  return(compexpaux(lambda,a,&fexp2U01));
}
/* end of compexp2 */

/*
 * compexpi2()
 *      Computes variance related to constraint i.
 * Return:
 *      Variance of the value of i.
 */

double
compexp2i(unsigned i,const BIP *b,const double *l,double s) {
  return(compexpiaux(i,b,l,s,&compexp2));
}
/* end of compexpi2 */

/*
 * fexp2U01()
 *      Computes expectation of E[x^2]-E[x]^2 (variance) under distribution 
 *      \propto exp(l*x) where x is in [0,1].
 * Return:
 *      The value of variance for a given l.
 */

double
fexp2U01(double l) {
  double t;
  static double eps = -1.0;
  if(eps<0.0) eps = pow(DBL_EPSILON,0.25);
  if(fabs(l)<eps) return(1.0/12.0);
  t = exp(l)-1.0;
  return(1.0/(l*l)-1.0/t-1.0/(t*t));
}
/* end of fexp2U01 */


/*
 * initgstack()
 *     Initializes gstack data structure.
 * Return:
 *     Nothing.
 */


void
initgstack(GSTACK *s,int size) {
  s->n = 0;
  s->stack = smalloc(size*sizeof(int));
  s->val = smalloc(size*sizeof(double));
  s->flag = smalloc(size*sizeof(bool));
  for(int i=0; i<size; i++) s->flag[i] = false;
}
/* end of initstack */

/* 
 * freegstack()
 *     Frees gstack data structure.
 * Return:
 *     Nothing.
 */

void freegstack(GSTACK *s) {
  free(s->val);
  free(s->flag);
  free(s->stack);
}
/* end of freegstack */

/* 
 * resetgstack()
 *     Resets gstack data structure to initial state.
 * Return:
 *     Nothing.
 */

void
resetgstack(GSTACK *s) {
  while(s->n>0) s->flag[s->stack[--s->n]] = false;
}
/* end of resetgstack */

/*
 * addgstack()
 *     Adds index i with value val to gstack data structure.
 * Return:
 *     Nothing.
 */

void
addgstack(GSTACK *s,int i,double val) {
  if(s->flag[i]) {
    s->val[i] += val;
  } else {
    s->flag[i] = true;
    s->val[i] = val;
    s->stack[s->n] = i;
    s->n++;
  }
}
/* end of addgstack */
 
