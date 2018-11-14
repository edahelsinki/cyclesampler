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


typedef struct {
  unsigned n, m, z; 
  unsigned **u, **v; 
  unsigned *nu, *nv;
} BIP;

typedef struct {
  unsigned n;
  double *lambdas;
  unsigned *selfcount;
  double s;
} LAMBDAS;

typedef struct {
  int n;
  int *stack;
  double *val;
  bool *flag;
} GSTACK;

void samplecycles2C(int *,int *,int *,int *,int *,
		    int *,int *,int *,int *,int *,
		    int *,int *,int *,int *,
		    double *,double *,int *,double *);
void samplecyclesC(int *,int *,int *,int *,int *,int *,int *,int *,int *,
		   double *);
void findlambdasC(int *,int *,int *,int *,double *,double *,
		  double *,double *,double *,int *,int *);
void *smalloc(size_t);
int uniform(int);
void initBIP(BIP *,unsigned,unsigned,unsigned **,unsigned *);
void freeBIP(BIP *);
void fswap(double *,double *);
void uswap(unsigned *,unsigned *);
double qinterpolate(double,double,double,
		    double,double,double);
double findroot(double (*)(double,const void *),const void *,
		double,double,double,double,double,int);
double brent(double (*)(double,const void *),const void *,
	     double,double,double,double,double,int,int *);
void permute(unsigned *,unsigned);
double fexpU01(double);
double rootfexpU01(double,const void *);
double findlambda(double);
double compexpaux(double,const LAMBDAS *,double (*)(double));
double compexp(double,const LAMBDAS *);
double findlambdai(unsigned,const BIP *b,const double *,double);
double compexpiaux(unsigned,const BIP *,const double *,double,
		   double (*)(double,const LAMBDAS *));
double compexpi(unsigned,const BIP *,const double *,double);
double optimlambdas(const BIP *,const double *,double *);
double rootcompexp(double,const void *);
double compexp2(double,const LAMBDAS *);
double fexp2U01(double);
double compexp2i(unsigned,const BIP *,const double *,double);
void initgstack(GSTACK *,int);
void freegstack(GSTACK *);
void resetgstack(GSTACK *);
void addgstack(GSTACK *,int,double);
bool allequal(int,int *);
