#ifndef _RANDOMH_
#define _RANDOMH_

#include <cstdlib>
#include <cstdio>

#include "hmcerrs.h"
#include "globaldefs.h"
#include "types.h"
#include "operations.h"
#include "geometry.h"

typedef unsigned long long int Ullong;
typedef unsigned int Uint;

//Seed for Random
const unsigned long long int seed = 500000;

//mod. NR-version

struct Ran {
	Ullong u,v,w;
	Ran(){}
	Ran(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}

	// init function (steo)
	void Init(Ullong j) {
	v=4101842887655102017LL;
	w=1;
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}

	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};

/*
Aus Numerical recipes -- black box
*/


struct Random {
	Ullong u,v,w;
	Random(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};

/*
Aus Numerical Recipes. Angeblich genauso gut, aber viel schneller -- könnte man ja mal testen
*/

struct Ranq1 {
	Ullong v;
	Ranq1(Ullong j) : v(4101842887655102017LL) {
		v ^= j;
		v = int64();
	}
	inline Ullong int64() {
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};

/*
Aus Numerical Recipes. Exponentielle Verteilung -- p(x)~ b*exp(bx)
*/

struct Expondev : Random {
	double beta;
	Expondev(double bbeta, Ullong i) : Random(i), beta(bbeta) {}
	double dev() {
		double u;
		do u = doub(); while (u == 0.);
		return -log(u)/beta;
	}
};

//Zufallszahl vom Typ Double zwischen 0 und 1
double random_01 ();

//Zufallszahl vom Typ Double aus ]0,1]
double random_01_halboffen ();

//Zufallszahl 1,2,3 vom Typ int 
int random_123 ();

//Zufallszahl 0,1,2,3 vom Typ int
int random_0123 ();

//Zufallszahl zwischen -1 und 1
double random_11 ();

//!! to go
//Zufälliger Gitterpunkt auf dem Gitter
//multi1d<int> random_gitter(const multi1d<int> & latt_size);

//Zufallszahl 0,1 vom Typ int 
int random_0_1 ();

//Gibt Zufallszahl zwischen 0 und 2*Pi
double random_02pi ();

//Gibt drei Zufallszahlen 1,2,3
void random_1_2_3 (int rand[3]);


#endif
