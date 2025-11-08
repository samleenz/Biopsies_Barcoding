/*******************************************************************************
 tumour_barcoding_simulation.cpp — adapted from TumourSimulator v1.2.2
 
 This program is a minor adaptation of TumourSimulator, originally developed by
 Dr Bartek Waclaw (University of Edinburgh), to support the tracking of genetic
 barcodes during tumour growth simulations.
 
 All credits for the original design, structure, and simulation mechanics go to
 Dr Waclaw and collaborators.
 
 Adaptation and barcode integration by:
 T.S. Weber (weber.ts@wehi.edu.au)
 For:
 Serrano & Weber et al.
 "Genetic barcoding uncovers the clonal makeup of solid and liquid biopsies
 and their ability to capture intra-tumoral heterogeneity" 
 
 This adaptation preserves the GPLv3 licensing terms of the original code.
 A copy of the GNU General Public License can be found at:
 <http://www.gnu.org/licenses/>
 
 Original reference:
 Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban,
 Bert Vogelstein, and Martin A. Nowak. "Spatial Model Predicts That 
 Dispersal and Cell Turnover Limit Intratumour Heterogeneity"
 Nature 525, no. 7568 (2015): 261–264. doi:10.1038/nature14971
 
 ------------------------------------------------------------------------------
 
 Compilation:
 Linux/macOS: g++ tumour_barcoding_simulation.cpp -w -O3 -o tumour_barcoding_simulation
 Windows (MinGW-w64): g++ tumour_barcoding_simulation.cpp -w -O3 -lpsapi -o tumour_barcoding_simulation.exe
 
 Execution:
 ./tumour_barcoding_simulation <output_folder> <n_samples> <random_seed>
 
 Example:
 ./tumour_barcoding_simulation output 1 1
 
 *******************************************************************************/
 
 //Key Parameters (adjust as needed):
 int max_size     = int(1e5);   // population size after which simulation stops
 int ini_barcodes = 200;        // number of initial barcoded cells
 

/*******************************************************************************
   TumourSimulator v.1.2.2 - a program that simulates a growing solid tumour.
   Based on the algorithm described in
   
   Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban, 
   Bert Vogelstein, and Martin A. Nowak. "Spatial Model Predicts That 
   Dispersal and Cell Turnover Limit Intratumour Heterogeneity" Nature 525, 
   no. 7568 (September 10, 2015): 261?64. doi:10.1038/nature14971.

   Contributing author:
   Dr Bartek Waclaw, University of Edinburgh, bwaclaw@staffmail.ed.ac.uk

   Copyright (2015) The University of Edinburgh.

    This file is part of TumourSimulator.

    TumourSimulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TumourSimulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file 
    License.txt or at <http://www.gnu.org/licenses/>.
*******************************************************************************/

// Previous versions:

// v1.2.2 - optional pushing added
// v1.2.1 - saves most_abundant
// v1.2.0 - new modes added, files reorganised
// v1.1.0 - new modes added
// v1.0.1 - drivers not saved anymore
// v1.0

#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp
#include <vector>
#include <iostream>
#include <random>
using namespace std;

char *NUM ; // name given as 1st argument from the command line

// ===== BEGIN INLINED HEADER: params.h =====
#line 1 "params.h"
/*******************************************************************************
   TumourSimulator v.1.2.2 - a program that simulates a growing solid tumour.
   Based on the algorithm described in
   
   Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban, 
   Bert Vogelstein, and Martin A. Nowak. "Spatial Model Predicts That 
   Dispersal and Cell Turnover Limit Intratumour Heterogeneity" Nature 525, 
   no. 7568 (September 10, 2015): 261-64. doi:10.1038/nature14971.

   Contributing author:
   Dr Bartek Waclaw, University of Edinburgh, bwaclaw@staffmail.ed.ac.uk

   Copyright (2015) The University of Edinburgh.

    This file is part of TumourSimulator.

    TumourSimulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TumourSimulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file 
    License.txt or at <http://www.gnu.org/licenses/>.
*******************************************************************************/

//---- exactly one method should be defined------
//#define GILLESPIE
//#define FASTER_KMC
#define NORMAL  // standard simulation method described in the paper (non-KMC)
//-----------------------------------------------

//#define PUSHING // if defined, cells can push away other cells as they grow

//#define CONST_BIRTH_RATE // if defined, birth rate does not depend on the no. of empty neighbours as long as there is at least one

// ----------exactly one of these should be defined as the replication neighbourhood -----------
//#define VON_NEUMANN_NEIGHBOURHOOD // 6 neighbours
//#define VON_NEUMANN_NEIGHBOURHOOD_QUADRATIC // 6 neighbours but non linear growth rate (12n-n^2)/36 where n = no. of empty neighbours
#define MOORE_NEIGHBOURHOOD // 26 neighbours
//--------------------------------------------------------


//#define MAKE_TREATMENT_N // if defined, simulate treatment after reaching given size
//#define MAKE_TREATMENT_T // if defined, simulate treatment after reaching given time

const float gama=1e-2, gama_res=5e-8 ; // these are rates per daughter cell. Rates per diploid exome will be 2x higher (these values are given in the paper)

//#define MIGRATION_MATRIX

// used when MIGRATION_MATRIX is not defined
const float migr=0 ;
// used only when MIGRATION_MATRIX is defined
//const float migr[2][2]={{0,0} // before treatment: WT/resistant
//                        ,{0,1e-5}}; // after treatment: WT/resistant

//#define DEATH_ON_SURFACE ; // if defined then cells die on surface only upon treatment, if not then cells die also in the volume
//#define CORE_IS_DEAD  // when set, core cells are set to dead 
#define SHOW_ONLY_DRIVERS
//#define NO_MECHANICS  // if defined, no mechanics is simulated (this speeds up everything but looses info on spatial distribution)
const float timescale=1./log(2.) ; // calculates the timescale factor from the division time [days] (first number, here 1.)
const float death0=0.95, growth0=1.0 ;   // before treatment

// if death on surface:
//float death1=/*0.1*/0.99, growth1=0.0 ;   // after treatment
// if death in volume
const float death1=1., growth1=0.5 ; // after treatment

const float driver_adv=0.05 ; 
const float driver_prob=2e-5 ; // driver probability per haploid genome (should be 2e-5)
const float driver_balance=1 ; // 1==drivers decrease death, 0==drivers increase growth, intermediate values are also possible
const float driver_migr_adv=0.0,  max_migr=1e-5 ; // maximum migration prob. is taken into account only when driver_migr_adv>0
const int driver_mode = 0 ; // 0== drivers affect bd only, 1==drivers affect simultaneously bd and migr, 2==drivers affect bd xor migr with equal prob.

const float cutoff=0.1 ;
#ifdef __MAIN
float time_to_treat=10*(365./12) ; // this is ignored when MAKE_TREATMENT_N is defined
#endif

const int _resol=1 ; // spatial resolution of sampling [cells]
const int _bins=10000 ; // max number of bins
//#define PAUSE_WHEN_MEMORY_LOW 10000 ; // if defined, the program pauses if there is less than PAUSE_WHEN_MEMORY_LOW MB of free memory

// ===== END INLINED HEADER: params.h =====
#line 1 "simulation.cpp"

// ===== BEGIN INLINED HEADER: classes.h =====
#line 1 "classes.h"
/*******************************************************************************
   TumourSimulator v.1.2.2 - a program that simulates a growing solid tumour.
   Based on the algorithm described in
   
   Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban, 
   Bert Vogelstein, and Martin A. Nowak. "Spatial Model Predicts That 
   Dispersal and Cell Turnover Limit Intratumour Heterogeneity" Nature 525, 
   no. 7568 (September 10, 2015): 261-64. doi:10.1038/nature14971.

   Contributing author:
   Dr Bartek Waclaw, University of Edinburgh, bwaclaw@staffmail.ed.ac.uk

   Copyright (2015) The University of Edinburgh.

    This file is part of TumourSimulator.

    TumourSimulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TumourSimulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file 
    License.txt or at <http://www.gnu.org/licenses/>.
*******************************************************************************/


#include <stdio.h>
#include <math.h>
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(x) (x)*(x)
#define SWAPD(x, y) tempd = (x); (x) = (y); (y) = tempd
#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp

#include <vector>
using namespace std;

double _drand48(void) ;
void _srand48(int a) ;
void err(char *reason) ;
void err(char *reason, int a) ;
void err(char *reason, double a);
void quicksort2(float *n, int *nums, int lower, int upper) ;
void init();
int main_proc(int exit_size, int save_size, double max_time, double wait_time) ;
void end() ;
void reset() ;
void save_data() ;
void save_spatial(int *) ;
void save_snps(char *name,int *n, int total, int mode, int *) ;
float average_distance_ij() ;

extern int L ;
extern const int _resol, _bins ;
extern int volume, max_size ; 

#ifndef classes_already_defined
#define classes_already_defined

#define __linux

#ifdef __linux
typedef unsigned int DWORD ;
typedef unsigned char BYTE ;
#else 
//#include <windows.h>
#endif

class vecd  // class of 3d vectors
{
public:
	double x, y, z;

	vecd(double x0, double y0, double z0) : x(x0), y(y0), z(z0) { }
	vecd() : x(0), y(0), z(0) {	}
	inline bool operator== (const vecd& V2) const  { return (x == V2.x && y == V2.y && z == V2.z) ; }
	inline vecd operator+ (const vecd& V2) const { return vecd(x+V2.x,y+V2.y,z+V2.z) ; }
	inline vecd operator- (const vecd& V2) const { return vecd(x-V2.x,y-V2.y,z-V2.z) ; }
	inline vecd operator- () const { return vecd(-x,-y,-z);	}
	inline vecd operator/ (double S ) const {
		double f = 1.0/S;
		return vecd (x*f,y*f,z*f);
	}
	inline vecd operator* (const vecd& V2) const { return vecd (x*V2.x,y*V2.y,z*V2.z); }
	inline vecd operator* (double S) const { return vecd (x*S,y*S,z*S);	}
	inline void operator+= ( const vecd& V2 ) {
		x+=V2.x; y+=V2.y; z+=V2.z;
	}
	inline void operator-= ( const vecd& V2 ) {
		x -= V2.x; y -= V2.y; z -= V2.z;
	}
	inline vecd operator*= (double S) {
		x*=S ; y*=S ; z*=S ;
	}
	inline vecd operator/= (double S) {
    double S2=1./S ;
    x*=S2 ; y*=S2 ; z*=S2 ;
	}
};

inline double norm(vecd &a)
{
  return sqrt(a.x*a.x+a.y*a.y+a.z*a.z) ;
}

inline double squared(vecd &a)
{
  return (a.x*a.x+a.y*a.y+a.z*a.z) ;
}

inline double scalar(vecd &a, vecd &b) {
  return a.x*b.x+a.y*b.y+a.z*b.z ; 
}

inline vecd cross(vecd &a,vecd &b)
{
  vecd c ; 
  c.x=a.y*b.z-a.z*b.y ;
  c.y=a.z*b.x-a.x*b.z ;
  c.z=a.x*b.y-a.y*b.x ;
  return c ;  
}

inline void normalize(vecd &a)
{
  double l=sqrt(a.x*a.x+a.y*a.y+a.z*a.z) ;
  a.x/=l ; a.y/=l ; a.z/=l ;  
}

class IVec  // class of 3d integer vectors
{
  public:
	int i, j, k;
	IVec( int Ini, int Inj, int Ink ) : i( Ini ), j( Inj ), k( Ink ) { }
	IVec( ) : i(0), j(0), k(0) { }
	inline bool operator== (const IVec& V2) const { return (i == V2.i && j == V2.j && k == V2.k); }
	inline void operator+= ( const IVec& V2 ) { i += V2.i; j += V2.j; k += V2.k; }
};


struct Cell {
  short unsigned int lesion ;
  short int x,y,z ;
  unsigned int gen ; 
};

#ifndef PUSHING
class Sites {
  public:
    DWORD *s ;
    Sites(int n0) { int n=1+(n0>>5) ; s=new DWORD[n] ; if (s==NULL) err("out of memory when allocating Sites") ; for (int i=0;i<n;i++) s[i]=0 ; }
    ~Sites() { delete [] s ; }
    inline void set(const unsigned int i) { s[(i>>5)]|=1<<(i&31) ; }
    inline void unset(const unsigned int i) { s[(i>>5)]&=~(1<<(i&31)) ; }
    inline int is_set(const unsigned int i) { return (s[(i>>5)]>>(i&31))&1 ; }
};
#else
typedef int Sites ;
#endif

extern vector <Cell> cells ;

#ifndef PUSHING
struct Lesion {
  int wx ; 
  vecd r,rold,rinit ; 
  double rad, rad0 ;
  int n,n0 ; 
  vector <int> closest ;
  Sites **p ;
  static int nl ;
  static double maxdisp ;
  Lesion(int cellno, int g, int x0, int y0, int z0) {
    // cellno is not used anywhere in this version of the method
    rad=rad0=1 ; 
    r=vecd(x0,y0,z0) ; rinit=rold=r ;
    closest.clear() ;
    wx=4 ; p=new Sites*[wx*wx] ;
    int i,j,k;
    for (i=0;i<wx*wx;i++) {
      p[i]=new Sites(wx) ;
    }
    Cell c ; c.x=c.y=c.z=0 ; c.gen=g ; c.lesion=nl++ ; if (nl>65000) err ("nl>65000") ;
    p[(wx/2)*wx+wx/2]->set(wx/2) ;
    cells.push_back(c) ; volume++ ; n=n0=1 ; 
  }  
  ~Lesion() {
    nl-- ; 
    for (int i=0;i<wx*wx;i++) delete p[i] ;
    delete [] p ;    
  }
  void update_wx() ;
  void find_closest() ;
  void one_move_step() ;  
  void reduce_overlap() ; 
  int no_free_sites(int x, int y, int z);  
  void choose_nn(int &x, int &y, int &z);  
};
#else
struct Lesion {
  int wx ; 
  vecd r,rold,rinit ; 
  double rad, rad0 ;
  int n,n0 ; 
  vector <int> closest ;
  Sites ***p ;
  static int nl ;
  static double maxdisp ;
  Lesion(int cellno, int g, int x0, int y0, int z0) {
    rad=rad0=1 ; 
    r=vecd(x0,y0,z0) ; rinit=rold=r ;
    closest.clear() ;
    wx=16 ; p=new Sites**[wx] ;
    int i,j,k;
    for (i=0;i<wx;i++) {
      p[i]=new Sites*[wx] ; if (p[i]==NULL) err("out of memory") ;
      for (j=0;j<wx;j++) {
        p[i][j]=new Sites[wx] ;
        for (k=0;k<wx;k++) p[i][j][k]=-1 ;
      }
    }
    Cell c ; c.x=c.y=c.z=0 ; c.gen=g ; c.lesion=nl++ ; if (nl>65000) err ("nl>65000") ;
    p[wx/2][wx/2][wx/2]=cellno ;
    cells.push_back(c) ; volume++ ; n=n0=1 ; 
  }  
  ~Lesion() {
    nl-- ; 
    for (int i=0;i<wx;i++) {
      for (int j=0;j<wx;j++) delete p[i][j] ;
      delete [] p[i] ;
    }
    delete [] p ;    
  }
  void update_wx() ;
  void find_closest() ;
  void one_move_step() ;  
  void reduce_overlap() ; 
  int no_free_sites(int x, int y, int z);  
  void find_dir_min_drag(int i, int j, int k, int &in, int &jn, int &kn);       // find direction of least drag  
};
#endif

const unsigned int RESISTANT_PM = 1<<31 ;
const unsigned int DRIVER_PM = 1<<30 ;
const unsigned int L_PM = (1<<30) - 1 ; 

struct Genotype {
  vector <unsigned int> sequence ;
  BYTE no_resistant, no_drivers ;
  float death[2], growth[2], m[2] ; // m = migration probability before/after treatment
  int number ; // number of cells of this genotype total/on the surface
  int prev_gen ;
  int index ; // this is set only when saving data
  Genotype(void) ;
  ~Genotype(void) { sequence.clear() ; }
  Genotype(Genotype *mother, int prevg, int no_snp) ;
};

class Hist {
  public:
  int x,n ;
  long long int x2 ;
  Hist() { x=n=0 ; x2=0 ; }
	inline void operator+= ( Hist& H2 ) { x+=H2.x ; x2+=H2.x2 ; n+=H2.n ; }
	inline void operator+= ( int val ) { x+=val ; x2+=val*val ; n++ ; }
  void r() { x=n=0 ; x2=0 ; }
};


#endif

extern vector<Genotype*> genotypes ;
extern vector<Lesion*> lesions ;

// ===== END INLINED HEADER: classes.h =====
#line 1 "simulation.cpp"

#if defined __linux
#include <unistd.h>
typedef unsigned int DWORD ;
int memory_taken() // return memory available in MB
{
  long long int rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	long long int mem=((size_t)rss * (size_t)sysconf( _SC_PAGESIZE)) ;
	return (int) (mem/(1<<20));
}
#include <sys/sysinfo.h>
unsigned int freemem() // returns available memory in MB
{
  struct sysinfo s ;
  sysinfo(&s) ;
  return ((s.freeram)>>20) ;
}
#elif defined __APPLE__
typedef unsigned int DWORD ;
int memory_taken()
{
  return 0 ; // not implemented
}
#else
#include <windows.h>
#include <psapi.h>
int memory_taken() // return memory available in MB
{
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (int) (info.WorkingSetSize/(1<<20));
}
#endif

void err(char *reason)
{
  cout <<reason<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif  
  exit(0) ;
}

void err(char *reason, int a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, char *a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

void err(char *reason, double a)
{
  cout <<reason<<": "<<a<<endl ; 
#ifdef __WIN32
  system("pause") ;
#endif    
  exit(0) ;
}

static long long unsigned int _x=0x000100010001LL, _mul=0x0005deece66dLL, _add=0xbLL ;
double _drand48(void)  // works only on compilers with long long int!
{
  _x=_mul*_x+_add ; _x&=0xffffffffffffLL ;
  return (_x/281474976710656.0) ;
}

void _srand48(int a) { _x=a ; }

void init();
void end() ;

double tt=0, tt_at_start ;
int start_clock ;

int L=0 ; // total number of SNPs
int volume ; // total volume of the tumor
vector <int> drivers ; // vector of driver mutations
FILE *drivers_file ;
int treatment=0, cells_at_start ;
FILE *times ; 
extern int sample ;
int RAND ; // random seed
char *timesbuffer ;

int poisson(void)  // generates k from P(k)=exp(-gamma) gamma^k / k!
{
  return 0;
  /*const double l=exp(-gama) ;
  double p=1. ;
  int k=0 ;
  do {
    k++ ;
    p*=_drand48() ;
  } while (p > l) ;
  return k - 1 ;*/
}

vector<Cell> cells ;

int Lesion::nl=0 ;
double Lesion::maxdisp=0 ;
double max_growth_rate ;

Genotype::Genotype(void) 
{ 
  death[0]=death0 ; death[1]=death1 ; growth[0]=growth0 ; growth[1]=growth1 ;
#ifdef MIGRATION_MATRIX
  m[0]=migr[0][0] ; m[1]=migr[1][0] ;
#else
  m[0]=m[1]=migr ; 
#endif
  number=1 ; no_resistant=no_drivers=0 ; sequence.clear() ; prev_gen=-1 ;
}

Genotype::Genotype(Genotype *mother, int prevg, int no_snp) { 
  death[0]=mother->death[0] ; growth[0]=mother->growth[0] ; m[0]=mother->m[0] ;
  death[1]=mother->death[1] ; growth[1]=mother->growth[1] ; m[1]=mother->m[1] ;
  prev_gen=prevg ;
  sequence=mother->sequence ; no_resistant=mother->no_resistant ; no_drivers=mother->no_drivers; 
  for (int i=0;i<no_snp;i++) {
    if ((driver_adv>0 || driver_migr_adv>0) && _drand48()<driver_prob/gama) { 
      float q=_drand48() ;
      if (driver_mode<2 || q<0.5) {
        death[0]*=1-driver_adv*driver_balance ; 
        growth[0]*=1+driver_adv*(1-driver_balance) ; if (max_growth_rate<growth[0]) max_growth_rate=growth[0] ;
      }
      if (driver_migr_adv>0 && ((q>=0.5 && driver_mode==2) || driver_mode==1)) {
        m[0]*=1+driver_migr_adv ; if (m[0]>max_migr) m[0]=max_migr ;
        m[1]*=1+driver_migr_adv ; if (m[1]>max_migr) m[1]=max_migr ;
      }
      // drivers decrease prob. of death or increase prob. of growth
      drivers.push_back(L) ; //fprintf(drivers_file,"%d ",L) ; fflush(drivers_file) ; 
      sequence.push_back((L++)|DRIVER_PM) ; no_drivers++ ;
    } else {
      if (_drand48()<gama_res/gama) {  
        sequence.push_back((L++)|RESISTANT_PM) ; no_resistant++ ; // resistant mutation
        death[1]=death0 ; growth[1]=growth0 ; 
#ifdef MIGRATION_MATRIX
        m[0]=migr[0][1] ; m[1]=migr[1][1] ;
#endif
      } 
      else sequence.push_back(L++) ;
    }
  }
  if (L>1e9) err("L too big") ;
  number=1 ;
}

vector<Genotype*> genotypes ;
vector<Lesion*> lesions ;

void Lesion::update_wx()
{
  int i,j,k;
  int nwx=int(wx*1.25) ;
  if (nwx%2==1) nwx++ ; // make sure it's even
  int dwx=(nwx-wx)/2 ;

#ifdef PUSHING
  for (i=0;i<wx;i++) {
    for (j=0;j<wx;j++) {
      Sites *np=new Sites[nwx] ;
      for (k=0;k<nwx;k++) np[k]=-1 ;
      for (k=dwx;k<wx+dwx;k++) np[k]=p[i][j][k-dwx] ;
      delete p[i][j] ;
      p[i][j]=np ;
    }
  }

  Sites ***np=new Sites**[nwx] ;
  for (i=0;i<nwx;i++) np[i]=new Sites*[nwx] ;
    
  for (i=0;i<nwx;i++) {
    for (j=0;j<nwx;j++) {
      if (i<dwx || i>=wx+dwx || j<dwx || j>=wx+dwx) {
        np[i][j]=new Sites[nwx] ;
        for (k=0;k<nwx;k++) np[i][j][k]=-1 ;
      } else {
        np[i][j]=p[i-dwx][j-dwx] ; 
      }
    }
  }

  for (i=0;i<wx;i++) delete [] p[i] ; 

#else
  for (i=0;i<wx*wx;i++) {
    Sites *np=new Sites(nwx) ;
    for (k=dwx;k<wx+dwx;k++) if (p[i]->is_set(k-dwx)) np->set(k) ;
    delete p[i] ;
    p[i]=np ;
  }

  if (float(nwx)*float(nwx)>2e9) err("nwx too large",nwx) ;
  Sites **np=new Sites*[nwx*nwx] ;
    
  for (i=0;i<nwx;i++) {
    for (j=0;j<nwx;j++) {
      if (i<dwx || i>=wx+dwx || j<dwx || j>=wx+dwx) {
        np[i*nwx+j]=new Sites(nwx) ;
      } else {
        np[i*nwx+j]=p[(i-dwx)*wx+j-dwx] ; 
      }
    }
  }
#endif

  delete [] p ;
  p=np ;
  wx=nwx ;

}

void Lesion::one_move_step() {
  int i,j;
  double mthis=this->n ;
  for (i=0;i<closest.size();i++) {
    vecd dr=lesions[closest[i]]->r - this->r ;
    double r2=squared(dr), sumrad2=SQR(this->rad+lesions[closest[i]]->rad) ;
    if (r2<sumrad2) {
      double mi=lesions[closest[i]]->n ;
      double disp=(sqrt(sumrad2/r2)-1) ;
      if (fabs(disp)>maxdisp) maxdisp=fabs(disp) ;
      dr*=disp*1.1 ;
      this->r-=dr*mi/(mi+mthis) ;
      lesions[closest[i]]->r+=dr*mthis/(mi+mthis) ;
    }    
  }
}

void Lesion::find_closest() 
{
  rold=r ;
  closest.clear() ;
  for (int i=0;i<lesions.size();i++) {
    vecd dr=this->r - lesions[i]->r ;
    double r2=squared(dr) ;
    if (r2>0 && r2<2*(SQR(this->rad+lesions[i]->rad))) {
      closest.push_back(i) ;
    }
  }  
}

void Lesion::reduce_overlap()
{
  int i,j,k,temp ;
  int *ind=new int[lesions.size()] ;
  for (j=0;j<lesions.size();j++) ind[j]=j ;
  do {
    maxdisp=0 ;
    for (j=0;j<lesions.size();j++) { k=_drand48()*lesions.size() ; SWAP(ind[j],ind[k]) ; }
    for (j=0;j<lesions.size();j++) {  // go through a random permutation
      i=ind[j] ; 
      lesions[i]->one_move_step() ; 
        
      vecd dr=lesions[i]->r - lesions[i]->rold ; 
      if (squared(dr)>SQR(lesions[i]->rad)) lesions[i]->find_closest() ; 
    }    
  } while (maxdisp>1e-2) ;  
  delete [] ind ;
}  

void reset() 
{
  tt=0 ; L=0 ; max_growth_rate=growth0 ;
  treatment=0 ; 
  for (int i=0;i<genotypes.size();i++) if (genotypes[i]!=NULL) delete genotypes[i] ;
  genotypes.clear() ; genotypes.push_back(new Genotype) ;  
  for (int i=0;i<lesions.size();i++) delete lesions[i] ;
  lesions.clear() ;
  cells.clear() ; volume=0 ;
  drivers.clear() ;
  lesions.push_back(new Lesion(0,0, 0,0,0)) ;
  
  // erase output buffer for "times"
#if defined __linux
  times->_IO_write_ptr = times->_IO_write_base ;
#elif defined __APPLE__
  // not defined yet
#else
  times->_ptr = times->_base ; // this operates on elements of _iobuf and is specific to Windows GNU C++
#endif
}

#ifdef MOORE_NEIGHBOURHOOD
const int _nonn=26 ;
const int kx[27]={0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1,0,1,1,0,-1,-1,-1,0,1},
          ky[27]={0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1,0,0,1,1,1,0,-1,-1,-1},
          kz[27]={0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1};
int kln[27] ; // this is filled with lengths of (kx,ky,kz)
#endif

#ifdef VON_NEUMANN_NEIGHBOURHOOD
const int _nonn=6 ;
const int kx[7]={0,1,-1,0,0,0,0},
          ky[7]={0,0,0,1,-1,0,0},
          kz[7]={0,0,0,0,0,1,-1};
int kln[7] ; // this is filled with lengths of (kx,ky,kz)
#endif



void init()
{
  int i,j,k;
  for (i=0;i<=_nonn;i++) kln[i]=sqrt(1.*SQR(kx[i])+1.*SQR(ky[i])+1.*SQR(kz[i])) ;

  char txt[256] ;
  sprintf(txt,"mkdir %s",NUM) ; system(txt) ;
  sprintf(txt,"%s/%s_%d.dat",NUM,NUM,RAND) ; times=fopen(txt,"w") ;
  timesbuffer=new char[(1<<16)] ;
  setvbuf (times , timesbuffer , _IOFBF , (1<<16));  // this is to prevent saving data if no fflush is attempted 
                                                  // (this e.g. allows one to discard N<256)
  start_clock=clock() ;
}

void end() {
  fclose(times) ; 
}

#ifdef PUSHING
inline int Lesion::no_free_sites(int x, int y, int z)
{
  int nfree=_nonn ;
  for (int n=1;n<=_nonn;n++) nfree-=p[(wx+z+kz[n])%wx][(wx+y+ky[n])%wx][(wx+x+kx[n])%wx]==-1?0:1 ;
  return nfree ;
}
#else
inline int Lesion::no_free_sites(int x, int y, int z)
{
  int nfree=_nonn ;
  for (int n=1;n<=_nonn;n++)
    nfree-=p[((wx+z+kz[n])%wx)*wx + (wx+y+ky[n])%wx]->is_set((wx+x+kx[n])%wx) ;
  return nfree ;
}
inline void Lesion::choose_nn(int &x, int &y, int &z)
{
  static int nns[_nonn] ;
  int no=0,n ;
  for (n=1;n<=_nonn;n++)
    if (p[((wx+z+kz[n])%wx)*wx + (wx+y+ky[n])%wx]->is_set((wx+x+kx[n])%wx)==0) nns[no++]=n ;
  if (no==0) { x=-1000000 ; return ; }
  n=nns[int(_drand48()*no)] ; 
  z=(wx+z+kz[n])%wx ; y=(wx+y+ky[n])%wx ; x=(wx+x+kx[n])%wx ;
}
#endif


inline int free_sites(int n)
{
  Lesion *ll=lesions[cells[n].lesion] ;
  int wx=ll->wx ; 
  int k=cells[n].x+wx/2, j=cells[n].y+wx/2, i=cells[n].z+wx/2 ; 
  return ll->no_free_sites(k,j,i) ;  
}


void quicksort2(float *n, int *nums, int lower, int upper)
{
	int i, m, temp ;
  float pivot, tempd;
	
	if (lower < upper)
	{
		SWAPD(n[lower], n[(upper + lower) / 2]); SWAP(nums[lower], nums[(upper + lower) / 2]);
		pivot = n[lower];
		m = lower;
		for (i = lower + 1; i <= upper; i++)
			if (n[i] > pivot)
			{
				m++;
				SWAPD(n[m], n[i]); SWAP(nums[m], nums[i]);
			}
		SWAPD(n[lower], n[m]); SWAP(nums[lower], nums[m]);
		quicksort2(n, nums, lower, m - 1);
		quicksort2(n, nums, m + 1, upper);
	}
}


void save_data()
{
  int i,j,ntot=cells.size(), nsurf=0 ;
  double raver=0, raver2=0 ;
  int no_on_surface=0 ;
  int no_resistant=0, no_resistant_surf=0 ;
  int cells_drv=0,cells_drv_surf=0 ;
  double drv_per_cell=0,drv_per_cell_surf=0, pms_per_cell=0 ;
  double aver_growth_rate=0,av_migr=0 ;

  int *snp_no=new int[L] ; // array of SNPs abundances
  for (i=0;i<L;i++) { snp_no[i]=0 ; }
  for (i=0;i<genotypes.size();i++) {
    if (genotypes[i]!=NULL && genotypes[i]->number>0) {
      for (int j=0;j<genotypes[i]->sequence.size();j++) snp_no[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;      
    }
  }  
  int snps_det=0 ;
  for (i=0;i<L;i++) if (snp_no[i]>cutoff*ntot) snps_det++ ;
  delete [] snp_no ;

  for (i=0;i<ntot;i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    int wx=ll->wx ; 
    double rr=SQR(cells[i].x+ll->r.x)+SQR(cells[i].y+ll->r.x)+SQR(cells[i].z+ll->r.x) ;
    raver+=sqrt(rr) ; raver2+=rr ;
#ifdef PUSHING
    if (ll->p[cells[i].z+wx/2][cells[i].y+wx/2][cells[i].x+wx/2]==-1) err("p[][][]=",i) ;
#else
    if (ll->p[(cells[i].z+wx/2)*wx+cells[i].y+wx/2]->is_set(cells[i].x+wx/2)==0) err("save: is_set: ",i) ;
#endif

    Genotype *g=genotypes[cells[i].gen] ; if (g==NULL) err("g=NULL)") ;
    int free_sites=ll->no_free_sites(cells[i].x+wx/2,cells[i].y+wx/2,cells[i].z+wx/2) ;
    int is_on_surface=(free_sites>0?1:0) ;    
    pms_per_cell+=g->sequence.size() ;

    if (g->no_resistant) {
      no_resistant++ ; 
      if (is_on_surface) no_resistant_surf++ ;
    }
    if (g->no_drivers>0) {
      cells_drv++ ; drv_per_cell+=g->no_drivers ; 
      if (is_on_surface) { cells_drv_surf++ ; drv_per_cell_surf+=g->no_drivers ; }
    }
    if (is_on_surface) nsurf++ ;    
    aver_growth_rate+=g->growth[treatment]*free_sites/float(_nonn) ;
    av_migr+=g->m[treatment] ;
  }
  raver/=ntot ; raver2/=ntot ; aver_growth_rate/=timescale ;
  drv_per_cell/=ntot ; drv_per_cell_surf/=nsurf ; pms_per_cell/=ntot ;
  av_migr/=ntot ;

  // 1.ntot 2.time  3.#genotypes  4.radius      
#ifndef CORE_IS_DEAD
  fprintf(times,"%d %lf %d %lf  ",ntot,tt,genotypes.size(),raver) ; //sqrt(raver2-raver*raver)) ;
#else
  fprintf(times,"%d %lf %d %lf  ",volume,tt,genotypes.size(),raver) ; //sqrt(raver2-raver*raver)) ;
#endif
  //  5.#cells_surf    6.#metas     7.#resistant  8.#resistant_surf
  fprintf(times,"%d %d %d %d   ",nsurf,lesions.size(),no_resistant,no_resistant_surf) ;
  // 9.#drivers   10.#cells_with_drv  11.#cells_with_drv_surf    12.#drv/cell   13.#der/cell_surf
  fprintf(times,"%d %d %d  %lf %lf  ",drivers.size(),cells_drv,cells_drv_surf,drv_per_cell,drv_per_cell_surf) ;
  // 14.growth_rate(n)   15.av_distance   16.pms_per_cell   17.snps_detected  18.<migr>
  fprintf(times,"%lf\t%f\t%lf %d\t %le\t",aver_growth_rate,average_distance_ij(),pms_per_cell,snps_det,av_migr) ;
  // #MBs   time_taken
  fprintf(times,"%d %f\n",memory_taken(),float(1.*(clock()-start_clock)/CLOCKS_PER_SEC)) ;
  if (treatment>0 || ntot>512 || ntot==max_size) fflush(times) ; // flush only when size big enough, this allows us to discard runs that died out

  //if (ntot>256) { printf("%d %lf   no.les.=%d  no.res=%d drv_cell=%lf max_growth=%lf\n",ntot,tt,lesions.size(),no_resistant, drv_per_cell,max_growth_rate) ; fflush(stdout) ; }
}

void snps_corr(Hist *snps) ;
void snps_corr_cutoff(Hist *snps,float cutoff,int *snp_no) ;
void snps_corr_cond_driver(Hist *snps) ;
void find_p_driver(Hist *pr1, Hist *pr2, Hist *pr3) ;
void save_snp_corr(char *name, Hist *snps);

void save_spatial(int *snp_no)
{
#ifndef NO_MECHANICS  
  printf("save spatial\n") ;
  char tmp[256] ;
  Hist *snp_corr, *snp_corr_cutoff ;  // this is for measuring correlations between PMs in different parts of the tumor
  Hist *p_driver1, *p_driver2, *p_driver3, *snp_corr_cd ;

  snp_corr=new Hist[_bins] ; 
  snp_corr_cutoff=new Hist[_bins] ; 
  snp_corr_cd=new Hist[_bins] ; 

  snps_corr(snp_corr) ;
  sprintf(tmp,"%s/corr_%d_%d.dat",NUM,RAND,sample) ; save_snp_corr(tmp, snp_corr) ;                
  snps_corr_cutoff(snp_corr,0.1,snp_no) ;
  sprintf(tmp,"%s/cutoff01corr_%d_%d.dat",NUM,RAND,sample) ; save_snp_corr(tmp, snp_corr) ;                

  if (driver_adv>0 || driver_migr_adv>0) {
    p_driver1=new Hist[_bins] ; p_driver2=new Hist[_bins]; p_driver3=new Hist[_bins]; 
    find_p_driver(p_driver1,p_driver2,p_driver3) ;
    sprintf(tmp,"%s/P_driver1_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, p_driver1) ;
    sprintf(tmp,"%s/P_driver2_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, p_driver2) ;
    sprintf(tmp,"%s/P_driver3_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, p_driver3) ;
    snps_corr_cond_driver(snp_corr_cd) ;
    sprintf(tmp,"%s/corr_cond_driver_%d_%d.dat",NUM,RAND,sample) ;  save_snp_corr(tmp, snp_corr_cd) ;
    delete [] p_driver1 ;    delete [] p_driver2 ;    delete [] p_driver3 ;  
  }

  delete [] snp_corr ;  delete [] snp_corr_cutoff ; delete [] snp_corr_cd ; 

  printf("done\n") ;
#endif
}

#ifdef PUSHING
void Lesion::find_dir_min_drag(int i, int j, int k, int &in, int &jn, int &kn)       // find direction of least drag
{
  int nn, in0,jn0,kn0 ;
  vector <IVec> vis ; // vector of visited sites. it will begin with (i,j,k) and end at empty site

rep:      
  vis.clear() ;
  vis.push_back(IVec(i,j,k)) ;
  in0=i ; jn0=j ; kn0=k ; 
      
  do {      // loop goes over subsequent pushing events
    float mind=wx ;
    nn=-1 ;
    for (int nnnn=0;nnnn<10;nnnn++) {
      int nnn=1+_drand48()*_nonn ;
      in=in0 ; jn=jn0 ; kn=kn0 ; 
      for (float drag=0;drag<mind;drag+=kln[nnn]) {
        in+=kz[nnn] ; jn+=ky[nnn] ; kn+=kx[nnn] ;
        for (int vv=0;vv<vis.size();vv++) if (vis[vv]==IVec(in,jn,kn)) goto brk ; // reject if trajectory passes through prev. visited sites
        if (p[in][jn][kn]==-1) { mind=drag ; nn=nnn ; break ; } 
      }
brk:  continue ;
    }
    if (nn==-1) goto rep ; 
        // now nn gives the direction of pushing
        
    in0+=kz[nn] ; jn0+=ky[nn] ; kn0+=kx[nn] ; // update position of the cell to be pushed
    vis.push_back(IVec(in0,jn0,kn0)) ; // and remember it...
  } while (p[in0][jn0][kn0]!=-1) ; // if the next position contains an empty site then exit

  // push all remembered cells except mother to make space for a single new daughter cell
  Sites sup=-1 ;        // sup is the elevated cell that needs to be inserted into new position
  for (int vv=1;vv<vis.size();vv++) {
    in0=vis[vv].i ; jn0=vis[vv].j ; kn0=vis[vv].k ;  
    Sites snew=p[in0][jn0][kn0] ;  
    p[in0][jn0][kn0]=sup ;
    if (sup!=-1) {
      Cell *c=&cells[sup] ; 
      c->x=kn0-wx/2 ; c->y=jn0-wx/2 ; c->z=in0-wx/2 ;
    }
    sup=snew ;
  }

  in=vis[1].i ; jn=vis[1].j ; kn=vis[1].k ; // the new cell will be the second position from the list (1st is the mother cell which is not pushed)
  if (p[in][jn][kn]!=-1) err("!!!") ;
}
#endif


//-----------------------------------------------------
#if defined(NORMAL)

int main_proc(int exit_size, int save_size, double max_time, double wait_time)
{
  int i,j=0,k,n,l,in,jn,kn,ntot;  
  int cc=0, timeout=0 ;
  double tt_old=tt ;
  
  std::random_device g;
  std::mt19937_64 r(g());
  std::uniform_int_distribution<int>I(100,300);
  //int ini_barcodes=200;//I(r);
  
  double next_time=0; int counter=0;
  
  for(int K=0;;K++) {      // main loop 
  #ifdef PAUSE_WHEN_MEMORY_LOW
    timeout++ ; if (timeout>1000000) {
      timeout=0 ; 
      while (freemem()<PAUSE_WHEN_MEMORY_LOW) { sleep(1) ; } 
    }    
#endif
    double tsc=0.01*cells.size() ; if (tsc>1./max_growth_rate) tsc=1./max_growth_rate ;
    tt+=tsc*timescale/cells.size() ; 
    n=_drand48()*cells.size() ;
    
    if(K%1000000==0) {cout<<K<<" "<<"pop size "<<cells.size()<<endl;}
    
    if(cells.size()==ini_barcodes)
    {
      for(int c=0; c<cells.size(); c++)
      {
        cells[c].gen=c;
        genotypes.push_back(new Genotype(genotypes[0],c,0)) ; //first entry is link to mother cell (link to first cell for now)
      }
      
      for (int i=0;i<genotypes.size();i++) {
          if (genotypes[i]!=NULL && genotypes[i]->number>0) genotypes[i]->index=j++ ;
         }
    }

    Lesion *ll=lesions[cells[n].lesion] ;
    int wx=ll->wx ; 
    k=cells[n].x+wx/2 ; j=cells[n].y+wx/2 ; i=cells[n].z+wx/2 ; 
    int need_wx_update=0 ;
    if (k<2 || k>=ll->wx-3 || j<2 || j>=ll->wx-3 || i<2 || i>=ll->wx-3) need_wx_update=1 ; 
#ifdef PUSHING
    if (ll->p[i][j][k]!=n) err("ll->p[i][j][k]!=n, p=",ll->p[i][j][k]) ;
#else
    if (ll->p[i*wx+j]->is_set(k)==0) err("ll->p[i][j][k]==0, wx=",wx) ;
#endif

    if (_drand48()<tsc*genotypes[cells[n].gen]->growth[treatment]) { // reproduction
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
      int nn=1+int(_drand48()*_nonn) ;
      int in=(wx+i+kz[nn])%wx, jn=(wx+j+ky[nn])%wx, kn=(wx+k+kx[nn])%wx ;
#ifdef VON_NEUMANN_NEIGHBOURHOOD_QUADRATIC // one more trial to find an empty site
      if (ll->p[in*wx+jn]->is_set(kn)==1) {
        nn=1+int(_drand48()*_nonn) ;
        in=(wx+i+kz[nn])%wx ; jn=(wx+j+ky[nn])%wx ; kn=(wx+k+kx[nn])%wx ;
      }
#endif
      if (ll->p[in*wx+jn]->is_set(kn)==0) {
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
      int in=i, jn=j, kn=k ;
      ll->choose_nn(kn,jn,in) ;
      if (kn!=-1000000) { // if there is at least one empty n.n., then.....
#elif defined(PUSHING)
      int in,jn,kn ;
      ll->find_dir_min_drag(i,j,k, in,jn,kn) ;
      {
#else
  #error inconsistent growth conditions 
#endif
        int no_SNPs=poisson() ; // newly produced cell mutants
        if (_drand48()>genotypes[cells[n].gen]->m[treatment]) { // make a new cell in the same lesion
          Cell c ; c.x=kn-wx/2 ; c.y=jn-wx/2 ; c.z=in-wx/2 ; c.lesion=cells[n].lesion ;
#ifndef PUSHING
          ll->p[in*wx+jn]->set(kn) ;
#else
          ll->p[in][jn][kn]=cells.size() ;
#endif
          if (no_SNPs>0) { 
            c.gen=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ; // mutate 
          } else { 
            c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
          }
          cells.push_back(c) ; volume++ ;
          ll->n++ ; 
#ifndef NO_MECHANICS
          double d=(c.x*c.x+c.y*c.y+c.z*c.z) ; if (d>SQR(ll->rad)) ll->rad=sqrt(d) ;
          if (ll->rad/ll->rad0>1.05) {
            ll->reduce_overlap() ;  
            ll->find_closest() ; 
            ll->rad0=ll->rad ;
            ll->n0=ll->n ; 
          }
#endif
        } else { // make a new lesion
          int x=kn-wx/2+ll->r.x, y=jn-wx/2+ll->r.y, z=in-wx/2+ll->r.z ;
          if (no_SNPs>0) { 
            genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
            lesions.push_back(new Lesion(cells.size(),genotypes.size()-1,x,y,z)) ;
          } else {
            genotypes[cells[n].gen]->number++ ; 
            lesions.push_back(new Lesion(cells.size(),cells[n].gen,x,y,z)) ;
          }        
#ifndef NO_MECHANICS
          lesions[lesions.size()-1]->find_closest() ; 
#endif
        }
// BOTH_MUTATE          
        no_SNPs=poisson() ; // old cell mutates
        if (no_SNPs>0) { 
          genotypes[cells[n].gen]->number-- ; 
          int pn=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
          cells[n].gen=genotypes.size()-1 ;
          if (genotypes[cells[n].gen]->number<=0) { 
            delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
          }
        }
      }
    }
#ifdef CORE_IS_DEAD
    if (ll->no_free_sites(k,j,i)==0) { // remove cell from the core but leave p[i,j,k] set
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      cells[n]=cells[cells.size()-1] ; cells.pop_back() ; 
    }
#endif

// now we implement death
//    if (tsc>1./max_growth_rate) tsc=1./max_growth_rate ; // this is an alternative way but it does not change anything so no need to use it
#ifdef DEATH_ON_SURFACE    
    if (genotypes[cells[n].gen]->death[treatment]>0 && _drand48()<tsc*genotypes[cells[n].gen]->death[treatment]*ll->no_free_sites(k,j,i)/float(_nonn))  { // death on the surface
#else
    if (_drand48()<tsc*genotypes[cells[n].gen]->death[treatment]) { // death in volume
#endif
#ifndef PUSHING
      ll->p[i*wx+j]->unset(k) ;
#else
      ll->p[i][j][k]=-1 ;
#endif
      ll->n-- ; 
//      if (ll->n<0) err("ll->n<0") ;
#ifndef NO_MECHANICS
      if (ll->n>1000 && 1.*ll->n/ll->n0<0.9) { // recalculate radius
        ll->rad=0 ; 
        for (i=0;i<wx;i++) for (j=0;j<wx;j++) for (k=0;k<wx;k++) {
#ifndef PUSHING
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i*wx+j]->is_set(k) && d>SQR(ll->rad)) ll->rad=sqrt(d) ; 
#else
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i][j][k]!=-1 && d>SQR(ll->rad)) ll->rad=sqrt(d) ; 
#endif
        }
        ll->rad0=ll->rad ; ll->n0=ll->n ;
      }
#endif
      if (ll->n==0) {
        int nn=cells[n].lesion ;
//        if (ll!=lesions[nn]) err("ll!") ;
        ll=NULL ; 
        delete lesions[nn] ; 
        if (nn!=lesions.size()-1) { // move lesion to a different index, and change cells->lesion correspondingly
          for (i=0;i<cells.size();i++) if (cells[i].lesion==lesions.size()-1) cells[i].lesion=nn ;
          lesions[nn]=lesions[lesions.size()-1] ; 
        }        
        lesions.pop_back() ;        
#ifndef NO_MECHANICS
        for (i=0;i<lesions.size();i++) {
          lesions[i]->find_closest() ;  
        }
#endif 
      }
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      if (n!=cells.size()-1) { 
        cells[n]=cells[cells.size()-1] ;
#ifdef PUSHING
        Lesion *ll2=lesions[cells[n].lesion] ;
        int ii=cells[n].z+ll2->wx/2, jj=cells[n].y+ll2->wx/2, kk=cells[n].x+ll2->wx/2 ;
        ll2->p[ii][jj][kk]=n ;
#endif
      }
      cells.pop_back() ; volume-- ;
      //if (lesions.size()==0) err("N=",int(cells.size())) ;
    }

    if (need_wx_update && ll!=NULL) ll->update_wx() ;    
      
#ifdef CORE_IS_DEAD
    ntot=volume ;
#else
    ntot=cells.size() ;
#endif

    if (wait_time>0 && tt>tt_old+wait_time) { tt_old=tt ; save_data(); }
    if (save_size>1 && ntot>=save_size) { save_size*=2 ; save_data() ; }

    if (cells.size()==0) {cout<<"no cells"<<endl; return 1 ;} 
    if (max_time>0 && tt>max_time) {cout<<"time over"<<endl; return 3 ;}
    if (exit_size>0 && ntot>=exit_size) {cout<<"max size reached"<<endl; return 4 ;}
    
   // cout<<cells.size()<<" "<<tt<<" "<<next_time<<endl;
    
    if(cells.size()>200 && tt>next_time)
    {
    
    ofstream myfile;
    string fn="output/tumour_sim_1e6_"+std::to_string(counter)+".txt";
    myfile.open (fn, ios::out); 
    
    myfile<<cells.size()<<endl;
    myfile<<tt<<endl<<endl;
    
    for (int i=0;i<cells.size();i++) {
      Lesion *ll=lesions[cells[i].lesion] ;
      Genotype *g=genotypes[cells[i].gen] ;
      myfile<<int(genotypes[cells[i].gen]->index)<<endl;
      myfile <<int(cells[i].x+ll->r.x)<<endl;
      myfile <<int(cells[i].y+ll->r.x)<<endl;
      myfile <<int(cells[i].z+ll->r.x)<<endl;
    }
    myfile.close();
    next_time+=0.7; counter++;
    }
    
  }

}
#endif // NORMAL


//-----------------------------------------------------------------------------

#if defined(FASTER_KMC) || defined(GILLESPIE)

int main_proc(int exit_size, int save_size, double max_time, double wait_time)
{
  int i,j,k,n,l,in,jn,kn,ntot;  
  int cc=0, timeout=0 ;
  double tt_old=tt ;

  for(;;) {      // main loop 
#ifdef PAUSE_WHEN_MEMORY_LOW
    timeout++ ; if (timeout>1000000) {
      timeout=0 ; 
      while (freemem()<PAUSE_WHEN_MEMORY_LOW) { sleep(1) ; } 
    }    
#endif

#ifdef GILLESPIE
  // Gillespie
    double tot_rate=0 ;
    for (n=0;n<cells.size();n++) {
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
      tot_rate+=genotypes[cells[n].gen]->growth[treatment] * free_sites(n)/float(_nonn) ;
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
      if (free_sites(n)>0) tot_rate+=genotypes[cells[n].gen]->growth[treatment] ;
#elif defined(PUSHING)
      tot_rate+=genotypes[cells[n].gen]->growth[treatment] ;
#else
      #error inconsistent growth conditions 
#endif
#ifdef DEATH_ON_SURFACE    
      tot_rate+=genotypes[cells[n].gen]->death[treatment] * free_sites(n)/float(_nonn) ; // death on the surface
#else
      tot_rate+=genotypes[cells[n].gen]->death[treatment] ;  // death in volume
#endif
    }        
    tt+=-log(1-_drand48())*timescale/tot_rate ; 
    double q=_drand48()*tot_rate,r=0 ;
    int mode=0 ;
    for (n=0;n<cells.size();n++) {
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
      r+=genotypes[cells[n].gen]->growth[treatment] * free_sites(n)/float(_nonn) ;
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
      if (free_sites(n)>0) r+=genotypes[cells[n].gen]->growth[treatment] ;
#elif defined(PUSHING)
      r+=genotypes[cells[n].gen]->growth[treatment] ;
#else
      #error inconsistent growth conditions 
#endif
      if (r>q) break ;
#ifdef DEATH_ON_SURFACE    
      r+=genotypes[cells[n].gen]->death[treatment] * free_sites(n)/float(_nonn)  ; // death on the surface
#else
      r+=genotypes[cells[n].gen]->death[treatment] ;  // death in volume
#endif
      if (r>q) { mode=1 ; break ; }
    }
    if (n==cells.size()) err("n==cells.size() at t=",tt) ;
#endif

#ifdef FASTER_KMC
    double max_death_rate=1 ;
    double tot_rate=cells.size()*(max_growth_rate+max_death_rate) ;
    tt+=-log(1-_drand48())*timescale/tot_rate ; 
    n=_drand48()*cells.size() ;
    double q=_drand48()*(max_growth_rate+max_death_rate), br,dr  ;
    int mode=0 ; 
#if !defined(CONST_BIRTH_RATE) && !defined(PUSHING)
    br=genotypes[cells[n].gen]->growth[treatment] * free_sites(n)/float(_nonn) ;
#elif (!defined(PUSHING)) && (defined(CONST_BIRTH_RATE))
    if (free_sites(n)>0) br=genotypes[cells[n].gen]->growth[treatment] ; else br=0 ;
#elif defined(PUSHING)
    br=genotypes[cells[n].gen]->growth[treatment] ;
#else
      #error inconsistent growth conditions 
#endif
#ifdef DEATH_ON_SURFACE    
    dr=genotypes[cells[n].gen]->death[treatment] * free_sites(n)/float(_nonn) ; // death on the surface
#else
    dr=genotypes[cells[n].gen]->death[treatment] ;  // death in volume
#endif
    if (q<br) mode=0 ;
    else if (q<br+dr) mode=1 ;
    else mode=2 ;
#endif    

    Lesion *ll=lesions[cells[n].lesion] ;
    int wx=ll->wx ; 
    k=cells[n].x+wx/2 ; j=cells[n].y+wx/2 ; i=cells[n].z+wx/2 ; 
    int need_wx_update=0 ;
    if (k<2 || k>=ll->wx-3 || j<2 || j>=ll->wx-3 || i<2 || i>=ll->wx-3) need_wx_update=1 ; 

    if (mode==0) { // reproduction
#if !defined(PUSHING)
      int in=i, jn=j, kn=k ;
      ll->choose_nn(kn,jn,in) ;
#else
      int in,jn,kn ;
      ll->find_dir_min_drag(i,j,k, in,jn,kn) ;      
#endif

      int no_SNPs=poisson() ; // newly produced cell mutants
      if (_drand48()>genotypes[cells[n].gen]->m[treatment]) { // make a new cell in the same lesion
        Cell c ; c.x=kn-wx/2 ; c.y=jn-wx/2 ; c.z=in-wx/2 ; c.lesion=cells[n].lesion ;
#ifdef PUSHING
        ll->p[in][jn][kn]=cells.size() ;
#else
        ll->p[in*wx+jn]->set(kn) ;
#endif
        if (no_SNPs>0) { 
          c.gen=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ; // mutate 
        } else { 
          c.gen=cells[n].gen ; genotypes[cells[n].gen]->number++ ; 
        }
        cells.push_back(c) ; volume++ ;

        ll->n++ ; 
#ifndef NO_MECHANICS
        double d=(c.x*c.x+c.y*c.y+c.z*c.z) ; if (d>SQR(ll->rad)) ll->rad=sqrt(d) ;
        if (ll->rad/ll->rad0>1.05) {
          ll->reduce_overlap() ;  
          ll->find_closest() ; 
          ll->rad0=ll->rad ;
          ll->n0=ll->n ; 
        }
#endif
      } else { // make a new lesion
        int x=kn-wx/2+ll->r.x, y=jn-wx/2+ll->r.y, z=in-wx/2+ll->r.z ;
        if (no_SNPs>0) { 
          genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
          lesions.push_back(new Lesion(cells.size(),genotypes.size()-1,x,y,z)) ;
        } else {
          genotypes[cells[n].gen]->number++ ; 
          lesions.push_back(new Lesion(cells.size(),cells[n].gen,x,y,z)) ;
        }        
#ifndef NO_MECHANICS
        lesions[lesions.size()-1]->find_closest() ; 
#endif
      }
// BOTH_MUTATE          
      no_SNPs=poisson() ; // old cell mutates
      if (no_SNPs>0) { 
        genotypes[cells[n].gen]->number-- ; 
        int pn=genotypes.size() ; genotypes.push_back(new Genotype(genotypes[cells[n].gen],cells[n].gen,no_SNPs)) ;
        cells[n].gen=genotypes.size()-1 ;
        if (genotypes[cells[n].gen]->number<=0) { 
          delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
        }
      }
    }
  
#ifdef CORE_IS_DEAD
    if (ll->no_free_sites(k,j,i)==0) { // remove cell from the core but leave p[i,j,k] set
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      cells[n]=cells[cells.size()-1] ; cells.pop_back() ; 
    }
#endif

// now we implement death
    if (mode==1) {
#ifdef PUSHING
      ll->p[i][j][k]=-1 ;
#else
      ll->p[i*wx+j]->unset(k) ;
#endif
      ll->n-- ; 
//      if (ll->n<0) err("ll->n<0") ;
#ifndef NO_MECHANICS
      if (ll->n>1000 && 1.*ll->n/ll->n0<0.9) { // recalculate radius
        ll->rad=0 ; 
        for (i=0;i<wx;i++) for (j=0;j<wx;j++) for (k=0;k<wx;k++) {
#ifdef PUSHING
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i][j][k]!=-1 && d>SQR(ll->rad)) ll->rad=sqrt(d) ; 
#else
          double d=SQR(i-wx/2)+SQR(j-wx/2)+SQR(k-wx/2) ; if (ll->p[i*wx+j]->is_set(k) && d>SQR(ll->rad)) ll->rad=sqrt(d) ; 
#endif
        }
        ll->rad0=ll->rad ; ll->n0=ll->n ;
      }
#endif
      if (ll->n==0) {
        int nn=cells[n].lesion ;
//        if (ll!=lesions[nn]) err("ll!") ;
        ll=NULL ; 
        delete lesions[nn] ; 
        if (nn!=lesions.size()-1) { // move lesion to a different index, and change cells->lesion correspondingly
          for (i=0;i<cells.size();i++) if (cells[i].lesion==lesions.size()-1) cells[i].lesion=nn ;
          lesions[nn]=lesions[lesions.size()-1] ; 
        }        
        lesions.pop_back() ;        
#ifndef NO_MECHANICS
        for (i=0;i<lesions.size();i++) {
          lesions[i]->find_closest() ;  
        }
#endif 
      }
      genotypes[cells[n].gen]->number-- ; if (genotypes[cells[n].gen]->number<=0) { 
        delete genotypes[cells[n].gen] ; genotypes[cells[n].gen]=NULL ; 
      }
      if (n!=cells.size()-1) { 
        cells[n]=cells[cells.size()-1] ;
#ifdef PUSHING
        Lesion *ll2=lesions[cells[n].lesion] ;
        int ii=cells[n].z+ll2->wx/2, jj=cells[n].y+ll2->wx/2, kk=cells[n].x+ll2->wx/2 ;
        ll2->p[ii][jj][kk]=n ;
#endif
      }
      cells.pop_back() ; volume-- ;
      //if (lesions.size()==0) err("N=",int(cells.size())) ;
    }

    if (need_wx_update && ll!=NULL) ll->update_wx() ;    
      
#ifdef CORE_IS_DEAD
    ntot=volume ;
#else
    ntot=cells.size() ;
#endif

    if (wait_time>0 && tt>tt_old+wait_time) { tt_old=tt ; save_data(); }
    if (save_size>1 && ntot>=save_size) { save_size*=2 ; save_data() ; }

    if (cells.size()==0) return 1 ; 
    if (max_time>0 && tt>max_time) return 3 ;
    if (exit_size>0 && ntot>=exit_size) return 4 ;

  }

}

#endif // FASTER_KMC or GILLESPIE

// ------- END FILE: simulations_2025/simulation.cpp -------

// ------- BEGIN FILE: simulations_2025/functions.cpp -------
/*******************************************************************************
   TumourSimulator v.1.2.2 - a program that simulates a growing solid tumour.
   Based on the algorithm described in
   
   Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban, 
   Bert Vogelstein, and Martin A. Nowak. "Spatial Model Predicts That 
   Dispersal and Cell Turnover Limit Intratumour Heterogeneity" Nature 525, 
   no. 7568 (September 10, 2015): 261-64. doi:10.1038/nature14971.

   Contributing author:
   Dr Bartek Waclaw, University of Edinburgh, bwaclaw@staffmail.ed.ac.uk

   Copyright (2015) The University of Edinburgh.

    This file is part of TumourSimulator.

    TumourSimulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TumourSimulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file 
    License.txt or at <http://www.gnu.org/licenses/>.
*******************************************************************************/


#include <math.h>
// (include "params.h" elided; inlined above)
// (include "classes.h" elided; inlined above)

void save_snps(char *name,int *n, int total, int mode, int *most_abund) 
{
  const float cutoff=0.01 ;
  FILE *f=fopen(name,"w") ;
  if (f==NULL) err(name) ;
  int i, j, nsnps=0, nsnpsc=0;
  for (i=0;i<L;i++) { 
    if (n[i]>(1e-4)*total) nsnps++ ;
    if (1.*n[i]/total>cutoff) nsnpsc++ ;
  }
  float *abund=new float[nsnps], tempd ;
  int *num=new int[nsnps], temp ;
  nsnps=0 ;
  for (i=0;i<L;i++) if (n[i]>(1e-4)*total) { num[nsnps]=i ; abund[nsnps]=float(1.*n[i]/total)*(1+0.000001*i/L) ; nsnps++ ; }
  quicksort2(abund,num,0,nsnps-1) ;

  if (mode) {
    for (i=0;i<nsnps;i++) fprintf(f,"%d %d %f\n",i,num[i],abund[i]) ;
  } else {
    for (i=0;i<nsnps;i++) if (abund[i]>cutoff || i<100) fprintf(f,"%d %d %f\n",i,num[i],abund[i]) ;
  }
  if (most_abund!=NULL) { // store first 100 most abundant PMs
    for (i=0;i<MIN(100,nsnps);i++) most_abund[i]=num[i] ; 
  }
  delete [] abund ; delete [] num ;
  fclose(f) ;
}

void save_snp_corr(char *name, Hist *snps)
{
  FILE *f=fopen(name,"w") ;
  int last ;
  for (last=_bins-1;snps[last].x==0 && last>0;last--) ;
  for (int i=0;i<=last;i++) if (snps[i].n>0) fprintf(f,"%d %f %f\n",i*_resol,
      float(1.*snps[i].x/snps[i].n),float(sqrt(( (1.*snps[i].x2/snps[i].n) - (1.*snps[i].x/snps[i].n)*(1.*snps[i].x/snps[i].n))/(snps[i].n-1)))) ; 
    else fprintf(f,"%d 0 0\n",i*_resol) ;
  fclose(f) ;    
}



int how_many_SNPs_identical(Genotype *a, Genotype *b)
{
  int i,j,n=0;
  int al=a->sequence.size(), bl=b->sequence.size() ;
  for (i=0;i<al;i++)
    for (j=0;j<bl;j++) 
      if (a->sequence[i]==b->sequence[j]) n++ ;
  return n ;
}

int how_many_SNPs_identical(Genotype *a, Genotype *b, float cutoff, int *snp_no)
{
  int i,j,n=0, ntot=cells.size();
  int al=a->sequence.size(), bl=b->sequence.size() ;
  for (i=0;i<al;i++)
    for (j=0;j<bl;j++) 
      if (a->sequence[i]==b->sequence[j] && snp_no[(a->sequence[i])&L_PM]>cutoff*ntot) n++ ;
  return n ;
}


void select_two_random_cells(int &i, int &j, int &n)
{
  int ntot=cells.size() ;
  double dist ;
  do {
    i=int(_drand48()*ntot) ; j=int(_drand48()*ntot) ;
    vecd ri=lesions[cells[i].lesion]->r, rj=lesions[cells[j].lesion]->r, rij=ri-rj ;
    dist=SQR(rij.x+cells[i].x-cells[j].x)+SQR(rij.y+cells[i].y-cells[j].y)+SQR(rij.z+cells[i].z-cells[j].z) ;
  } while (_drand48()<0.5 && dist>_drand48()*100) ; // choose preferentially cells which are close to one another
  dist=sqrt(dist) ;
  
  n=int(dist/_resol) ;
  if (n<0 || n>=_bins) err("n",n);
}

float average_distance_ij()
{
  int ntot=cells.size() ;
  double avdist=0 ;
  for (int n=0;n<ntot;n++) {
    int i=int(_drand48()*ntot), j=int(_drand48()*ntot) ; 
    vecd ri=lesions[cells[i].lesion]->r, rj=lesions[cells[j].lesion]->r, rij=ri-rj ;
    avdist+=sqrt(SQR(rij.x+cells[i].x-cells[j].x)+SQR(rij.y+cells[i].y-cells[j].y)+SQR(rij.z+cells[i].z-cells[j].z)) ;
  } 
  return (avdist/ntot) ;
}

void snps_corr(Hist *snps)
{
  int i,j,k,n;
  int ntot=cells.size() ;
    
  for (i=0;i<_bins;i++) snps[i].x=snps[i].x2=snps[i].n=0 ;
  for (k=0;k<ntot;k++) {
    select_two_random_cells(i,j,n) ;
    snps[n]+=how_many_SNPs_identical(genotypes[cells[i].gen],genotypes[cells[j].gen]) ;
  }  
}

void snps_corr_cutoff(Hist *snps,float cutoff,int *snp_no) 
{
  int i,j,k,n;
  int ntot=cells.size() ;
    
  for (i=0;i<_bins;i++) snps[i].x=snps[i].x2=snps[i].n=0 ;
  for (k=0;k<ntot;k++) {
    select_two_random_cells(i,j,n) ;
    snps[n]+=how_many_SNPs_identical(genotypes[cells[i].gen],genotypes[cells[j].gen],cutoff,snp_no) ;
  }  
  
}


void snps_corr_cond_driver(Hist *snps)
{
  int ntot=cells.size() ;
  int i,j,k,n,d,ai,aj;
    
  for (k=0;k<ntot;k++) {
    select_two_random_cells(i,j,n) ;
    Genotype *gi=genotypes[cells[i].gen], *gj=genotypes[cells[j].gen]  ;    
    d=0 ;
    for (ai=0;ai<gi->sequence.size();ai++) if (gi->sequence[ai]<0) {
      for (aj=0;aj<gj->sequence.size();aj++) if (gi->sequence[ai]==gj->sequence[aj]) { d++ ; break ; }
    }
    if (d>0) {
      snps[n]+=how_many_SNPs_identical(genotypes[cells[i].gen],genotypes[cells[j].gen]) ;
    }
  }  
}


void find_p_driver(Hist *pr1, Hist *pr2, Hist *pr3) 
{
//Probability of finding the same driver in two cells separated by some distance x.
// 1) prob. that two cells at distance x have at least one common driver
// 2) prob. that two cells at distance x have the same last driver
// 3) prob. that two cells at distance x have the same first driver
  int i,j,k,n,ai,aj,d,di,dj;
  int ntot=cells.size() ;
    
  for (k=0;k<ntot;k++) {
    select_two_random_cells(i,j,n) ;
    Genotype *gi=genotypes[cells[i].gen], *gj=genotypes[cells[j].gen]  ;    
    
    d=0 ; 
    for (ai=0;ai<gi->sequence.size();ai++) if (gi->sequence[ai]<0) {
      di++ ;
      for (aj=0;aj<gj->sequence.size();aj++) if (gi->sequence[ai]==gj->sequence[aj]) { d++ ; break ; }
    }
    pr1[n]+=(d>0?1:0) ;

    d=0 ;
    for (ai=gi->sequence.size()-1;ai>=0;ai--) if (gi->sequence[ai]<0) {
      for (aj=0;aj<gj->sequence.size();aj++) if (gi->sequence[ai]==gj->sequence[aj]) { d++ ; break ; }
      break ;
    }
    pr2[n]+=(d>0?1:0) ;

    d=0 ;
    for (ai=0;ai<gi->sequence.size();ai++) if (gi->sequence[ai]<0) {
      for (aj=0;aj<gj->sequence.size();aj++) if (gi->sequence[ai]==gj->sequence[aj]) { d++ ; break ; }
      break ;
    }
    pr3[n]+=(d>0?1:0) ;
        
  }  
}

// ------- END FILE: simulations_2025/functions.cpp -------

// ------- BEGIN FILE: simulations_2025/main.cpp -------
/*******************************************************************************
   TumourSimulator v.1.2.2 - a program that simulates a growing solid tumour.
   Based on the algorithm described in
   
   Bartlomiej Waclaw, Ivana Bozic, Meredith E. Pittman, Ralph H. Hruban, 
   Bert Vogelstein, and Martin A. Nowak. ?A Spatial Model Predicts That 
   Dispersal and Cell Turnover Limit Intratumour Heterogeneity.? Nature 525, 
   no. 7568 (September 10, 2015): 261?64. doi:10.1038/nature14971.

   Contributing author:
   Dr Bartek Waclaw, University of Edinburgh, bwaclaw@staffmail.ed.ac.uk

   Copyright (2015) The University of Edinburgh.

    This file is part of TumourSimulator.

    TumourSimulator is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TumourSimulator is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
    See the GNU General Public License for more details.

    A copy of the GNU General Public License can be found in the file 
    License.txt or at <http://www.gnu.org/licenses/>.
*******************************************************************************/


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
// (include "classes.h" elided; inlined above)

#define __MAIN
// (include "params.h" elided; inlined above)

#if (defined(GILLESPIE) && defined(FASTER_KMC)) || (defined(GILLESPIE) && defined(NORMAL)) || (defined(NORMAL) && defined(FASTER_KMC))
  #error too many methods defined!
#endif

#if !defined(GILLESPIE) && !defined(FASTER_KMC) && !defined(NORMAL) 
  #error no method defined!
#endif

#if (defined(VON_NEUMANN_NEIGHBOURHOOD) && defined(MOORE_NEIGHBOURHOOD))
  #error both VON_NEUMANN_NEIGHBOURHOOD and MOORE_NEIGHBOURHOOD defined!
#endif

#if (!defined(VON_NEUMANN_NEIGHBOURHOOD) && !defined(MOORE_NEIGHBOURHOOD))
  #error neither VON_NEUMANN_NEIGHBOURHOOD nor MOORE_NEIGHBOURHOOD defined!
#endif

extern char *NUM ; 
extern int RAND, sample, treatment, max_size ;
extern double tt ;
extern float time_to_treat ;

int sample=0 ;

void save_positions(char *name, float dz) 
{
  ofstream myfile;
  myfile.open (name, ios::out); 
  for (int i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    Genotype *g=genotypes[cells[i].gen] ;
    myfile<<int(genotypes[cells[i].gen]->index)<<endl;
    myfile <<int(cells[i].x+ll->r.x)<<endl;
    myfile <<int(cells[i].y+ll->r.x)<<endl;
    myfile <<int(cells[i].z+ll->r.x)<<endl;
    }
  myfile.close();

  /*  
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    Genotype *g=genotypes[cells[i].gen] ;
    //if (abs(int(cells[i].z+ll->r.z))<dz || cells.size()<2e5) //
      //for max 182 colors use %182
    //fprintf(data,"%d %d %d %d \n",int(genotypes[cells[i].gen]->index%182), int(cells[i].x+ll->r.x), int(cells[i].y+ll->r.y),int(cells[i].z+ll->r.z)) ;//genotypes[cells[i].gen]->index
      fprintf(data,"%d %d %d %d \n",int(genotypes[cells[i].gen]->index), int(cells[i].x+ll->r.x), int(cells[i].y+ll->r.y), int(cells[i].z+ll->r.z)) ;
  }        
  fclose(data) ; */
}

float save_2d_image(char *name, vecd li)
{
  int i,j,k;
  int minx=1<<20,maxx=-minx,miny=minx,maxy=-minx,minz=minx,maxz=-minx ;  
  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    if (cells[i].x+ll->r.x<minx) minx=int(cells[i].x+ll->r.x) ;
    if (cells[i].x+ll->r.x>maxx) maxx=int(cells[i].x+ll->r.x) ;
    if (cells[i].y+ll->r.y<miny) miny=int(cells[i].y+ll->r.y) ;
    if (cells[i].y+ll->r.y>maxy) maxy=int(cells[i].y+ll->r.y) ;
    if (cells[i].z+ll->r.z<minz) minz=int(cells[i].z+ll->r.z) ;
    if (cells[i].z+ll->r.z>maxz) maxz=int(cells[i].z+ll->r.z) ;
  }
  maxx++ ; maxy++ ; 
  float density=1.*cells.size()/(float(maxx-minx)*float(maxy-miny)*float(maxz-minz)) ;
  if (cells.size()<1e3) return density ;
  float diam=pow(float(maxx-minx)*float(maxy-miny)*float(maxz-minz),1./3) ;
  printf("density=%f\n",density) ;

  int nnn=(maxx-minx)*(maxy-miny) ;
  if (float(maxx-minx)*float(maxy-miny)>2e9) err("(maxx-minx)*(maxy-miny) too large",float(maxx-minx)*float(maxy-miny)) ;
  int *types=new int[nnn] ; // 2d array of cells' types (foremost ones only)
  short int *zbuf=new short int[nnn] ; // 2d array of cells' z positions (foremost ones only)
  BYTE *br=new BYTE[nnn] ; // 2d array of cells' brightness (foremost ones only)
  Sites **bit=new Sites*[nnn] ; // 3d array of bits: cell empty/occupied
  for (i=0;i<nnn;i++) { zbuf[i]=minz ; types[i]=-1 ; br[i]=255 ; bit[i]=new Sites(maxz-minz) ; }
  printf("%d %d\t %d %d\n",minx,maxx,miny,maxy) ;
  printf("%d x %d = %d\n",maxx-minx,maxy-miny,nnn) ;

  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    int z=int(cells[i].z+ll->r.z) ;
    int adr=int(cells[i].y+ll->r.y-miny)*(maxx-minx)+int(cells[i].x+ll->r.x-minx) ;
    if (adr<0 || adr>=nnn) err("adr",adr) ;
    bit[adr]->set(z-minz) ;
    if (z>zbuf[adr]) { zbuf[adr]=z ; types[adr]=genotypes[cells[i].gen]->index;}
  }
  
  normalize(li) ;
  float dmul=0.93 ; 
  float range=-(0.916291/(density*log(dmul))) ; 
  if (range>diam) { range=diam ; dmul=pow(0.4,1/(density*range)) ; }
  for (i=0;i<maxy-miny;i++) 
    for (j=0;j<maxx-minx;j++) {
      float d=1 ; 
      int adr=i*(maxx-minx)+j ;
      k=zbuf[adr]-minz ;
      for (float o=1;o<range;o++) {
        int il=int(i-li.y*o) ; int jl=int(j-li.x*o) ; int kl=int(k-li.z*o) ; 
        if (il>0 && jl>0 && il<maxy-miny && jl<maxx-minx && kl>0 && kl<maxz-minz) {          
          if (bit[il*(maxx-minx)+jl]->is_set(kl)) d*=dmul ;
        } else break ;
        if (d<0.4) { d=0.4 ; break ; }
      }
      br[adr]*=d ;
    }

  FILE *data=fopen(name,"w") ;    
  for (i=0;i<maxy-miny;i++) {
    for (j=0;j<maxx-minx;j++) fprintf(data,"%d %d ",types[i*(maxx-minx)+j],br[i*(maxx-minx)+j]) ;
    fprintf(data,"\n") ;
  }
  fclose(data) ; 
  delete [] types ; delete [] zbuf ; delete [] br ; delete [] bit ;
  return density ;
}
  

void save_genotypes(char *name)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *g=genotypes[i] ;
    if (g!=NULL && g->number>0) {
      fprintf(data,"%d  %d  %d %d  %d\t",i, g->prev_gen,g->no_resistant,g->no_drivers, g->number) ;
      for (int j=0;j<g->sequence.size();j++) fprintf(data," %u",g->sequence[j]) ; 
      fprintf(data,"\n") ;
    } 
  }

  fclose(data) ;  
}

void save_most_abund_gens(char *name, int *most_abund)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *gg=genotypes[i] ;
    if (gg!=NULL && gg->number>0) {
      int r=0,g=0,b=0 ;
      for (int j=0;j<gg->sequence.size();j++) {
        if ((gg->sequence[j]&L_PM)==most_abund[0]) r=1 ;
        if ((gg->sequence[j]&L_PM)==most_abund[1]) g=1 ; 
        if ((gg->sequence[j]&L_PM)==most_abund[2]) b=1 ;
      }
      if (r || g || b) fprintf(data,"%d %d %d\t%d\n",r,g,b,gg->index) ;
    }
  }
  fclose(data) ;  

}

int main(int argc, char *argv[])
{
#if defined(GILLESPIE)   
  cout <<"method: GILLESPIE\n" ;
#endif
#if defined(FASTER_KMC)   
  cout <<"method: FASTER_KMC\n" ;
#endif
#if defined(NORMAL)   
  cout <<"method: NORMAL\n" ;
#endif

    int nsam ; int filename;
  if (argc!=4) { err(" Error:: arguments needed: name, no_samples, RAND. Program terminated. \n"); }
  else { 
    NUM=argv[1] ;
    nsam=atoi(argv[2]) ;
    RAND=atoi(argv[3]) ;
}
  cout <<NUM<<" "<<" "<<nsam<<" "<<RAND<<endl ;
  _srand48(RAND) ;
  init();
  char name[256] ;
  //sprintf(name,"%s/each_run_%d.dat",NUM,max_size) ;
  //FILE *er=fopen(name,"w") ; fclose(er) ;
  for (sample=0;sample<nsam;sample++) { 
    reset() ;
#ifdef MAKE_TREATMENT_N
    int s=0 ; while (main_proc(max_size,-1,-1, 10)==1) { s++ ; reset() ; } ; // initial growth until max size is reached, saved every 10 days
    if (s>0) printf("resetted %d times\n",s) ;
    save_data() ; 
    treatment=1 ;       
    double max_time=2*tt ;
    main_proc(1.25*max_size,-1,max_time, 10) ; // treatment
#elif defined MAKE_TREATMENT_T
    int s=0 ; while (main_proc(-1,-1,time_to_treat, 10)==1) { s++ ; reset() ; } ; // initial growth until max time is reached, saved every 10 days
    if (s>0) printf("resetted %d times\n",s) ;
    save_data() ; 
    treatment=1 ;       
    double max_time=2*tt ;
    max_size=cells.size()*1.25 ;
    main_proc(max_size,-1,max_time, 10) ; // treatment

#else    
    int s=0 ; while (main_proc(max_size,2,-1, -1)==1) { s++ ; reset() ; } // initial growth until max size is reached
    if (s>0) printf("resetted %d times\n",s) ;
    fflush(stdout) ;
    if (nsam==1) {  // do this only when making images of tumours
      printf("saving images...\n") ;
      int j=0 ;
     sprintf(name,"%s/tumour_%d_%d_%d.dat",NUM,max_size,sample,RAND) ; save_positions(name,1.0) ;
      } 
   #endif
  } 
  end() ;
	return 0 ;
}


