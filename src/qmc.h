// Frances Kuo
// last updated 28 August 2009

#ifndef _QMC
#define _QMC

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

using namespace std;


// -----------------------------------------------------------------------------------------
// Pseudo-random number generator
// -----------------------------------------------------------------------------------------
static double ran2(long *idum) 
{
  // Long period (> 2e18) random number generator of L'Ecuyer 
  // with Bays-Durham shuffle and added safeguards. 
  // Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). 
  // Call with idum a negative integer to initialize; 
  // thereafter, do not alter idum between successive deviates in a sequence. 
  // RNMX should approximate the largest floating value that is less than 1. 
  //
  // change into a double version

  static const long IM1 = 2147483563;    
  static const long IM2 = 2147483399;    
  static const double AM = 1.0L / IM1;    
  static const long IMM1 = (IM1-1);      
  static const long IA1 = 40014;          
  static const long IA2 = 40692;         
  static const long IQ1 = 53668;        
  static const long IQ2 = 52774;        
  static const long IR1 = 12211;         
  static const long IR2 = 3791;         
  static const long NTAB = 32;          
  static const long NDIV =(1+IMM1/NTAB); 
  //static const double EPS = 4.4e-16; //1.2e-7;    
  //static const double RNMX = (1.0-EPS); 

  int j; 
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  if (*idum <= 0) {              // Initialize. 
    if (-(*idum) < 1) *idum=1;   // Be sure to prevent idum = 0. 
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {    // Load the shuffle table (after 8 warm-ups).
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1; 
      if (*idum < 0) *idum += IM1; 
      if (j < NTAB) iv[j] = *idum;
    } 
    iy=iv[0]; 
  } 
  k=(*idum)/IQ1;                 // Start here when not initializing.
  *idum=IA1*(*idum-k*IQ1)-k*IR1; // Compute idum=(IA1*idum) % IM1 without overflows by Schrage's 
                                 // method. 
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2; 
  idum2=IA2*(idum2-k*IQ2)-k*IR2; // Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;                     // Will be in the range 0..NTAB-1.
  iy=iv[j]-idum2;        // Here idum is shuffled, idum and idum2 are combined to generate output. 
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;

  double temp=AM*iy;
  if (temp >= 1) cout << "WARNING: random 1" << endl;
  //if ( temp > RNMX) return RNMX; // Because users don't expect endpoint values.
  //else return temp; 
  return temp;
} 

// -----------------------------------------------------------------------------------------
// inverse normal
//   see Peter J. Acklam  http://home.online.no/~pjacklam/notes/invnorm/
//   this is supposed to give full machine precision
// -----------------------------------------------------------------------------------------
static double inverse_normal(double p)
{
  static const double a1 = -3.969683028665376e+01;
  static const double a2 =  2.209460984245205e+02;
  static const double a3 = -2.759285104469687e+02;
  static const double a4 =  1.383577518672690e+02;
  static const double a5 = -3.066479806614716e+01;
  static const double a6 =  2.506628277459239e+00;
 
  static const double b1 = -5.447609879822406e+01;
  static const double b2 =  1.615858368580409e+02;
  static const double b3 = -1.556989798598866e+02;
  static const double b4 =  6.680131188771972e+01;
  static const double b5 = -1.328068155288572e+01;
  
  static const double c1 = -7.784894002430293e-03;
  static const double c2 = -3.223964580411365e-01;
  static const double c3 = -2.400758277161838e+00;
  static const double c4 = -2.549732539343734e+00;
  static const double c5 =  4.374664141464968e+00;
  static const double c6 =  2.938163982698783e+00;

  static const double d1 =  7.784695709041462e-03;
  static const double d2 =  3.224671290700398e-01;
  static const double d3 =  2.445134137142996e+00;
  static const double d4 =  3.754408661907416e+00;
  
  static const double p_low = 0.02425;
  static const double p_high = 1 - p_low;

  double q, r, x, e, u;

  if (p <= 0) {
    cout << "WARNING: inverse normal <=0 " << endl;
    return -HUGE_VAL; /* - infinity */
  }
  if (p >= 1) {
    cout << "WARNING: inverse normal >=1 " << endl;
    return HUGE_VAL;  /*   infinity */
  }
  
  /* Rational approximation for lower region */
  if (0 < p && p < p_low) {
    q = sqrt(-2*log(p));
    x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  }

  /* Rational approximation for central region */
  else if (p_low <= p && p <= p_high) {
    q = p - 0.5;
    r = q*q;
    x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
  }

  /* Rational approximation for upper region */
  else {
    q = sqrt(-2*log(1-p));
    x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  }

  /* One iteration of Halley's rational method */
  e = 0.5 * erfc(-x/sqrt(2.0)) - p;
  u = e * sqrt(2*M_PI) * exp(x*x/2);
  x = x - u/(1 + x*u/2);

  return x;
}

static double *inverse_normal(double *X, unsigned size)
{
  // index from 1 to size
  double *Z = new double [size+1];
  for (unsigned j=1;j<=size;j++) Z[j] = inverse_normal(X[j]);
  return Z;
}


// -----------------------------------------------------------------------------------------------------------------
// Other
// -----------------------------------------------------------------------------------------------------------------
static unsigned pow2(unsigned m)
{
  unsigned p = 1;
  for (unsigned i=1;i<=m;i++) p *= 2;
  return p;
}

static unsigned log2(unsigned x)
{
  unsigned p = 0;
  while (x != 1) {
    x /= 2;
    p++;
  }
  return p;
}

// ---------------------------------------------------------------------------------------------------------------------
// Points
//           index of M goes from 0 to M-1
//           index of N goes from 0 to N-1
//           index of D goes from 1 to D
// ---------------------------------------------------------------------------------------------------------------------
class Points {
protected:
  unsigned N;  // index of N goes from 0 to N-1
  unsigned D;  // index of D goes from 1 to D
  unsigned W;
  double **PNTS;

public:
  Points(unsigned n, unsigned d, unsigned w)
  {
    N = n;
    D = d;
    W = w;
    PNTS = NULL;
  }
  virtual ~Points()
  {
    if (PNTS != NULL) Points::free_points(PNTS,N,D);
  }

  virtual void init_point()
  {
  }

  virtual double *next_point()
  { 
  }

  virtual void generate_points()
  {
  }

  //double **get_points()
  //{
  //  if (PNTS == NULL) generate_points();
  //  return Points::copy_points(PNTS,N,D);
  //}
  double **get_points(double *S, bool digital = false)
  {
    if (PNTS == NULL) generate_points();
    return Points::shift_points(PNTS,S,N,D,digital);
  }
  double ***get_points(double **S, unsigned m, bool digital = false)
  {
    if (PNTS == NULL) generate_points();
    return Points::shift_points(PNTS,S,m,N,D,digital);
  }

  // -- STATIC MEMBER FUNCTIONS -----------------------------------------------------------
  static double *get_points(unsigned d, long &seed)   // the seed is changed afterward
  {
    //cout << "get_points: " << d << ", " << seed << endl;
    double *X = new double [d+1];
    for (unsigned j=1;j<=d;j++) X[j] = ran2(&seed);   
    return X;
  }
  static double **get_points(unsigned n, unsigned d, long &seed) 
  {
    double **X = new double * [n];
    for (unsigned i=0;i<=n-1;i++) {
      X[i] = get_points(d,seed);
    }
    return X;
  }
  static double ***get_points(unsigned m, unsigned n, unsigned d, long &seed) 
  {
    double ***X = new double ** [m];
    for (unsigned q=0;q<=m-1;q++) X[q] = get_points(n,d,seed);
    return X;
  }

  static void free_points(double *&X, unsigned d=0)
  {
    delete [] X;
    X = NULL;
  }
  static void free_points(double **&X, unsigned n, unsigned d=0)
  {
    for (unsigned i=0;i<=n-1;i++) delete [] X[i];
    delete [] X;
    X = NULL;
  }
  static void free_points(double ***&X, unsigned m, unsigned n, unsigned d=0)
  {
    for (unsigned p=0;p<=m-1;p++) {
      for (unsigned i=0;i<=n-1;i++) delete [] X[p][i];
      delete [] X[p];
    }
    delete [] X;
    X = NULL;
  }
  static void display_points(double *X, unsigned d)
  {
    for (unsigned j=1;j<=d;j++) cout << X[j] << " ";
    cout << endl;
  }
  static void display_points(double **X, unsigned n, unsigned d)
  {
    for (unsigned i=0;i<=n-1;i++) display_points(X[i],d);
  }
  static void display_points(double ***X, unsigned m, unsigned n, unsigned d)
  {
    for (unsigned p=0;p<=m-1;p++) {
      cout << "*** Set " << p+1 << " ***" << endl;
      display_points(X[p],n,d);
      cout << endl;
    }
  }
  static double *copy_points(double *X, unsigned d)
  {
    double *Y = new double [d+1];
    for (unsigned j=1;j<=d;j++) Y[j] = X[j];
    return Y;
  }
  static double **copy_points(double **X, unsigned n, unsigned d) 
  {
    double **Y = new double * [n];
    for (unsigned i=0;i<=n-1;i++) Y[i] = copy_points(X[i],d);
    return Y;
  }
  static double ***copy_points(double ***X, unsigned m, unsigned n, unsigned d) 
  {
    double ***Y = new double ** [m];
    for (unsigned q=0;q<=m-1;q++) Y[q] = copy_points(X[q],n,d);
    return Y;
  }
  static double *shift_points(double *X, double *S, unsigned d, bool digital=false)
  {
    double *Y = new double [d+1];
    if (digital) {
      for (unsigned j=1;j<=d;j++) {
	unsigned x = (unsigned)(X[j]*pow(2.0,31));
	unsigned s = (unsigned)(S[j]*pow(2.0,31));
	Y[j] = (double)(x ^ s)/pow(2.0,31);
	if (Y[j] == 0) {
	  //cout << "WARNING: digitally-shifting to zero! " 
	  //     << " ... relocate to 1e-32" 
	  //     << endl;
	  Y[j] = 1e-32L;  // inverse normal(1e-32) = -12 
	}
      }
    }
    else {
      for (unsigned j=1;j<=d;j++) {
	double tmp = X[j] + S[j];
	Y[j] = tmp - floor(tmp);
	if (Y[j] == 0) {
	  //cout << "WARNING: shifting to zero! " << X[j] << " + " << S[j] << " = " << tmp 
	  //     << " ... relocate to 1e-32" 
	  //     << endl;
	  Y[j] = 1e-32L;  // inverse normal(1e-32) = -12 
	}
      }
    }   
    return Y;
  }
  static double **shift_points(double **X, double *S, unsigned n, unsigned d, bool digital=false)
  {
    double **Y = new double * [n];
    for (unsigned i=0;i<=n-1;i++) Y[i] = shift_points(X[i],S,d,digital);
    return Y;
  }
  static double ***shift_points(double **X, double **S, unsigned m, unsigned n, unsigned d, bool digital=false)
  {
    double ***Y = new double ** [m];
    for (unsigned q=0;q<=m-1;q++) Y[q] = shift_points(X,S[q],n,d,digital);
    return Y;
  }

};

// ------------------------------------------------------------------------------------------
// Monte Carlo points
// ------------------------------------------------------------------------------------------
class MonteCarlo : public Points {
protected:
  long seed;

public:
  MonteCarlo(unsigned n, unsigned d, long s) : Points(n,d,12345)
  {
    seed = s;
  }
  virtual double *next_point()
  { 
    return Points::get_points(D,seed);
  }
  virtual void generate_points()
  {
    PNTS = Points::get_points(N,D,seed);
  }
};

// ------------------------------------------------------------------------------------------
// Lattice points
// ------------------------------------------------------------------------------------------
class Lattice : public Points {
//protected:
    public:
  unsigned *Z;
  unsigned counter;
  unsigned *Perm;

public:
  Lattice(unsigned n, unsigned d, unsigned w, unsigned *P = NULL) : Points(n,d,w)
  {
    Z = NULL;
    Perm = P;
  }  

  Lattice(unsigned n, unsigned d, unsigned* z) : Points(n,d,10000)
  {
    unsigned *gen = new unsigned [d+1];
    cout << "Setting Z directly" << endl;
    for (int i=0; i<d; i++) gen[i+1] = z[i];

    Z = gen;
    Perm = NULL;
  }  

  virtual ~Lattice()
  {
    if (Z != NULL) {
      delete [] Z;
      Z = NULL;
    }
    // do not delete Perm since it is passed in and could be shared
  }
  
  static unsigned *get_generator(unsigned n, unsigned d, unsigned w, unsigned *perm = NULL)
  {
    // Open input file
    ostringstream f;
    if ((w >= 20000 && w < 30000) || (w >= 40000 && w < 50000)) // standard lattice rules
      f << "/Users/james/projects/qmc_dtrw/of_QMC_DTRW/lattice-" << w << "-" << n << ends;
    else if ((w >= 30000 && w < 40000) || (w >= 50000 && w < 60000)) 
      f << "/Users/james/projects/qmc_dtrw/of_QMC_DTRW/lattice-" << w << "-1024-1048576" << ends; // extensible lattice rules up to 2^20
    else if ((w >= 60000 && w < 70000) ) 
      f << "porous_rules/lattice-pp-" << w << "." << n << ends; // porous flow weights 
    ifstream infile(f.str().c_str(),ios::in);
    if (!infile) {
      cout << "WARNING: No matching file " << f.str().c_str() << endl;
      exit(1);
    }

    // Read in generating vector
    unsigned *gen = new unsigned [d+1];
    unsigned j;
    unsigned z;
    double e;
    while (infile >> j >> z) {
      gen[j] = z;
      if (j == d) break;
    }
    if (j != d) {
      cout << "WARNING: No matching entry in the file!" << endl;
      exit(1);
    }
    infile.close();

    // Permute if necessary
    if (perm != NULL) {
      unsigned *tmp = gen;
      gen = new unsigned [d+1];
      for (unsigned j=1;j<=d;j++) gen[j] = tmp[perm[j]];
      delete [] tmp;
    }

    return gen;
  }

  virtual void generate_points()
  {
    if (Z == NULL) Z = get_generator(N,D,W,Perm);

    PNTS = new double * [N];
    for (unsigned i=0;i<=N-1;i++) {
      PNTS[i] = new double [D+1];
      for (unsigned j=1;j<=D;j++) {
	double x = (double)i * (double)Z[j] / N;
	PNTS[i][j] = x - floor(x);
      }
    }
  }

  virtual void init_point()
  { 
    if (Z == NULL) Z = get_generator(N,D,W,Perm);
    counter = 0;
  }
  virtual double *next_point()
  { 
    double *X = new double [D+1];
    for (unsigned j=1;j<=D;j++) {
      double x = (double)counter * (double)Z[j] / N;
      X[j] = x - floor(x);
    }
    counter++;
    return X;
  }
};

// ------------------------------------------------------------------------------------------
// Extensible Lattice points
// ------------------------------------------------------------------------------------------
class ExtensibleLattice : public Lattice {
protected:
  unsigned graycode;
public:
  ExtensibleLattice(unsigned n, unsigned d, unsigned w, unsigned *perm = NULL) : Lattice(n,d,w,perm)
  {
    // extensiblity has nothing to do with the choice of w
    // i.e. generator determined by w
  }  
  virtual void generate_points()
  {
    if (Z == NULL) Z = get_generator(N,D,W,Perm);
    PNTS = new double * [N];
    PNTS[0] = new double [D+1];
    for (unsigned j=1;j<=D;j++) PNTS[0][j] = 0;  // 0th point is the origin
    unsigned graycode = 0;
    for (unsigned i=1;i<=N-1;i++) {
      PNTS[i] = new double [D+1];
      unsigned index = 1;
      unsigned mask = 1 << 31;
      while ((i & index) == 0) {
	index <<= 1;
	mask >>= 1;
      }
      if ((graycode & mask) == 0) graycode |= mask;
      else graycode &= ~(mask);
      double code = (double)graycode/pow(2.0,32);
      //cout << code << endl;
      for (unsigned j=1;j<=D;j++) {
	double x = code * Z[j];
	PNTS[i][j] = x - floor(x);
      }
    }
  }

  virtual void init_point()
  { 
    if (Z == NULL) Z = get_generator(N,D,W,Perm);
    counter = 0;
    graycode = 0;
  }
  virtual double *next_point()
  { 
    double *X = new double [D+1];
    if (counter == 0) 
      for (unsigned j=1;j<=D;j++) X[j] = 0;
    
    else {
      unsigned index = 1;
      unsigned mask = 1 << 31;
      while ((counter & index) == 0) {
	index <<= 1;
	mask >>= 1;
      }
      if ((graycode & mask) == 0) graycode |= mask;
      else graycode &= ~(mask);
      double code = (double)graycode/pow(2.0,32);
      //cout << code << endl;
      for (unsigned j=1;j<=D;j++) {
	double x = code * Z[j];
	X[j] = x - floor(x);
      }
    }

    counter++;
    return X;
  }

};

// ------------------------------------------------------------------------------------------
// Sobol Points
// ------------------------------------------------------------------------------------------
class Sobol : public Points {
protected:
  unsigned maxN;
  unsigned **V;
  unsigned *XV;
  unsigned counter;
  unsigned *Perm;
  
public:
  Sobol(unsigned n, unsigned d, unsigned w, unsigned *perm = NULL) : Points(n,d,w)
  {
    maxN = n;
    //N = maxN + skip(maxN); // N <= pow(2,32)-1  -- This can be enlarged if we use "unsigned long"
    //generate_points();
    V = NULL;
    XV = NULL;
    Perm = perm;
  }
  virtual ~Sobol()
  {
    if (V != NULL) {
      for (unsigned j=1;j<=D;j++) delete [] V[j];
      delete [] V;
      V = NULL;
    }
    if (XV != NULL) delete [] XV;
  }

  unsigned skip(unsigned n)
  {
    // the number of points skipped is the largest power of 2 smaller than n
    // unsigned s = (unsigned)floor(log((double)n)/log(2.0));
    // return pow2(s);
    return 0;    // no skipping at the moment
  }
  
  static unsigned **get_directions(unsigned n, unsigned d, unsigned w, unsigned *perm = NULL)
  {
    ostringstream f;
    f << w << ends;
    ifstream infile(f.str().c_str(),ios::in);
    if (!infile) {
      cout << "WARNING: No matching file " << f.str().c_str() << endl;
      exit(1);
    }
    char buffer[1000];
    infile.getline(buffer,1000,'\n');

    unsigned L = (unsigned)ceil(log((double)n)/log(2.0));  // max number of bits needed 
    unsigned **v = new unsigned * [d+1]; 
    // first index is dimension 1 <= j <= d
    // second index is bits     1 <= i <= L
    // These v's are different from the Algorithm. They have been scaled by pow(2,32)
    
    v[1] = new unsigned [L+1];
    for (unsigned i=1;i<=L;i++) v[1][i] = 1 << (32-i);  // divide by pow(2,i) then scale by pow(2,32)
    for (unsigned j=2;j<=d;j++) {
      // read in parameters from file
      unsigned i, s, a;
      infile >> i >> s >> a;
      unsigned *m = new unsigned [s+1];
      for (unsigned i=1;i<=s;i++) infile >> m[i];
      
      // Evaluate v[j][1] to v[j][L]
      v[j] = new unsigned [L+1];
      if (L <= s) {
	for (unsigned i=1;i<=L;i++) v[j][i] = m[i] << (32-i);  // divide by pow(2,i) then scale by pow(2,32)
      }
      else {
	for (unsigned i=1;i<=s;i++) v[j][i] = m[i] << (32-i);  // divide by pow(2,i) then scale by pow(2,32)
	for (unsigned i=s+1;i<=L;i++) {
	  v[j][i] = v[j][i-s] ^ (v[j][i-s] >> s);  // divide by pow(2,s)
	  for (unsigned k=1;k<=s-1;k++)  
	    v[j][i] ^= (((a >> (s-1-k)) & 1) * v[j][i-k]); // A[k] = ((a >> (s-1-k)) & 1);  
	}
      }
      delete [] m;
    }

    // permute if necessary
    if (perm != NULL) {
      unsigned **tmp = v;
      v = new unsigned * [d+1];
      for (unsigned j=1;j<=d;j++) v[j] = tmp[perm[j]];
      delete [] tmp;
    }

    return v;
  }
  
  virtual void init_point()
  { 
    if (V == NULL) V = get_directions(N,D,W,Perm);
    XV = new unsigned [D+1];
    for (unsigned j=1;j<=D;j++) XV[j] = 0;
    counter = 0;
  }
  virtual double *next_point()
  { 
    double *X = new double [D+1];
    if (counter == 0) 
      for (unsigned j=1;j<=D;j++) X[j] = 0;
    
    else {
      unsigned C=1;
      unsigned maskC = 1;      
      unsigned value = counter-1;
      while (value & maskC) {
	value >>= 1;
	C++;
      }
      double scale = pow(2.0,32);
      for (unsigned j=1;j<=D;j++) {
	XV[j] = XV[j] ^ V[j][C];  // update XV
   	X[j] = (double)XV[j]/scale;
      }
    }

    counter++;
    return X;
  }

  virtual void generate_points()
  {
    // no permutation possible

    // ------------- initialise ------------------------------------------------------------------------
    ostringstream f;
    f  << "../shares/datafiles/" << W << ends;
    ifstream infile(f.str().c_str(),ios::in);
    if (!infile) {
      cout << "WARNING: No matching file " << f.str().c_str() << endl;
      exit(1);
    }
    char buffer[1000];
    infile.getline(buffer,1000,'\n');

    // Work out L = max number of bits needed 
    unsigned L = (unsigned)ceil(log((double)N)/log(2.0));  // here "log" is the natural logorithm

    // Work out C[1] to C[N-1]
    // C[i] = index from the right of the first zero bit of i = unsigned 1,2,3,..,L 
    unsigned *C = new unsigned [N];
    unsigned maskC = 1;
    for (unsigned i=1;i<=N-1;i++) {
      C[i] = 1;
      unsigned value = i;
      while (value & maskC) {
	value >>= 1;
	C[i]++;
      }
    };

    // Initialise X[i][j]
    // These X's are different from the Algorithm. They have been scaled by pow(2,32).
    //   ie. x_{i,j} = X[i][j]/pow(2,32) 
    unsigned **X = new unsigned * [N];
    for (unsigned i=0;i<=N-1;i++) {
      X[i] = new unsigned [D+1];
      for (unsigned j=1;j<=D;j++) X[i][j] = 0;
    }


    // ------------------ compute the first dimension -------------------------------------------------

    // Evaluate V[1] to V[L] assuming all m's are 1
    // These V's are different from the Algorithm. They have been scaled by pow(2,32)
    unsigned *V = new unsigned [L+1]; 
    for (unsigned i=1;i<=L;i++)
      V[i] = 1 << (32-i);  // divide by pow(2,i) then scale by pow(2,32)
    
    // Evalulate X[1][1] to X[N-1][1]
    // These X's are different from the Algorithm. They have been scaled by pow(2,32)
    X[1][1] = V[1];
    for (unsigned i=2;i<=N-1;i++)
      X[i][1] = X[i-1][1] ^ V[C[i-1]];

    // Clean up
    delete [] V;


    // ------------------ compute the remaining dimensions -------------------------------------------------
    for (unsigned d=2;d<=D;d++) {
 
      // Read in parameters from file 
      unsigned s;
      unsigned a;
      infile >> d >> s >> a;
      unsigned *m = new unsigned [s+1];
      for (unsigned i=1;i<=s;i++) infile >> m[i];
      
      // Work out A[1] to A[s-1]
      // A[k] = unsigned 0 or 1 
      unsigned *A = new unsigned [s];
      unsigned maskA = 1 << (s-2);
      for (unsigned k=1;k<=s-1;k++) {
	A[k] = (a & maskA ? 1 : 0);
	maskA >>= 1;
      }
    
      // Evaluate V[1] to V[L]
      //   These V's are different from the Algorithm. They have been scaled by pow(2,32) 
      unsigned *V = new unsigned [L+1];
      if (L <= s) {
	for (unsigned i=1;i<=L;i++) 
	  V[i] = m[i] << (32-i);          // divide by pow(2,i) then scale by pow(2,32)
      }
      else {
	for (unsigned i=1;i<=s;i++) 
	  V[i] = m[i] << (32-i);          // divide by pow(2,i) then scale by pow(2,32)
	for (unsigned i=s+1;i<=L;i++) {
	  V[i] = V[i-s] ^ (V[i-s] >> s);  // divide by pow(2,s)
	  for (unsigned k=1;k<=s-1;k++) V[i] ^= A[k] * V[i-k];  
	}
      }
      
      // Evalulate X[1][d] to X[N-1][d]
      // These X's are different from the Algorithm. They have been scaled by pow(2,32) 
      X[1][d] = V[1];
      for (unsigned i=2;i<=N-1;i++)
	X[i][d] = X[i-1][d] ^ V[C[i-1]];
      
      // Clean up
      delete [] m;
      delete [] A;
      delete [] V;
    }

    // ----------------------------- rescale points and clean up ------------------------------------------
    double scale = pow(2.0,32);
    PNTS = new double * [maxN];
    for (unsigned i=0;i<=maxN-1;i++) {
      PNTS[i] = new double [D+1];
      for (unsigned j=1;j<=D;j++)
	PNTS[i][j] = (double)X[i+skip(maxN)][j]/scale;
    }

    delete [] C;
    for (unsigned i=0;i<=N-1;i++) delete [] X[i];
    delete [] X;
  }
};


#endif
