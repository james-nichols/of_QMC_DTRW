#pragma once

#include "ofMain.h"
#include "ofMath.h"
#include "ofSerial.h"
#include "ofEvents.h"
#include "ofThread.h"
#include "ofxGui.h"

#include "qmc.h"

#define H_MARGIN 100
#define W_MARGIN 100
#define TH_MARGIN 20
#define TW_MARGIN 20

#define FRAME_RATE 50

#define _T 0.2
#define _L 1.0
#define _DX 0.01

#define SCALE 20.0

// -----------------------------------------------------------------------------------------
// Pseudo-random number generator
// -----------------------------------------------------------------------------------------
static double ran2local(long *idum) 
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


class dtrwThread : public ofThread
{
  public:

    double* U;
    int n;

    double* K;
    double r;
    double alpha;
    int nx;
    int nt;
    double dX;
    double dT;

    bool centred;

    /// Create a ThreadedObject and initialize the member
    /// variable in an initialization list.
    dtrwThread(int _nx, int _nt, double _dX, double _dT, double _alpha, double _r, bool _centred) 
    {
        nx=_nx; nt=_nt; dX=_dX; dT=_dT; alpha=_alpha; r=_r; centred=_centred;
        n = 0;

        U = new double[nx * nt];
        K = new double[nt];
        
        for (int i=0; i<nx; i++)
            U[i] = 0;
        U[nx/2] = 1.0;

        K[0] = 0;
        K[1] = alpha;
        K[2] = alpha * (alpha - 1.0) / 2.0;
        for (int i=3; i<nt; i++) {
           K[i] = K[i-1] * ((double)i + alpha - 2)/((double)i);
        }
        
    }
    ~dtrwThread() {
       // delete [] U;
       // delete [] K;
    }

    void start();
    void stop();
    void draw();

  protected:

    void threadedFunction();
};

class qmcThread : public ofThread
{
  public:

    double* density;
    int n;

    long seed;

   double* survival; 

    int nx;
    int nt;
    double dX;
    double dT;
    int nMC;
    double alpha;
    double r;
    
    bool centred;
    
    Points* qmc_rule;
    double* shift;

    // Create a ThreadedObject and initialize the member
    // variable in an initialization list.
    qmcThread(int _nx, int _nt, double _dX, double _dT, int _nMC, double _alpha, double _r, bool _centred) 
    {
        nx=_nx; nt=_nt; dX=_dX; dT=_dT; nMC=_nMC; alpha=_alpha; r=_r; centred=_centred;
        n = 0;
        seed = -2;
        
        density = new double[nx];
        for (int i=0; i<nx; i++)
            density[i] = 0.0;
        
        survival = new double[nt];
        survival[0] = 1.0 - alpha;
        for (int i=1; i<nt; i++)
            survival[i] = survival[i-1] * (1. - alpha/((double)(i+1)));
        
        qmc_rule = new Lattice(1048576, nt, 39102);
        qmc_rule->init_point();
         
        shift = Points::get_points(3600,seed); // the seed is changed after this
    }
    ~qmcThread() {
        //delete [] density;
        //delete [] survival;
    }
        
    // Start the thread.
    void start();
    void stop();
    void draw();

  protected:
    void threadedFunction();
    int inv_cum_dist(double x, double* dist, int len) {
        int n = 1;
        while (x < dist[n-1] && n < len) n++;
        return n;
    }

};

class mcThread : public ofThread
{
  public:

    double* density;
    int n;

    long seed;

   double* survival; 

    int nx;
    int nt;
    double dX;
    double dT;
    int nMC;
    double alpha;
    double r;

    bool centred;
    
    // Create a ThreadedObject and initialize the member
    // variable in an initialization list.
    mcThread(int _nx, int _nt, double _dX, double _dT, int _nMC, double _alpha, double _r, bool _centred) 
    {
        nx=_nx; nt=_nt; dX=_dX; dT=_dT; nMC=_nMC; alpha=_alpha; r=_r; centred=_centred;
        n = 0;
        seed = -2;
        
        density = new double[nx];
        for (int i=0; i<nx; i++)
            density[i] = 0.0;
        
        survival = new double[nt];
        survival[0] = 1.0 - alpha;
        for (int i=1; i<nt; i++)
            survival[i] = survival[i-1] * (1. - alpha/((double)(i+1)));;
    }
    ~mcThread() {
        //delete [] density;
        //delete [] survival;
    }
        
    // Start the thread.
    void start();
    void stop();
    void draw();

  protected:
    void threadedFunction();
    int inv_cum_dist(double x, double* dist, int len) {
        int n = 1;
        while (x < dist[n-1] && n < len) n++;
        return n;
    }

};



class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);

        void restartSims();

        // The raw data
        dtrwThread* dtrw_solver;
        mcThread* mc_solver;
        qmcThread* qmc_solver;
    
        double T, L, dT, dX;

        double alpha, D_alpha;
        int nx, nt, nMC;
        double r;

        // Parameters setting GUI
    	ofxPanel gui;
    	ofxFloatSlider T_slider;
    	ofxFloatSlider D_alpha_slider;
    	ofxFloatSlider alpha_slider;
    	ofxIntSlider num_mc;
        ofxButton restart;

        bool bRun;
};
