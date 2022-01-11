// code to solve small electron-phonon problem with two electrons on two sites
// add option to compute spectral function from photoemission formula
// with retarded spectrum
// Phonon-photon polariton + nonlinear Q^2*U coupling relevant for organic kappa salts
// Optimized code with new BLAS routines. Work for F != 0
// Particle-Hole Symmetry 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <assert.h>
#include <math.h>
#include "matrix.h"
#include <mkl.h>
#ifndef NO_MPI
    #include <mpi.h>
#endif


using namespace std;

typedef complex<double> cdouble;
typedef matrix<cdouble> cmat;
cdouble II(0,1);

#define NSITES        2
#define NPHON         30 // Number of phonon
#define NHILB         3600 // NHILBPHON * 4 (Total Hilbert space)
#define NHILBPHON     900 // NPHON^2 (Number of phonon * Number of photon)
#define NHILBREDUCED  1800 // NHILB / 2 (or NHILBPHON*2)

// here we use same bosonic Hilbert space size for photons and phonons!

#define OMEGA       0.20 // phonon freq
#define QUARTIC     0.00 // quartic term for phonon potential
#define OMEGAPHOT   0.20 // photon freq

#define WPLASMA     0.10  // plasma frequency for phonon-photon coupling
#define G2          -0.100 // Q^2*U nonlinear electron-phonon coupling; this should be around 0.1*OMEGA

#define JHOP        0.25
#define U1          1.25   // Hubbard U

#define OMEGAPUMP   0.55
#define FMAX        0.00
#define CWDRIVE     // if defined use cw drive; else use pulse
#define FWIDTH      3.0
#define FCENTER     15.0
#define TMAX        50.
#define NUMT        100

// for spectrum
#define TARPES      25.0
#define SIGARPES    8.0  // sigma for ARPES = probe duration
#define NTARPES     24  // 2*NTARPES+1 = number of time points for ARPES around TARPES = 6 sigma window
#define NOMEGAARPES 1600
#define OMEGAMIN    -6.
#define OMEGAMAX    0.

// Constants that CBLAS routines need regularly
const cdouble ZERO = 0.0;
const cdouble ONE = 1.0;

void invert(cdouble* G, MKL_INT msize)
{
    MKL_INT LWORK = msize, INFO=0;
    cdouble work[NHILB]; // use maximal size here
    assert(NHILB>=msize);
    MKL_INT permutations[msize];
    
    zgetrf_(&msize, &msize,reinterpret_cast<MKL_Complex16*>(G), &msize, permutations, &INFO);
    if(INFO != 0)
    {
        cout << "Complex matrix inverse: Error at zgetrf. INFO: " << INFO<< endl;
    }
    
    zgetri_(&msize, reinterpret_cast<MKL_Complex16*>(G), &msize, permutations,reinterpret_cast<MKL_Complex16*>(work), &LWORK, &INFO);
    if(INFO != 0)
    {
        cout << "Complex matrix inverse: Error at zgetri. INFO: " << INFO<< endl;
    }
    
}

void diagonalize(cmat &H, vector<double> &evals, int N)
{
    // MKL_INT     matsize = N;
    // MKL_INT     lwork = 2*(2*N-1);
    // double  rwork[3*N-2];
    // cdouble work[2*(2*NHILB-1)]; // use maximal size here
    char    jobz = 'V';
    char    uplo = 'U';
    // MKL_INT info;

    lapack_int info = LAPACKE_zheev(
        LAPACK_COL_MAJOR,
        jobz, uplo, N,
        reinterpret_cast<MKL_Complex16*>(H.ptr()),
        N, &evals[0]
    );

    
	// zheev_(&jobz, &uplo, &matsize,reinterpret_cast<MKL_Complex16*>(H.ptr()), &matsize, &evals[0], reinterpret_cast<MKL_Complex16*>(&work[0]), &lwork, &rwork[0], &info);
	assert(!info);
}

void times(cmat &A, cmat &B, cmat &C, int N)
{
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                N, N, N,
                &ONE,
                &A(0, 0), N,
                &B(0, 0), N,
                &ZERO,
                &C(0, 0), N);
   
}

void set_U_step(cmat &U, cmat &A1, cmat &A2, vector<double> &evals, int N, double dtstep)
{
    cmat tmp1(N, N);
    cmat tmp2(N, N);
    
    cmat B1(N, N);
    cmat B2(N, N);

    cmat B_work(N, N);
    vector<cdouble> exp_ite(N); 
    
    double c1 = (3.-2.*sqrt(3.))/12.;
    double c2 = (3.+2.*sqrt(3.))/12.;
    
    for(int i=0; i < N; i++)
        for(int j=0; j < N; j++)
        {
            B1(i,j) = c1*A1(i,j) + c2*A2(i,j);
            B2(i,j) = c2*A1(i,j) + c1*A2(i,j);
        }
    
    
    diagonalize(B1, evals, N);
    for(int k = 0; k < N; k++) {
        exp_ite[k] = exp(II * dtstep * evals[k]);
    }
    for (int j = 0; j < N; j++) {
        for(int k = 0; k < N; k++) {
            B_work(j, k) = exp_ite[k] * B1(j, k);
        }
    }
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, N, N, N,
                &ONE, B1.ptr(), N, B_work.ptr(), N, &ZERO, tmp1.ptr(), N);


    diagonalize(B2, evals, N);
    for(int k = 0; k < N; k++) {
        exp_ite[k] = exp(II * dtstep * evals[k]);
    }
    for (int j = 0; j < N; j++) {
        for(int k = 0; k < N; k++) {
            B_work(j, k) = exp_ite[k] * B2(j, k);
        }
    }
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, N, N, N,
                &ONE, B2.ptr(), N, B_work.ptr(), N, &ZERO, tmp2.ptr(), N);


    times(tmp1, tmp2, U, N);
}

inline double shape(double time)
{
    return exp(-(time-FCENTER)*(time-FCENTER)/(2.*FWIDTH*FWIDTH));
}

inline double arpesshape(double time1, double time2, double time0, double sigma)
{
    return exp(-(time1-time0)*(time1-time0)/(2.*sigma*sigma)) * exp(-(time2-time0)*(time2-time0)/(2.*sigma*sigma)) / (sigma*sqrt(M_PI)*2.*M_PI); // only one sigma in denominator for norm
}

int ed_ssh_main(int argc, char * argv[])
{
    //************** MPI INIT ***************************
  	int numprocs=1, myrank=0, namelen;

  //  ofstream csvtime;
  //  ofstream csvspin2;
 //   ofstream csvarpes2;
    ofstream csvGS_carac;
  //  ofstream outspecweight;

  //  ofstream csvPhot_distrib;
  //  ofstream csvPhon_distrib;
  //  ofstream csvPhotPhon_distrib;

    /* if(myrank==0)
    {
        csvspin2.open("Spin_PHS_MN40_J0_U125_Omega02_G002.csv");
        csvspin2 << setw(4) << setprecision(15);
        csvspin2 << "# omega" << ';' << "spin" << endl;
    }*/

    /* if(myrank==0)
    {
        csvarpes2.open("ARPES_PHS_MN60_J0_U125_Omega02_G005_gam1_WP12_SIG8.csv");
        csvarpes2 << setw(4) << setprecision(15);
        csvarpes2 << "# omega" << ';' << "Arpes occ" << ';' << "Arpes tot" << endl;
    }*/

    /*if(myrank==0)
    {
        csvtime.open("time_PHS_MN40_J0_U125_Omega02_G002.csv");
        csvtime << setw(4) << setprecision(15);
        csvtime << "# time" << ';' << "ekin" << ';' << "epot" << ';' << "ekin+epot" << ';' << "ephon" << ';' << "etot" << ';' << "x(t)" << ';' << "field" << ';' << "docc" << endl;
    }
*/
     if(myrank==0)
    {
        csvGS_carac.open("Carac_PHS_MN30_Omega02_G008_WP08.csv");
        csvGS_carac << setw(4) << setprecision(15);
        csvGS_carac << "#WP" << ';' << "g coupling" << ';' <<"Energy GS" << ';' << "Doublon" << ';' << "Ueff" << ';' << "nphot"<< ';' << "nphon"<< ';' << "Variance phonon" << ';' << "OMEGA eff (phon)" << endl;
    }/*
    if(myrank==0)
    {
        outspecweight.open("Spec_weight_PHS_MN60_J0_U125_Omega02_G005_gam1_WP12_SIG8.dat");
        outspecweight << setw(4) << setprecision(15);
        outspecweight << "# Spectral weight  " << endl;
    }
*/
  /*        if(myrank==0)
    {
        csvPhot_distrib.open("Photdistrib_PHS_MN40_J0_U125_Omega02_G002.csv");
        csvPhot_distrib << setw(4) << setprecision(15);
        csvPhot_distrib << "#NPHON" << ';' << "g coupling" << ';' <<"wP" << ';' << "Photon distribution" << ';' << "Photon distribution coherent state " << ';' << "Number of phot for distrib" << ';' << "Photon distribution squeezed state" << ';' << "Number of phot for squeezed" << endl;
    }


          if(myrank==0)
    {
        csvPhon_distrib.open("Phondistrib_PHS_MN40_J0_U125_Omega02_G002.csv");
        csvPhon_distrib << setw(4) << setprecision(15);
        csvPhon_distrib << "#NPHON" << ';' << "g coupling" << ';' <<"wP" << ';' << "Phonon distribution" << ';' << "Phonon distribution coherent state " << ';' << "Number of phon for distrib" << ';' << "Phonon distribution squeezed state" << ';' << "Number of phon for squeezed" << endl;
    }

               if(myrank==0)
    {
        csvPhotPhon_distrib.open("TMdistrib_PHS_MN40_J0_U125_Omega02_G002.csv");
        csvPhotPhon_distrib << setw(4) << setprecision(15);
        csvPhotPhon_distrib << "#NPHON" << ';' << "g coupling" << ';' <<"wP" << ';' << "Photon Phonon distribution" << ';' << "Photon Phonon distribution TM squeezed state" << ';' << "Number of phon for squeezed" << endl;
    }

*/
    
#ifndef NO_MPI
  	char processor_name[MPI_MAX_PROCESSOR_NAME];
  	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  	MPI_Get_processor_name(processor_name, &namelen);
    
	MPI_Barrier(MPI_COMM_WORLD);
    
	MPI_Datatype ttype  = MPI_DOUBLE;
    
#endif
	if(myrank==0) cout << "\n\tProgram running on " << numprocs << " processors." << endl;
    
#ifndef NO_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    //************** Build the Hilbert space ***************************
    
    // 2 electrons with spin on 2 sites; I can add Hubbard U if I want
    
    // label all possible states in occupation number representation
    vector<double> nph0(NHILB);
    vector<double> nph1(NHILB);
    vector<double> nup0(NHILB);
    vector<double> ndn0(NHILB);
    vector<double> nup1(NHILB);
    vector<double> ndn1(NHILB);
    vector<double> nel0(NHILB);
    vector<double> nel1(NHILB);
    
    // fermionic ordering to clarify signs: site 0 left, site 1 right; up left, down right
    // electrons |up down ; 0 0>
    for(int n0=0; n0<NPHON; ++n0)
        for(int n1=0; n1<NPHON; ++n1)
        {
            nph0[n0 + n1*NPHON] = n0;
            nph1[n0 + n1*NPHON] = n1;
            nup0[n0 + n1*NPHON] = 1;
            ndn0[n0 + n1*NPHON] = 1;
            nup1[n0 + n1*NPHON] = 0;
            ndn1[n0 + n1*NPHON] = 0;
            nel0[n0 + n1*NPHON] = 2;
            nel1[n0 + n1*NPHON] = 0;
        }
    
    // electrons |up 0 ; 0 down>
    for(int n0=0; n0<NPHON; ++n0)
        for(int n1=0; n1<NPHON; ++n1)
        {
            nph0[NHILBPHON + n0 + n1*NPHON] = n0;
            nph1[NHILBPHON + n0 + n1*NPHON] = n1;
            nup0[NHILBPHON + n0 + n1*NPHON] = 1;
            ndn0[NHILBPHON + n0 + n1*NPHON] = 0;
            nup1[NHILBPHON + n0 + n1*NPHON] = 0;
            ndn1[NHILBPHON + n0 + n1*NPHON] = 1;
            nel0[NHILBPHON + n0 + n1*NPHON] = 1;
            nel1[NHILBPHON + n0 + n1*NPHON] = 1;
        }
    
    // electrons |0 down ; up 0>
    for(int n0=0; n0<NPHON; ++n0)
        for(int n1=0; n1<NPHON; ++n1)
        {
            nph0[2*NHILBPHON + n0 + n1*NPHON] = n0;
            nph1[2*NHILBPHON + n0 + n1*NPHON] = n1;
            nup0[2*NHILBPHON + n0 + n1*NPHON] = 0;
            ndn0[2*NHILBPHON + n0 + n1*NPHON] = 1;
            nup1[2*NHILBPHON + n0 + n1*NPHON] = 1;
            ndn1[2*NHILBPHON + n0 + n1*NPHON] = 0;
            nel0[2*NHILBPHON + n0 + n1*NPHON] = 1;
            nel1[2*NHILBPHON + n0 + n1*NPHON] = 1;
        }
    
    // electrons |0 0 ; up down>
    for(int n0=0; n0<NPHON; ++n0)
        for(int n1=0; n1<NPHON; ++n1)
        {
            nph0[3*NHILBPHON + n0 + n1*NPHON] = n0;
            nph1[3*NHILBPHON + n0 + n1*NPHON] = n1;
            nup0[3*NHILBPHON + n0 + n1*NPHON] = 0;
            ndn0[3*NHILBPHON + n0 + n1*NPHON] = 0;
            nup1[3*NHILBPHON + n0 + n1*NPHON] = 1;
            ndn1[3*NHILBPHON + n0 + n1*NPHON] = 1;
            nel0[3*NHILBPHON + n0 + n1*NPHON] = 0;
            nel1[3*NHILBPHON + n0 + n1*NPHON] = 2;
        }
    
    
    // build static Hamiltonian
    cmat Ham0(NHILB, NHILB);
    cmat Hamkin(NHILB, NHILB);
    cmat Hamu(NHILB, NHILB);
    cmat Hamdocc(NHILB, NHILB);
    cmat Hamdocc_PHS(NHILB, NHILB);
    //cmat Hamu0(NHILB, NHILB); // initial U0 for U quench
    cmat Hamph(NHILB, NHILB);
    cmat Hamx(NHILB, NHILB); // X term for photons!, to be multiplied by f(t)
    //cmat Hamg1(NHILB, NHILB); // linear interaction
    cmat Hamg2(NHILB, NHILB); // nonlinear interaction
    
    cmat Hamxsquared(NHILB, NHILB); //phonon (b+bdag)^2
    cmat Hamnphot(NHILB, NHILB); //photon number
    cmat Hamnphon(NHILB, NHILB); //phonon number

    cmat Sz0(NHILB, NHILB); // for spin correlations
    
    for(int i=0; i<NHILB; ++i)
    {
        Sz0(i,i) = 0.5*(nup0[i]-ndn0[i]);
        // TESTING:
        
        int ielec = int(i/NHILBPHON); // whether we are in sector 0, 1, 2, 3
        Hamph(i,i) = OMEGAPHOT * nph0[i]+ OMEGA* nph1[i] + WPLASMA*WPLASMA/(2.*OMEGAPHOT) * nph0[i]; // automatically diagonal in electrons, since we are in diagonal part
        Hamxsquared(i,i) = (2.*nph1[i]+1.) ; // from 2n+1 term in (b+bdag)^2

       // Hamg2(i,i) = G2 * ((nup0[i] * ndn0[i] + nup1[i] * ndn1[i])  * (2.*nph1[i]+1.)); // from 2n+1 term in (b+bdag)^2
        // (Below for Particle-Hole Symmetry)

        Hamnphot(i,i) = nph0[i];
        Hamnphon(i,i) = nph1[i];
        
        Hamg2(i,i) = -G2/2. *  (2.*nph1[i]+1.); // Offset for PHS. From 2n+1 term in (b+bdag)^2
        Hamu(i,i) = -U1/2;  // Offset for Particle-Hole transformation
        Hamdocc_PHS(i,i) = -1/2; 

        if(ielec == 0 || ielec == 3)
        {
            Hamu(i,i) += U1;
            Hamdocc(i,i) = 1.;
            Hamdocc_PHS(i,i) += 1.;
            Hamg2(i,i) += G2 * (2.*nph1[i]+1.); // from 2n+1 term in (b+bdag)^2
        }
        
        for(int j=0; j<NHILB; ++j)
        {
            int jelec = int(j/NHILBPHON);
            
            // add here the phonon-photon coupling terms G1*(adag * b + bdag * a)
            // diagonal in electronic sector
            //
            if(ielec == jelec)
            {
                if(nph0[i]==nph0[j]+1 && nph1[i]==nph1[j]-1) // this is the adag*b term
                    Hamph(i,j) = II*WPLASMA/2. * sqrt(OMEGA/OMEGAPHOT) * sqrt(nph0[j]+1) * sqrt(nph1[j]); // larger of the two numbers for adag, smaller of the two numbers for b
                if(nph0[i]==nph0[j]-1 && nph1[i]==nph1[j]+1) // this is the a*bdag term
                    Hamph(i,j) = -II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT) * sqrt(nph0[j]) * sqrt(nph1[j]+1);
                
                // now counter-rotating terms
                if(nph0[i]==nph0[j]-1 && nph1[i]==nph1[j]-1) // this is the a*b term
                    Hamph(i,j) = II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT) * sqrt(nph0[j]) * sqrt(nph1[j]);
                
                if(nph0[i]==nph0[j]+1 && nph1[i]==nph1[j]+1) // this is the adag*bdag term
                    Hamph(i,j) = -II*WPLASMA/2. * sqrt(OMEGA/OMEGAPHOT) * sqrt(nph0[j]+1) * sqrt(nph1[j]+1);
                
                if(nph0[i]==nph0[j]+2 && nph1[i]==nph1[j]) // wplasma terms to photon
                    Hamph(i,j) = WPLASMA*WPLASMA/(4.*OMEGAPHOT) * sqrt((nph0[j]+2)*(nph0[j]+1));
                
                if(nph0[i]==nph0[j]-2 && nph1[i]==nph1[j]) // wplasma terms to photon
                    Hamph(i,j) = WPLASMA*WPLASMA/(4.*OMEGAPHOT) * sqrt((nph0[j])*(nph0[j]-1));
                
                
                // add here quartic terms for phonon
                if(nph0[i]==nph0[j])
                {
                    if(nph1[i]==nph1[j]-4)
                        Hamph(i,j) += QUARTIC * sqrt(nph1[j]*(nph1[j]-1)*(nph1[j]-2)*(nph1[j]-3));
                    if(nph1[i]==nph1[j]+4)
                        Hamph(i,j) += QUARTIC * sqrt((nph1[j]+1)*(nph1[j]+2)*(nph1[j]+3)*(nph1[j]+4));
                    if(nph1[i]==nph1[j]+2)
                        Hamph(i,j) += 4.*QUARTIC * sqrt(nph1[j]*nph1[j]*(nph1[j]+1)*(nph1[j]+2));
                    if(nph1[i]==nph1[j]-2)
                        Hamph(i,j) += 4.*QUARTIC * sqrt(nph1[j]*(nph1[j]-1)*(nph1[j]-2)*(nph1[j]-2));
                    if(nph1[i]==nph1[j])
                        Hamph(i,j) += 6.*QUARTIC * sqrt(nph1[j]*(nph1[j]-1)*(nph1[j]-1)*nph1[j]);
		    if(nph1[i]==nph1[j]+2)
                        Hamph(i,j) += 6.*QUARTIC * sqrt((nph1[j]+1)*(nph1[j]+2));
                    if(nph1[i]==nph1[j]-2)
                        Hamph(i,j) += 6.*QUARTIC * sqrt((nph1[j])*(nph1[j]-1));
                    if(nph1[i]==nph1[j])
                        Hamph(i,j) += 12.*QUARTIC * nph1[j];
                }
                
            }
            
            if( ((ielec == 1 && jelec == 0) || (ielec == 0 && jelec == 1) || (ielec == 1 && jelec == 3) || (ielec == 3 && jelec == 1)) && nph0[i]==nph0[j] && nph1[i]==nph1[j])
            {
                Hamkin(i,j) = -JHOP;
            }
            
            if( ((ielec == 2 && jelec == 0) || (ielec == 0 && jelec == 2) || (ielec == 2 && jelec == 3) || (ielec == 3 && jelec == 2)) && nph0[i]==nph0[j] && nph1[i]==nph1[j])
            {
                Hamkin(i,j) = JHOP;
            }
            
            
            if(nph0[i]==nph0[j]+1 && nph1[i]==nph1[j] && ielec == jelec)
            {  
                Hamx(i,j) = sqrt(nph0[j]+1);
            }
            
            if(nph0[i]==nph0[j]-1 && nph1[i]==nph1[j]  && ielec == jelec)
            {
                Hamx(i,j) = sqrt(nph0[j]); // automatically 0 for nph=0 state
            }
            
            if(nph1[i]==nph1[j]+2 && nph0[i]==nph0[j] && ielec == jelec) // phonon is nph1; photon is nph0!!!
            {
               // Hamg2(i,j) = G2 * (nup0[j] * ndn0[j] + nup1[j] * ndn1[j]) * sqrt((nph1[j]+2)*(nph1[j]+1)); // No Particle-Hole Symmetry
                Hamg2(i,j) -= G2/2. *  sqrt((nph1[j]+2)*(nph1[j]+1)); // PHS offset
                Hamxsquared(i,j) = sqrt((nph1[j]+2)*(nph1[j]+1));
            }
            
            if(nph1[i]==nph1[j]-2 && nph0[i]==nph0[j]  && ielec == jelec)
            {
                //Hamg2(i,j) = G2 * (nup0[j] * ndn0[j] + nup1[j] * ndn1[j]) * sqrt((nph1[j])*(nph1[j]-1)); // automatically 0 for nph=0 state // No Particle-Hole Symmetry
                Hamg2(i,j) -= G2/2. * sqrt((nph1[j])*(nph1[j]-1)); // PHS offset
                Hamxsquared(i,j) = sqrt((nph1[j])*(nph1[j]-1));
            }

             // G2 term for PHS symmetry


            if(nph1[i]==nph1[j]+2 && nph0[i]==nph0[j] && ielec == jelec && ielec==0) // phonon is nph1; photon is nph0!!!
            {
                Hamg2(i,j) += G2 * sqrt((nph1[j]+2)*(nph1[j]+1));
            }
              if(nph1[i]==nph1[j]+2 && nph0[i]==nph0[j] && ielec == jelec && ielec==3) // phonon is nph1; photon is nph0!!!
            {
                Hamg2(i,j) += G2 *  sqrt((nph1[j]+2)*(nph1[j]+1));
            }

            if(nph1[i]==nph1[j]-2 && nph0[i]==nph0[j]  && ielec == jelec && ielec==0)
            {
                Hamg2(i,j) += G2 *  sqrt((nph1[j])*(nph1[j]-1)); // automatically 0 for nph=0 state
            }
            if(nph1[i]==nph1[j]-2 && nph0[i]==nph0[j]  && ielec == jelec && ielec==3)
            {
                Hamg2(i,j) += G2 * sqrt((nph1[j])*(nph1[j]-1)); // automatically 0 for nph=0 state
            }

            
        }
        
    }
    
    for(int i=0; i<NHILB; ++i)
        for(int j=0; j<NHILB; ++j)
        {
            Ham0(i,j) = Hamph(i,j) + Hamkin(i,j) + Hamu(i,j) + Hamg2(i,j);
           // Ham0(i,j) =  Hamkin(i,j) + Hamu(i,j);
        }
    
    // REDUCED HILBERT SPACE FOR N-1 PARTICLES (1 electron, spin down, to fix it)
    // ADDED HILBERT SPACE FOR N+1 PARTICLES (2 up spins, 1 down spin remaining)
    /*
    vector<double> nel0red(NHILBREDUCED);
    vector<double> nel1red(NHILBREDUCED);
    
    vector<double> nel0added(NHILBREDUCED);
    vector<double> nel1added(NHILBREDUCED);
    
    // fermionic ordering to clarify signs: site 0 left, site 1 right; up left, down right
    // electrons |0 down ; 0 0>
    for(int n0=0; n0<NPHON; ++n0)
        for(int n1=0; n1<NPHON; ++n1)
        {
            nel0red[n0 + n1*NPHON] = 1;
            nel1red[n0 + n1*NPHON] = 0;
        }
    
    // electrons |0 0 ; 0 down>
    for(int n0=0; n0<NPHON; ++n0)
        for(int n1=0; n1<NPHON; ++n1)
        {
            nel0red[NHILBPHON + n0 + n1*NPHON] = 0;
            nel1red[NHILBPHON + n0 + n1*NPHON] = 1;
        }
    
    // fermionic ordering to clarify signs: site 0 left, site 1 right; up left, down right
    // electrons |up down ; up 0>
    for(int n0=0; n0<NPHON; ++n0)
        for(int n1=0; n1<NPHON; ++n1)
        {
            nel0added[n0 + n1*NPHON] = 2;
            nel1added[n0 + n1*NPHON] = 1;
        }
    
    // electrons |up 0 ; up down>
    for(int n0=0; n0<NPHON; ++n0)
        for(int n1=0; n1<NPHON; ++n1)
        {
            nel0added[NHILBPHON + n0 + n1*NPHON] = 1;
            nel1added[NHILBPHON + n0 + n1*NPHON] = 2;
        }
    
    cmat Ham0red(NHILBREDUCED, NHILBREDUCED);
    cmat Hamxred(NHILBREDUCED, NHILBREDUCED);
    
    // need to add Hubbard U term in 3-electron sector!
    
    cmat Ham0added(NHILBREDUCED, NHILBREDUCED);
    cmat Hamxadded(NHILBREDUCED, NHILBREDUCED);
    
    for(int i=0; i<NHILBREDUCED; ++i)
    {
        
        int ielec = int(i/NHILBPHON); // whether we are in sector 0, 1
        Ham0red(i,i) += OMEGAPHOT * nph0[i] + OMEGA* nph1[i] + WPLASMA*WPLASMA/(2.*OMEGAPHOT) * nph0[i]; // automatically diagonal in electrons, since we are in diagonal part; phonon occupations are independent of electron sector and also work in reduced space! :-)
        Ham0added(i,i) += OMEGAPHOT * nph0[i] +  OMEGA * nph1[i] + WPLASMA*WPLASMA/(2.*OMEGAPHOT) * nph0[i];
        
        //Ham0red(i,i) += G2 * (nel0red[i] * (2.*nph0[i]+1.) + nel1red[i] * (2.*nph1[i]+1.)) ; // from 2n+1 term in (b+bdag)^2
        // G2 term is not active in reduced sector because there is no double occupancy
        // Ham0added(i,i) += G2 * (2.*nph1[i]+1.) ; // all electronic states in 3-electron-sector have one double occupancy
        //Ham0added(i,i) += U1; // all electronic states in 3-electron-sector have one double occupancy
       
       // No U or G2 term in Added/Reduced sectors for Particle-Hole Symmetry
        
        for(int j=0; j<NHILBREDUCED; ++j)
        {
            int jelec = int(j/NHILBPHON); // whether we are in sector 0, 1
            
            // add here the phonon-photon coupling terms G1*(adag * b + bdag * a)
            // diagonal in electronic sector
            //
            if(ielec == jelec)
            {
                if(nph0[i]==nph0[j]+1 && nph1[i]==nph1[j]-1) // this is the adag*b term
                {
                    Ham0red(i,j) += II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT)* sqrt(nph0[j]+1) * sqrt(nph1[j]);
                    Ham0added(i,j) += II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT)* sqrt(nph0[j]+1) * sqrt(nph1[j]); // larger of the two numbers for adag, smaller of the two numbers for b
                }
                if(nph0[i]==nph0[j]-1 && nph1[i]==nph1[j]+1) // this is the a*bdag term
                {
                    Ham0red(i,j) -= II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT)* sqrt(nph0[j]) * sqrt(nph1[j]+1);
                    Ham0added(i,j) -= II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT)* sqrt(nph0[j]) * sqrt(nph1[j]+1);
                }
                
                // now counter-rotating terms
                if(nph0[i]==nph0[j]-1 && nph1[i]==nph1[j]-1) // this is the a*b term
                {
                    Ham0red(i,j) += II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT)* sqrt(nph0[j]) * sqrt(nph1[j]);
                    Ham0added(i,j) += II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT)* sqrt(nph0[j]) * sqrt(nph1[j]);
                }
                
                if(nph0[i]==nph0[j]+1 && nph1[i]==nph1[j]+1) // this is the adag*bdag term
                {
                    Ham0red(i,j) -= II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT)* sqrt(nph0[j]+1) * sqrt(nph1[j]+1);
                    Ham0added(i,j) -= II*WPLASMA/2.  * sqrt(OMEGA/OMEGAPHOT)* sqrt(nph0[j]+1) * sqrt(nph1[j]+1);
                }
                
                if(nph0[i]==nph0[j]+2 && nph1[i]==nph1[j]) // wplasma terms to photon
                {
                    Ham0red(i,j) += WPLASMA*WPLASMA/(4.*OMEGAPHOT) * sqrt((nph0[j]+2)*(nph0[j]+1));
                    Ham0added(i,j) += WPLASMA*WPLASMA/(4.*OMEGAPHOT) * sqrt((nph0[j]+2)*(nph0[j]+1));
                }
                
                if(nph0[i]==nph0[j]-2 && nph1[i]==nph1[j]) // wplasma terms to photon
                {
                    Ham0red(i,j) += WPLASMA*WPLASMA/(4.*OMEGAPHOT) * sqrt((nph0[j])*(nph0[j]-1));
                    Ham0added(i,j) += WPLASMA*WPLASMA/(4.*OMEGAPHOT) * sqrt((nph0[j])*(nph0[j]-1));
                }
                
                // add here quartic terms for phonon
                if(nph0[i]==nph0[j])
                {
                    if(nph1[i]==nph1[j]-4)
                    {
                        Ham0red(i,j) += QUARTIC * sqrt(nph1[j]*(nph1[j]-1)*(nph1[j]-2)*(nph1[j]-3));
                        Ham0added(i,j) += QUARTIC * sqrt(nph1[j]*(nph1[j]-1)*(nph1[j]-2)*(nph1[j]-3));
                    }
                    if(nph1[i]==nph1[j]+4)
                    {
                        Ham0red(i,j) += QUARTIC * sqrt((nph1[j]+1)*(nph1[j]+2)*(nph1[j]+3)*(nph1[j]+4));
                        Ham0added(i,j) += QUARTIC * sqrt((nph1[j]+1)*(nph1[j]+2)*(nph1[j]+3)*(nph1[j]+4));
                    }
                    if(nph1[i]==nph1[j]+2)
                    {
                        Ham0red(i,j) += 4.*QUARTIC * sqrt(nph1[j]*nph1[j]*(nph1[j]+1)*(nph1[j]+2));
                        Ham0added(i,j) += 4.*QUARTIC * sqrt(nph1[j]*nph1[j]*(nph1[j]+1)*(nph1[j]+2));
                    }
                    if(nph1[i]==nph1[j]-2)
                    {
                        Ham0red(i,j) += 4.*QUARTIC * sqrt(nph1[j]*(nph1[j]-1)*(nph1[j]-2)*(nph1[j]-2));
                        Ham0added(i,j) += 4.*QUARTIC * sqrt(nph1[j]*(nph1[j]-1)*(nph1[j]-2)*(nph1[j]-2));
                    }
                    if(nph1[i]==nph1[j])
                    {
                        Ham0red(i,j) += 6.*QUARTIC * sqrt(nph1[j]*(nph1[j]-1)*(nph1[j]-1)*nph1[j]);
                        Ham0added(i,j) += 6.*QUARTIC * sqrt(nph1[j]*(nph1[j]-1)*(nph1[j]-1)*nph1[j]);
                    }

		     if(nph1[i]==nph1[j]+2)
                    {
                            Ham0red(i,j) += 6.*QUARTIC * sqrt((nph1[j]+1)*(nph1[j]+2));
                            Ham0added(i,j) += 6.*QUARTIC * sqrt((nph1[j]+1)*(nph1[j]+2));
                    }
                    if(nph1[i]==nph1[j]-2)
                    {
                            Ham0red(i,j) += 6.*QUARTIC * sqrt((nph1[j])*(nph1[j]-1));
                            Ham0added(i,j) += 6.*QUARTIC * sqrt((nph1[j])*(nph1[j]-1));
                    }
                }
                
            }
            
            if( ((ielec == 1 && jelec == 0) || (ielec == 0 && jelec == 1)) && nph0[i]==nph0[j] && nph1[i]==nph1[j])
            {
                Ham0red(i,j) -= JHOP;
                Ham0added(i,j) += JHOP; // CAREFUL WITH SIGN - HERE WE ARE HOPPING BETWEEN OTHER STATES ...
            }
            
            
            if(nph0[i]==nph0[j]+1 && nph1[i]==nph1[j] && ielec == jelec)
            {
                Hamxred(i,j) = sqrt(nph0[j]+1);
                Hamxadded(i,j) = sqrt(nph0[j]+1);
            }
            
            if(nph0[i]==nph0[j]-1 && nph1[i]==nph1[j]  && ielec == jelec)
            {
                Hamxred(i,j) = sqrt(nph0[j]); // automatically 0 for nph=0 state
                Hamxadded(i,j) = sqrt(nph0[j]); // automatically 0 for nph=0 state
            }
         
         // No g2 nor U terms in added/reduced sectors for Particle-Hole Symmetry
*/
          /*  if(nph1[i]==nph1[j]+2 && nph0[i]==nph0[j]  && ielec == jelec)
            {
                //Ham0red(i,j) += G2 * sqrt((nph1[j]+2)*(nph1[j]+1)); // no double occ in 1-electron sector
                Ham0added(i,j) += G2 * sqrt((nph1[j]+2)*(nph1[j]+1)); // all electronic states in 3-electron-sector have one double occupancy
            }
            
            if(nph1[i]==nph1[j]-2 && nph0[i]==nph0[j]  && ielec == jelec)
            {
                //Ham0red(i,j) += G2 * sqrt((nph1[j])*(nph1[j]-1)); // no double occ in 1-electron sector
                Ham0added(i,j) += G2 * sqrt((nph1[j])*(nph1[j]-1)); // all electronic states in 3-electron-sector have one double occupancy
            }*/
            
            
  //      }
        
  //  }
    
    /*
    
    cmat Annihilator(NHILBREDUCED, NHILB); // operator that produces N-1 particle wave function from N particle wave function, annihilating up-electron on site 0
    for(int i=0; i<NHILBREDUCED; ++i)
    {
        //int ielec = int(i/NHILBPHON); // whether we are in sector 0, 1
        for(int j=0; j<NHILB; ++j)
        {
            Annihilator(i,j) = 0;
            //int jelec = int(j/NHILBPHON); // whether we are in sector 0, 1, 2, 3
            if(i==j)
                Annihilator(i,j) = 1; // really that simple ... note that this is not a diagonal matrix
        }
    }
    
    
    // check sign of hopping for added term above
    // sectors are 21, 12
    // sign is "-" ==> +JHOP
    
    
    cmat Creator(NHILBREDUCED, NHILB); // operator that produces N+1 particle wave function from N particle wave function, creating up-electron on site 0
    for(int i=0; i<NHILBREDUCED; ++i)
    {
        //int ielec = int(i/NHILBPHON); // whether we are in sector 0, 1
        for(int j=0; j<NHILB; ++j)
        {
            Creator(i,j) = 0;
            //int jelec = int(j/NHILBPHON); // whether we are in sector 0, 1, 2, 3
            if(i==j-2*NHILBPHON)
                Creator(i,j) = 1; // really that simple ... note that this is not a diagonal matrix // we are mapping 0 => null, 1 => null, 2 => 0, 3 => 1
        }
    }
    */
    // END CONSTRUCTION FOR REDUCED HILBERT SPACE
    
    
    
    // diagonalize Ham0 to find ground state
    cmat Dum(NHILB, NHILB);
    vector<double> evals(NHILB);
    for(int i=0; i<NHILB; ++i)
        for(int j=0; j<NHILB; ++j)
            Dum(i,j) = Ham0(i,j);
    
    diagonalize(Dum, evals, NHILB);
    
    double groundstate_energy = evals[0];
    
    // print ground state energy
    if(myrank==0)
        cout << "lowest energies = " << evals[0] << ", " << evals[1] << endl;
    
    // check if ground state is correct
    vector<cdouble> groundstate(NHILB);
    for(int i=0; i<NHILB; ++i)
    {
        groundstate[i] = Dum(i,0);
    }
    
    vector<cdouble> Hgroundstate(NHILB);
    for(int i=0; i<NHILB; ++i)
        for(int j=0; j<NHILB; ++j)
            Hgroundstate[i] += Ham0(i,j) * groundstate[j];
    
    double energy = 0;
    for(int i=0; i<NHILB; ++i)
        energy += (conj(groundstate[i]) * Hgroundstate[i]).real();
    
    if(myrank==0)
        cout << "ground state energy checked = " << energy << endl;
    
    //double occ.
    vector<cdouble> Hdoccpsi0init(NHILB);
    for(int i=0; i<NHILB; ++i)
        for(int j=0; j<NHILB; ++j)
            Hdoccpsi0init[i] += Hamdocc(i,j) * groundstate[j];
    
    double doccinit = 0;
    for(int i=0; i<NHILB; ++i)
        doccinit += (conj(groundstate[i]) * Hdoccpsi0init[i]).real();
   

	//double occ. PHS
    vector<cdouble> Hdoccpsi0init_PHS(NHILB);
    for(int i=0; i<NHILB; ++i)
        for(int j=0; j<NHILB; ++j)
            Hdoccpsi0init_PHS[i] += Hamdocc_PHS(i,j) * groundstate[j];

    double doccinit_PHS = 0;
    for(int i=0; i<NHILB; ++i)
        doccinit_PHS += (conj(groundstate[i]) * Hdoccpsi0init_PHS[i]).real(); 
    
    // photon occ
    vector<cdouble> Hnphot(NHILB);
    for(int i=0; i<NHILB; ++i)
        for(int j=0; j<NHILB; ++j)
            Hnphot[i] += Hamnphot(i,j) * groundstate[j];
    
    double nphot = 0;
    for(int i=0; i<NHILB; ++i)
        nphot += (conj(groundstate[i]) * Hnphot[i]).real();
 /*   
    //calculation of photon distribution
    vector<double> photdistrib(NPHON);
    for (int i=0; i<NHILB;i++)
    {
        photdistrib[int(Hamnphot(i,i).real())]+=((conj(groundstate[i]))*groundstate[i]).real();
    }

      //calculation of coherent state distribution
    vector<double> photdistrib_coherent(NPHON);
    for (int j=0; j<NPHON;j++)
    {
        if (j==0)
        {
            photdistrib_coherent[j]=exp(-nphot);
        }
        else
        {
            photdistrib_coherent[j]=photdistrib_coherent[j-1]*(nphot/j);
        }
    }

  
  //Calculation of squeezed state distribution
    double rphot = log(sqrt(nphot)+sqrt(nphot+1));
    vector<long double> photdistrib_squeezed(NPHON);
    for (int j=0; j<NPHON; j+=2)
    {
        if (j==0)
        {
            photdistrib_squeezed[j]= 1./sqrt(cosh(rphot));
        }
        else
        {
            photdistrib_squeezed[j]=(photdistrib_squeezed[j-2]*tanh(rphot)*sqrt(j*(j-1)))/(2*(j/2));
        }
    }

     for (int i=0; i<NPHON; i++)
    {
         if (myrank==0) csvPhot_distrib << NPHON << ';' << G2 << ';' << WPLASMA << ';' << photdistrib[i] << ';'<< photdistrib_coherent[i] <<';' << i  << ';' << photdistrib_squeezed[i] << ';' << i << endl;
    }
*/
    // phonon occ
    vector<cdouble> Hnphon(NHILB);
    for(int i=0; i<NHILB; ++i)
        for(int j=0; j<NHILB; ++j)
            Hnphon[i] += Hamnphon(i,j) * groundstate[j];
    
    double nphon = 0;
    for(int i=0; i<NHILB; ++i)
        nphon += (conj(groundstate[i]) * Hnphon[i]).real();
   /* 
    //calculation of phonon distribution
    vector<double> phondistrib(NPHON);

    for (int i=0; i<NHILB;i++)
    {
        phondistrib[int(Hamnphon(i,i).real())]+=((conj(groundstate[i]))*groundstate[i]).real();
    }

    //calculation of coherent state distribution
    vector<double> phondistrib_coherent(NPHON);
    for (int j=0; j<NPHON;j++)
    {
        if (j==0)
        {
            phondistrib_coherent[j]=exp(-nphon);
        }
        else
        {
            phondistrib_coherent[j]=phondistrib_coherent[j-1]*(nphon/j);
        }
    }

     //Calculation of squeezed state distribution
    double rphon = log(sqrt(nphon)+sqrt(nphon+1));
    vector<long double> phondistrib_squeezed(NPHON);
    for (int j=0; j<NPHON; j+=2)
    {
        if (j==0)
        {
            phondistrib_squeezed[j]= 1./sqrt(cosh(rphon));
        }
        else
        {
            phondistrib_squeezed[j]=(phondistrib_squeezed[j-2]*tanh(rphon)*sqrt(j*(j-1)))/(2*(j/2));
        }
    }


    for (int i=0; i<NPHON; i++)
    {
         if (myrank==0) csvPhon_distrib << NPHON << ';' << G2 << ';' << WPLASMA << ';' << phondistrib[i] << ';'<< phondistrib_coherent[i] <<';' << i  << ';' << phondistrib_squeezed[i] << ';' << i << endl;
    }

         //calculation of photon/phonon distribution TWO MODES
    vector<double> TMdistrib(NPHON);

    for (int i=0; i<NHILB;i++)
    {
        if (int(Hamnphot(i,i).real())==int(Hamnphon(i,i).real()))
        {
                TMdistrib[int(Hamnphon(i,i).real())]+=((conj(groundstate[i]))*groundstate[i]).real();
        }
    }

         //Calculation of TWO MODE squeezed state distribution
    double n_average=(nphot+nphon)/2.;
    double r_TM = log(sqrt(n_average)+sqrt(n_average+1));
    vector<long double> TMdistrib_squeezed(NPHON);
    for (int j=0; j<NPHON; j++)
    {
        if (j==0)
        {
            TMdistrib_squeezed[j]= 1./(cosh(r_TM));
        }
        else
        {
            TMdistrib_squeezed[j]=TMdistrib_squeezed[j-1]*tanh(r_TM);
        }
    }

        for (int i=0; i<NPHON; i++)
    {
        if (myrank==0) csvPhotPhon_distrib << NPHON << ';' << G2 << ';' << WPLASMA << ';' << TMdistrib[i] << ';'<< TMdistrib_squeezed[i] <<';' << i  << ';'<< endl;
    }

*/
    // phonon <X^2>
    vector<cdouble> Hxsquared(NHILB);
    for(int i=0; i<NHILB; ++i)
        for(int j=0; j<NHILB; ++j)
            Hxsquared[i] += Hamxsquared(i,j) * groundstate[j];
    
    double xsquared = 0;
    for(int i=0; i<NHILB; ++i)
        xsquared += (conj(groundstate[i]) * Hxsquared[i]).real();
    
    // print here first groundstate properties
    
    if(myrank==0)
    {
        cout << "\n\n" << endl;
        cout << "/************ PARAMS ************/" << endl;
        cout << "Number phon/phot:         " << NPHON <<endl;
        cout << "Hopping J:                " << JHOP <<endl;
        cout << "Hubbard U:                " << U1 <<endl;
        cout << "Photon Omega:             " << OMEGAPHOT <<endl;
        cout << "Phonon Omega:             " << OMEGA <<endl;
        cout << "Phonon Quartic term:      " << QUARTIC <<endl;
        cout << "Photon-phonon coupling:   " << WPLASMA <<endl;
        cout << "Phonon-electron coupling: " << G2 << endl;
        cout << "/********************************/" << endl;
        cout << "Double occupancy:         " << doccinit << endl;
        cout << "Ueff:                     " << U1+G2*xsquared <<endl;
        cout << "OMEGAeff:                 " << OMEGA+(2*G2*doccinit_PHS) <<endl;
        cout << "nphot:                    " << nphot <<endl;
        cout << "nphon:                    " << nphon <<endl;
        cout << "xsquared:                 " << xsquared <<endl;
        cout << "/********************************/" << endl;

        csvGS_carac << WPLASMA << ';' << G2 << ';' << evals[0] << ';' << doccinit << ';' << U1+G2*xsquared << ';' << nphot << ';' << nphon << ';' << xsquared << ';' << OMEGA+(2*G2*doccinit_PHS) << endl;
    }
/*
    //// Construction of a photonic coherent state 

    vector<cdouble> coherent(NHILB);
    double nphot_coherent = 0.5; // New <n> for the coherent state,<n> = |alpha|^2
    vector<double> photdistrib_coherent_newstate(NPHON);

    //Calculate the new coherent distribution
    for (int j=0; j<NPHON;j++) 
    {
        if (j==0)
        {
            photdistrib_coherent_newstate[j]=exp(-nphot_coherent);
        }
        else
        {
            photdistrib_coherent_newstate[j]=photdistrib_coherent_newstate[j-1]*(nphot_coherent/j);
        }
    }

    // Construct the new coherent state with the new coherent distribution. Respect the ratio of the groundstate, not to be too different from initial GS

    for (int i=0; i<NHILB;i++)
    {
        int photnumber=int(Hamnphot(i,i).real());
        double oldproba_phot=photdistrib[photnumber];
        double newproba_phot=photdistrib_coherent_newstate[photnumber];
        coherent[i]=groundstate[i]*sqrt(newproba_phot/oldproba_phot);
    }

    // Not necessary. New photonic distribution of the coherent state
    // photdistrib_newstate should be = photdistrib_coherent_newstate

    vector<double> photdistrib_newstate(NPHON);
    for (int i=0; i<NHILB;i++)
    {
        photdistrib_newstate[int(Hamnphot(i,i).real())]+=((conj(coherent[i]))*coherent[i]).real();
    }
*/
    /////////// Initialize time stepping
  /*  
    double dtstep = TMAX/NUMT;
    double eta = 1./(2.*SIGARPES); // width for broadening
    double dom = (OMEGAMAX-OMEGAMIN)/(NOMEGAARPES);
    
    vector<cdouble> psi0(NHILB);
    for(int i=0; i<NHILB; ++i)
        psi0[i] = groundstate[i];
    
    vector<cdouble> psi1(NHILB);
    cmat A1(NHILB,NHILB);
    cmat A2(NHILB,NHILB);
    
    vector<double> tlist(NUMT+1);
    for(int it=0; it<=NUMT; ++it)
        tlist[it] = it*dtstep;
    
    // FOR ARPES PREPS
    // need to store wavefunction at ARPES time points
    vector<double> arpestlist(2*NTARPES+1);
    vector<int> tmapping(2*NTARPES+1); // maps to actual propagation time points; careful, this requires overlapping grids
    double dtarpes = 6.*SIGARPES/(2*NTARPES);
    for(int it=-NTARPES; it<=NTARPES; ++it)
    {
        arpestlist[it+NTARPES] = double(it)*dtarpes + TARPES; // TARPES is the central time
        tmapping[it+NTARPES] = int(arpestlist[it+NTARPES]/dtstep);
        if(myrank==0) cout << "ARPES step: " << it+NTARPES << " , actual time step: " << tmapping[it+NTARPES] << endl;
    }
    
    // wavefunction storage // only needed in reduced space!!!
    cmat psireduced(2*NTARPES+1, NHILBREDUCED);
    // this is for the addition of a particle in order to compute the addition spectrum
    cmat psiadded(2*NTARPES+1, NHILBREDUCED);
    // this is for spin-spin, full Hilbert space:
    cmat psispin(2*NTARPES+1, NHILB);
    // U time evolution operator
    cmat U(NHILB,NHILB);
    // END FOR ARPES PREPS
    
    if (FMAX==0)
    {
        for(int i=0; i<NHILB; ++i)
            for(int j=0; j<NHILB; ++j)
            {
                A1(i,j) = Ham0(i,j);
                A2(i,j) = Ham0(i,j);
            }

        
        set_U_step(U, A1, A2, evals, NHILB, dtstep);
    }
    
    // time stepping
    for(int it=1; it<=NUMT; ++it)
    {
        
        double time1 = tlist[it-1] + (0.5-sqrt(3.)/6.) * dtstep;
        double time2 = tlist[it-1] + (0.5+sqrt(3.)/6.) * dtstep;
        

        double field1 = FMAX * sin(OMEGAPUMP * time1);
        double field2 = FMAX * sin(OMEGAPUMP * time2);
        double field = FMAX * sin(OMEGAPUMP * tlist[it]);

        if (FMAX !=0)
        {
            for(int i=0; i<NHILB; ++i)
                for(int j=0; j<NHILB; ++j)
                {
                    A1(i,j) = Ham0(i,j) + field1 * Hamx(i,j);
                    A2(i,j) = Ham0(i,j) + field2 * Hamx(i,j);
                }

        
            set_U_step(U, A1, A2, evals, NHILB, dtstep);
        }
    */    
	/*for(int i=0; i<NHILB; ++i)
            for(int j=0; j<NHILB; ++j)
                psi1[i] += U(i,j) * psi0[j];*/
      /*  
 	cblas_zgemv(CblasColMajor, CblasNoTrans,
                        NHILB, NHILB,
                        &ONE, U.ptr(), NHILB,
                        psi0.data(), 1, &ZERO, psi1.data(), 1);
 
        for(int i=0; i<NHILB; ++i)
        {
            psi0[i] = psi1[i];
            psi1[i] = 0.;
        }
        
        
        // store ARPES wavefunction; this is very bad programming since I go thru the whole loop every time
        for(int ita=-NTARPES; ita<=NTARPES; ++ita)
        {
            if(tmapping[ita+NTARPES] == it) // times matching at end of propagation step
            {
                for(int i=0; i<NHILBREDUCED; ++i)
                {
                    psireduced(ita+NTARPES, i) = 0;
                    psiadded(ita+NTARPES, i) = 0;
                    for(int j=0; j<NHILB; ++j)
                    {
                        psireduced(ita+NTARPES, i) += Annihilator(i,j) * psi0[j];
                        psiadded(ita+NTARPES, i) += Creator(i,j) * psi0[j];
                    }
                }

                for(int i=0; i<NHILB; ++i)
                {
                    psispin(ita+NTARPES, i) = 0;
                    for(int j=0; j<NHILB; ++j)
                        psispin(ita+NTARPES, i) += Sz0(i,j) * psi0[j];
                }
            }
        }
         
        
        vector<cdouble> Hkinpsi0(NHILB);
        for(int i=0; i<NHILB; ++i)
            for(int j=0; j<NHILB; ++j)
                Hkinpsi0[i] += Hamkin(i,j) * psi0[j];
        
        double ekin = 0;
        for(int i=0; i<NHILB; ++i)
            ekin += (conj(psi0[i]) * Hkinpsi0[i]).real();
        
        vector<cdouble> Hpotpsi0(NHILB);
        for(int i=0; i<NHILB; ++i)
            for(int j=0; j<NHILB; ++j)
                Hpotpsi0[i] += Hamu(i,j) * psi0[j];
        
        double epot = 0;
        for(int i=0; i<NHILB; ++i)
            epot += (conj(psi0[i]) * Hpotpsi0[i]).real();
        
        //double occ.
        vector<cdouble> Hdoccpsi0(NHILB);
        for(int i=0; i<NHILB; ++i)
            for(int j=0; j<NHILB; ++j)
                Hdoccpsi0[i] += Hamdocc(i,j) * psi0[j];
        
        double docc = 0;
        for(int i=0; i<NHILB; ++i)
            docc += (conj(psi0[i]) * Hdoccpsi0[i]).real();
        
        vector<cdouble> Hphonpsi0(NHILB);
        for(int i=0; i<NHILB; ++i)
            for(int j=0; j<NHILB; ++j)
                Hphonpsi0[i] += Hamph(i,j) * psi0[j];
        
        double ephon = 0;
        for(int i=0; i<NHILB; ++i)
            ephon += (conj(psi0[i]) * Hphonpsi0[i]).real();
        
        vector<cdouble> Htotpsi0(NHILB);
        for(int i=0; i<NHILB; ++i)
            for(int j=0; j<NHILB; ++j)
                Htotpsi0[i] += Ham0(i,j) * psi0[j];
        
        double etot = 0;
        for(int i=0; i<NHILB; ++i)
            etot += (conj(psi0[i]) * Htotpsi0[i]).real();
        
        vector<cdouble> Hxpsi0(NHILB);
        for(int i=0; i<NHILB; ++i)
            for(int j=0; j<NHILB; ++j)
                Hxpsi0[i] += Hamx(i,j) * psi0[j];
        
        double xt = 0;
        for(int i=0; i<NHILB; ++i)
            xt += (conj(psi0[i]) * Hxpsi0[i]).real();
        
        if(myrank==0)
        {
            cout << "time = " << tlist[it];
            cout << ", kinetic energy = " << ekin;
            cout << ", Hubbard energy = " << epot;
            cout << ", electronic energy = " << ekin+epot;
            cout << ", phononic energy = " << ephon;
            cout << ", x(t) = " << xt;
            cout << ", total energy = " << etot;
            cout << ", double occ = " << docc << endl;
            
 //           csvtime << tlist[it] << ';' << ekin << ';' << epot << ';' << ekin+epot << ';' << ephon << ';' << etot << ';' << xt << ';' << field << ';' << docc << endl;
        }
        
    }
    
    
    // construction for spectrum: need to store wave functions at some predetermined number of time steps grid (reduced sampling?!) - maybe for 50 times or so, depending on probe width
    // then apply Annihilator and propagate using reduced Hilbert space Hamiltonian, accumulate directly to spectrum
    // fix the center of the probe window!!!
    // propagate only for half the times
    
    // now propagate the reduced wavefunction in reduced Hilbert space
    
    
    vector<cdouble> psi0reduced(NHILBREDUCED);
    vector<cdouble> psi1reduced(NHILBREDUCED);
    cmat A1reduced(NHILBREDUCED,NHILBREDUCED);
    cmat A2reduced(NHILBREDUCED,NHILBREDUCED);
    vector<double> evalsreduced(NHILBREDUCED);
    
    vector<cdouble> psi0added(NHILBREDUCED);
    vector<cdouble> psi1added(NHILBREDUCED);
    cmat A1added(NHILBREDUCED,NHILBREDUCED);
    cmat A2added(NHILBREDUCED,NHILBREDUCED);
    vector<double> evalsadded(NHILBREDUCED);

    // need to also propagate spin-measured wavefunction
    vector<cdouble> psi0spin(NHILB);
    vector<cdouble> psi1spin(NHILB);
    cmat A1spin(NHILB,NHILB);
    cmat A2spin(NHILB,NHILB);
    vector<double> evalsspin(NHILB);
    
    double shaper = 0.;
    vector<double> arpessignal(NOMEGAARPES+1);
    vector<double> spinsignal(NOMEGAARPES+1);
    vector<double> retardedsignal(NOMEGAARPES+1);
    double domega = (OMEGAMAX-OMEGAMIN)/(NOMEGAARPES);

    cmat Ureduced(NHILBREDUCED, NHILBREDUCED);
    cmat Uadded(NHILBREDUCED, NHILBREDUCED);
    cmat Uspin(NHILB, NHILB);
    
    if(myrank==0) cout << "Starting ARPES postprocessing" << endl;

    if (FMAX ==0)
            {
            
                for(int i=0; i<NHILBREDUCED; ++i)
                    for(int j=0; j<NHILBREDUCED; ++j)
                    {
                        A1reduced(i,j) = Ham0red(i,j);
                        A2reduced(i,j) = Ham0red(i,j);
                    
                        A1added(i,j) = Ham0added(i,j);
                        A2added(i,j) = Ham0added(i,j);
                    }

                for(int i=0; i<NHILB; ++i)
                    for(int j=0; j<NHILB; ++j)
                    {
                        A1spin(i,j) = Ham0(i,j); 
                        A2spin(i,j) = Ham0(i,j);
                    }
            
            
                set_U_step(Ureduced, A1reduced, A2reduced, evalsreduced, NHILBREDUCED, dtstep);
                set_U_step(Uadded, A1added, A2added, evalsadded, NHILBREDUCED, dtstep);      
                set_U_step(Uspin, A1spin, A2spin, evalsspin, NHILB, dtstep);

            }
    
    for(int it=-NTARPES+myrank; it<=NTARPES; it+=numprocs)
    //for(int it=-NTARPES; it<=NTARPES; it+=1)
    {
        // start propagation from respective reference state at (t,t)
        // only propagate forward in time
        
      if(myrank==0) cout << "ARPES at time step " << it+NTARPES << " out of " << 2*NTARPES+1 << endl;
        
        for(int i=0; i<NHILBREDUCED; ++i)
        {
            psi0reduced[i] = psireduced(it+NTARPES, i);
            psi0added[i] = psiadded(it+NTARPES, i);
        }
        for(int i=0; i<NHILB; ++i)
        {
            psi0spin[i] = psispin(it+NTARPES, i);
        }
        

        for(double time=arpestlist[it+NTARPES]; time<arpestlist[2*NTARPES]; time += dtstep)
        {
            double time1 = time + (0.5-sqrt(3.)/6.) * dtstep;
            double time2 = time + (0.5+sqrt(3.)/6.) * dtstep;
            double field1 = FMAX * sin(OMEGAPUMP * time1);
            double field2 = FMAX * sin(OMEGAPUMP * time2);
            
            if (FMAX !=0)
            {
            
                for(int i=0; i<NHILBREDUCED; ++i)
                    for(int j=0; j<NHILBREDUCED; ++j)
                    {
                        A1reduced(i,j) = Ham0red(i,j) + field1 * Hamxred(i,j);
                        A2reduced(i,j) = Ham0red(i,j) + field2 * Hamxred(i,j);
                    
                        A1added(i,j) = Ham0added(i,j) + field1 * Hamxadded(i,j);
                        A2added(i,j) = Ham0added(i,j) + field2 * Hamxadded(i,j);
                    }

                for(int i=0; i<NHILB; ++i)
                    for(int j=0; j<NHILB; ++j)
                    {
                        A1spin(i,j) = Ham0(i,j) + field1 * Hamx(i,j); 
                        A2spin(i,j) = Ham0(i,j) + field2 * Hamx(i,j);
                    }
            
            
                set_U_step(Ureduced, A1reduced, A2reduced, evalsreduced, NHILBREDUCED, dtstep);
                set_U_step(Uadded, A1added, A2added, evalsadded, NHILBREDUCED, dtstep);      
                set_U_step(Uspin, A1spin, A2spin, evalsspin, NHILB, dtstep);

            }

            cblas_zgemv(CblasColMajor, CblasNoTrans,
                        NHILBREDUCED, NHILBREDUCED,
                        &ONE, Ureduced.ptr(), NHILBREDUCED,
                        psi0reduced.data(), 1, &ZERO, psi1reduced.data(), 1);
            cblas_zgemv(CblasColMajor, CblasNoTrans,
                        NHILBREDUCED, NHILBREDUCED,
                        &ONE, Uadded.ptr(), NHILBREDUCED,
                        psi0added.data(), 1, &ZERO, psi1added.data(), 1);
            cblas_zgemv(CblasColMajor, CblasNoTrans,
                        NHILB, NHILB,
                        &ONE, Uspin.ptr(), NHILB,
                        psi0spin.data(), 1, &ZERO, psi1spin.data(), 1);
            
            for(int i=0; i<NHILBREDUCED; ++i)
            {
                psi0reduced[i] = psi1reduced[i];
                psi1reduced[i] = 0.;
                
                psi0added[i] = psi1added[i];
                psi1added[i] = 0.;
            }

             for(int i=0; i<NHILB; ++i)
            {
                psi0spin[i] = psi1spin[i];
                psi1spin[i] = 0.;
            }
            
            
            
            
            // accumulate to ARPES signal at t'-t when congruent with arpestlist
            // use symmetry for other half of the t'-t triangle
            for(int ittest=-NTARPES; ittest<=NTARPES; ++ittest)
            {
                if(arpestlist[ittest+NTARPES] == time+dtstep) // the time to which we have just propagated
                {
                    // compute overlap
                    cdouble overlap = 0.;
                    cdouble overlapadded = 0.;
                    cdouble overlapspin = 0.;
                    for(int i=0; i<NHILBREDUCED; ++i)
                    {
                        overlap += conj(psireduced(ittest+NTARPES, i)) * psi0reduced[i];
                        overlapadded += conj(psiadded(ittest+NTARPES, i)) * psi0added[i];   
                    }
                    
                    for(int i=0; i<NHILB; ++i)
                    {
                        overlapspin += conj(psispin(ittest+NTARPES, i)) * psi0spin[i];
                    }

                    // add to ARPES
                    
                    for(int iw = 0; iw <= NOMEGAARPES; ++iw)
                    {
                        double omegaarpes = OMEGAMIN + domega * iw;
                        shaper = arpesshape(arpestlist[ittest+NTARPES], arpestlist[it+NTARPES], TARPES, SIGARPES);

                        arpessignal[iw] += 2.*real(exp(II * omegaarpes * (arpestlist[it+NTARPES] - arpestlist[ittest+NTARPES])) * overlap) *  dtarpes * dtarpes * shaper;
                        
                        retardedsignal[iw] += 2.*real(exp(II * omegaarpes * (arpestlist[ittest+NTARPES] - arpestlist[it+NTARPES])) * overlapadded) *  dtarpes * dtarpes * shaper;
                        
                        retardedsignal[iw] += 2.*real(exp(II * omegaarpes * (arpestlist[it+NTARPES] - arpestlist[ittest+NTARPES])) * overlap) *  dtarpes * dtarpes * shaper;
                        
                        spinsignal[iw] += 1.*real(exp(II * omegaarpes * (arpestlist[it+NTARPES] - arpestlist[ittest+NTARPES])) * overlapspin) *  dtarpes * dtarpes * shaper;
                        
                        spinsignal[iw] += 1.*real(exp(-II * omegaarpes * (arpestlist[it+NTARPES] - arpestlist[ittest+NTARPES])) * overlapspin) *  dtarpes * dtarpes * shaper;
                    }
                }
            }
            
            
            
        }
            
        cdouble overlapdiag = 0.; // not necessarily = 1 in general because it contains info about particle number
        
        cdouble overlapdiagadded = 0.;

        cdouble overlapdiagspin = 0.;
        
        for(int i=0; i<NHILBREDUCED; ++i)
        {
            overlapdiag += conj(psireduced(it+NTARPES, i)) * psireduced(it+NTARPES, i);
            overlapdiagadded += conj(psiadded(it+NTARPES, i)) * psiadded(it+NTARPES, i);
        }

        for(int i=0; i<NHILB; ++i)
            overlapdiagspin += conj(psispin(it+NTARPES, i)) * psispin(it+NTARPES, i);

        
        // add diagonal t-t term at the end
        shaper = arpesshape(arpestlist[it+NTARPES], arpestlist[it+NTARPES], TARPES, SIGARPES);
    
        for(int iw = 0; iw <= NOMEGAARPES; ++iw)
        {
            arpessignal[iw] += real(overlapdiag) * dtarpes * dtarpes * shaper ;
            retardedsignal[iw] += real(overlapdiag) * dtarpes * dtarpes * shaper ;
            retardedsignal[iw] += real(overlapdiagadded) * dtarpes * dtarpes * shaper ;
            spinsignal[iw] += real(overlapdiagspin) * dtarpes * dtarpes * shaper ;
        }
        
    }
    
#ifndef NO_MPI
    MPI_Allreduce(MPI_IN_PLACE, &arpessignal[0], NOMEGAARPES+1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &retardedsignal[0], NOMEGAARPES+1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    
    for(int iw = 0; iw <= NOMEGAARPES; ++iw)
    {
        //if(myrank==0) outarpes << OMEGAMIN + domega * iw << '\t' << arpessignal[iw] << '\t' << retardedsignal[iw] << endl;
        if(myrank==0) csvarpes2 << OMEGAMIN + domega * iw << ';' << arpessignal[iw] << ';' << retardedsignal[iw] << endl;
        if(myrank==0) cout << OMEGAMIN + domega * iw << ' ' << arpessignal[iw] << ' ' << retardedsignal[iw] << endl;
    }
*/
  /*   for(int iw = 0; iw <= NOMEGAARPES; ++iw)
    {
        if(myrank==0) csvspin2 << OMEGAMIN + domega * iw << ';' << spinsignal[iw]  << endl;
        if(myrank==0) cout << OMEGAMIN + domega * iw << ' ' << spinsignal[iw] <<  endl;
    }*/

    
     /*Spec weight*/
  /*  double sum=0;
    double spec_weight=0;
    for (int iw=1;iw<NOMEGAARPES; ++iw)
        sum+=retardedsignal[iw];

    spec_weight=(domega/2)*(retardedsignal[0]+retardedsignal[NOMEGAARPES]+(2*sum));
    if(myrank==0) cout << " Spec weight" << spec_weight << endl;
   // if(myrank==0) outspecweight << spec_weight << endl;
*/
    if(myrank==0)
        cout << "\n ... done! \n"<<flush;
    
#ifndef NO_MPI
	MPI_Finalize();
#endif
    
    return 0;

}

int main(int argc, char* argv[])
{
    return ed_ssh_main(argc, argv);
}

