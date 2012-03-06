#include <stdio.h>
#include <math.h>
#include "complex.h"

#define MAXLENGTH 30

double hbar = 6.626e-34 / 2 / 3.141592653589793;
double m0 = 9.109e-31;
double e0 = 1.602e-19;

void psiFn(double Eq, int startpoint, int xPsiSize, double xres, double *xVc, double *xEg, 
           double *xF, double *xEp, double *xESO, double *xMc, double *xMcE, 
           double *xPsi)
{
  int q=0;
  for(q=0;q<xPsiSize;q++)
  { if(1)
        *(xMcE+q) = m0 / ((1+2* *(xF+q)) + *(xEp+q)/3 * (2 / ((Eq-*(xVc+q))+*(xEg+q)) 
                       + 1 / ((Eq-*(xVc+q))+*(xEg+q)+*(xESO+q)) ));
    else
        *(xMcE+q) = m0 * *(xMc+q) * (1 - (*(xVc+q) - Eq) / *(xEg+q));
    if(q>1)
      *(xMcE+q-1) = 0.5 * (*(xMcE+q) + *(xMcE+q-1));
  }
  for(q=0;q<startpoint;q++)
    *(xPsi+q) = 0;
  *(xPsi+startpoint) = 1;
  //*xPsi=0;
  //*(xPsi+1)=1;
  for(q=startpoint;q<xPsiSize-1;q++)
  {
    *(xPsi+q+1) = (
                   (2*xres*1e-10*xres*1e-10 /hbar/hbar * (*(xVc+q) - Eq)*e0 
                     + 1 / *(xMcE+q) + 1 / *(xMcE+q-1)) * *(xPsi+q) 
                   - *(xPsi+q-1) / *(xMcE+q-1)
                  ) * *(xMcE+q);
  }
  return;
}

int psiFnEnd(double *eEq, int eEqSize, int xPsiSize, double xres, double EField,
             double *xVc, double *xEg, double *xF, double *xEp, double *xESO, 
             double *xMc, double *xMcE, double *xPsi, double *xPsiEnd)
{
  double Eq=0;
  int startpoint=1;
  int q=0;
  for(q=0;q<eEqSize;q++)
  {
    Eq = *(eEq+q);
    startpoint = xPsiSize - ceil((Eq - *eEq)/(EField * xres)*1e5 + 200/xres);
    if(startpoint<1)
      startpoint = 1;

    psiFn(Eq, startpoint, xPsiSize, xres, xVc, xEg, xF, xEp, xESO, xMc, xMcE, xPsi);
    *(xPsiEnd+q) = *(xPsi+xPsiSize-1);
    //printf("%d: %g %d        ", q, *(eEq+q), startpoint);
    //printf("%d  ", startpoint);
  }

  return 1;
}

int inv_quadratic_interp(double *xnew, double *ynew, double *idxs, int EigenESize, double *EigenE)
{
  double x0=0, fx0=0,
         x1=0, fx1=0,
         x2=0, fx2=0, x3=0;
  int idx=0;
  int q=0;
  for(q=0;q<EigenESize;q++)
  {
    idx = *(idxs+q);
    printf("%d\n",idx);
    x0=*(xnew+idx-1); fx0=*(ynew+idx-1);
    x1=*(xnew+idx);   fx1=*(ynew+idx);
    x2=*(xnew+idx+1); fx2=*(ynew+idx+1);
    x3 = x0*fx1*fx2/(fx0-fx1)/(fx0-fx2) + x1*fx0*fx2/(fx1-fx0)/(fx1-fx2) + x2*fx0*fx1/(fx2-fx0)/(fx2-fx1);
    *(EigenE+q) = x3;
  }
  return 1;
}


int returnme()
{return 42;}

int psiFill(int xPsiSize, double xres, int EigenESize, double *EigenE, double *xVc, double *xEg, 
            double *xF, double *xEp, double *xESO, double *xMc, double *xMcE, double *xyPsi)
{
  int col=0; //column
  double Eq=0;
  for(col=0;col<EigenESize;col++)
  {

     Eq = *(EigenE+col);
     int q=0;
      for(q=0;q<xPsiSize;q++)
      { if(1)
            *(xMcE+q) = m0 / ((1+2* *(xF+q)) + *(xEp+q)/3 * (2 / ((Eq-*(xVc+q))+*(xEg+q)) 
                           + 1 / ((Eq-*(xVc+q))+*(xEg+q)+*(xESO+q)) ));
        else
            *(xMcE+q) = m0 * *(xMc+q) * (1 - (*(xVc+q) - Eq) / *(xEg+q));
        if(q>1)
          *(xMcE+q-1) = 0.5 * (*(xMcE+q) + *(xMcE+q-1));
      }
      *(xyPsi+0+col*xPsiSize) = 0.;
      *(xyPsi+1+col*xPsiSize) = 1.;
      double PsiInt = 1; //one for xPsi[1]
      for(q=1;q<xPsiSize-1;q++)
      {
        *(xyPsi+q+1+col*xPsiSize) = (
                       (2*xres*1e-10*xres*1e-10 /hbar/hbar * (*(xVc+q) - Eq)*e0 
                         + 1 / *(xMcE+q) + 1 / *(xMcE+q-1)) * *(xyPsi+q+col*xPsiSize)
                       - *(xyPsi+q-1+col*xPsiSize) / *(xMcE+q-1)
                      ) * *(xMcE+q);
        PsiInt += *(xyPsi+q+col*xPsiSize) * *(xyPsi+q+col*xPsiSize) * (1+(Eq-*(xVc+q))/(Eq-*(xVc+q)+*(xEg+q)));
      }
      double A = 1 / sqrt(xres * 1e-10 * PsiInt);
                
      for(q=0;q<xPsiSize;q++)
        *(xyPsi+q+col*xPsiSize) *= A;



  }
  return 1;
}


typedef struct {complex aa; complex ab; complex ba; complex bb;} matrix;

matrix identity()
{
    matrix c;
    c.aa = cmplx(1,0);
    c.ab = cmplx(0,0);
    c.ba = cmplx(0,0);
    c.bb = cmplx(1,0);
    return c;
}

matrix mmult(matrix m1, matrix m2)
{
    matrix c;
    c.aa = cxadd(cxmul(m1.aa,m2.aa) , cxmul(m1.ab,m2.ba));
    c.ab = cxadd(cxmul(m1.aa,m2.ab) , cxmul(m1.ab,m2.bb));
    c.ba = cxadd(cxmul(m1.ba,m2.aa) , cxmul(m1.bb,m2.ba));
    c.bb = cxadd(cxmul(m1.ba,m2.ab) , cxmul(m1.bb,m2.bb));
    return c;
}

void chiImag_array(double wavelength, const double *thicknesses, const double *indexesReal, 
                const double *indexesImag, int numLayers, double *betaInReal, 
                double *betaInImag, int numBetas, double *chiImag)
{

    double pi = 3.1415926535897932385;
    double k = 2 * pi / wavelength;
    double z0 = 0.003768;

    int q=0;
    for (q=0; q<numBetas; q++)
    {
        int j = 0;
    
        complex beta = cmplx(*(betaInReal+q), *(betaInImag+q));

        complex index[MAXLENGTH];
        for (j=0; j < numLayers; j++)
            index[j] = cmplx(indexesReal[j],indexesImag[j]);
        
        //alpha = sqrt(self.stratumRIndexes**2-beta**2)
        complex alpha[MAXLENGTH];
        for (j=0; j < numLayers; j++)
            alpha[j] = cxsqrt(cxsub(cxsqr(index[j]),cxsqr(beta)));
        //make sure correct sign of alphac and alphas are chosen, see Chilwell
        if(imag(alpha[0]) < 0)
            alpha[0] = cxneg(alpha[0]);
        if(imag(alpha[numLayers-1]) < 0)
            alpha[numLayers-1] = cxneg(alpha[numLayers-1]);
    
        //gamma = z0*alpha/self.stratumRIndexes**2
        complex gamma[MAXLENGTH];
        for (j=0; j < numLayers; j++)
            gamma[j] = cmuld(cxdiv(alpha[j],cxsqr(index[j])),z0);
        
        //phi = k*self.stratumThicknesses*alpha
        complex phi[MAXLENGTH];
        for (j=0; j < numLayers; j++)
            phi[j] = cmuld(alpha[j],k*thicknesses[j]);
    
        double zeta[MAXLENGTH];
        for (j=0; j < numLayers; j++)
            zeta[j] = k*thicknesses[j]/z0;
        
        complex chi;
        matrix m=identity(), mt=identity();
        complex ni=cmplx(0,-1);
    
        for (j=numLayers-1; j>-1; j--)
        //array([[cos(phi[q]), -1j/gamma[q]*sin(phi[q])],[-1j*gamma[q]*sin(phi[q]), cos(phi[q])]])
        {
            if (real(index[j]) == real(beta) && imag(index[j]) == imag(beta))
                mt.ab = cxmul(cmuld(ni,zeta[j]),cxsqr(index[j]));  // -i*k*thickness*n^2/z0
            else
                mt.ab = cxdiv(cxmul(ni,cxsin(phi[j])),gamma[j]);  // -i*sin(phi)/gamma
    
            mt.aa = cxcos(phi[j]);
            mt.ba = cxmul(cxmul(ni,cxsin(phi[j])),gamma[j]);   // -i*sin(phi)*gamma
            mt.bb = mt.aa;
            
            m = mmult(mt,m);
        }
        chi = cxadd4(cxmul(gamma[numLayers-1],m.aa) , 
                           cxmul(gamma[0],cxmul(gamma[numLayers-1],m.ab)) , 
                           m.ba , 
                           cxmul(gamma[0],m.bb) );
        //chi = gammac*M[0,0] + gammac*gammas*M[0,1] + M[1,0] + gammas*M[1,1]
        
        *(chiImag+q) = imag(chi);
    }
}

double abschi_find(double wavelength, const double *thicknesses, const double *indexesReal, 
                const double *indexesImag, int numLayers, double betaInReal, double betaInImag)
{

    double pi = 3.1415926535897932385;
    double k = 2 * pi / wavelength;
    double z0 = 0.003768;
    int j = 0;

    complex beta = cmplx(betaInReal, betaInImag);
    complex index[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        index[j] = cmplx(indexesReal[j],indexesImag[j]);
    
    //alpha = sqrt(self.stratumRIndexes**2-beta**2)
    complex alpha[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        alpha[j] = cxsqrt(cxsub(cxsqr(index[j]),cxsqr(beta)));
    //make sure correct sign of alphac and alphas are chosen, see Chilwell
    if(imag(alpha[0]) < 0)
        alpha[0] = cxneg(alpha[0]);
    if(imag(alpha[numLayers-1]) < 0)
        alpha[numLayers-1] = cxneg(alpha[numLayers-1]);

    //gamma = z0*alpha/self.stratumRIndexes**2
    complex gamma[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        gamma[j] = cmuld(cxdiv(alpha[j],cxsqr(index[j])),z0);
    
    //phi = k*self.stratumThicknesses*alpha
    complex phi[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        phi[j] = cmuld(alpha[j],k*thicknesses[j]);

    double zeta[MAXLENGTH];
    for (j=0; j < numLayers; j++)
        zeta[j] = k*thicknesses[j]/z0;
    
    complex chi;
    matrix m=identity(), mt=identity();
    complex ni=cmplx(0,-1);

    for (j=numLayers-1; j>-1; j--)
    //array([[cos(phi[q]), -1j/gamma[q]*sin(phi[q])],[-1j*gamma[q]*sin(phi[q]), cos(phi[q])]])
    {
        if (real(index[j]) == real(beta) && imag(index[j]) == imag(beta))
            mt.ab = cxmul(cmuld(ni,zeta[j]),cxsqr(index[j]));  // -i*k*thickness*n^2/z0
        else
            mt.ab = cxdiv(cxmul(ni,cxsin(phi[j])),gamma[j]);  // -i*sin(phi)/gamma

        mt.aa = cxcos(phi[j]);
        mt.ba = cxmul(cxmul(ni,cxsin(phi[j])),gamma[j]);   // -i*sin(phi)*gamma
        mt.bb = mt.aa;
        
        m = mmult(mt,m);
    }
    chi = cxadd4(cxmul(gamma[numLayers-1],m.aa) , 
                       cxmul(gamma[0],cxmul(gamma[numLayers-1],m.ab)) , 
                       m.ba , 
                       cxmul(gamma[0],m.bb) );
    //chi = gammac*M[0,0] + gammac*gammas*M[0,1] + M[1,0] + gammas*M[1,1]
    return cxabs(chi);  // return a double value
}

int argmin(const double *values, int numElements)
{
    int idx=0;
    double minValue = values[0];
    int j=0;
    for (j=1; j<numElements; j++)
        idx = values[j] < minValue ? j : idx;
    return idx;
}

void beta_find(double wavelength, const double *thicknesses, const double *indexesReal, 
               const double *indexesImag, int numLayers, double betaInReal, double betaInImag,
               double beta_find_precision, double *betaOut)
{
    //initialize beta0
    complex beta0 = cmplx(betaInReal, betaInImag);

    complex rInc = cmplx(0.0001,0);
    complex iInc = cmplx(0,1.0e-6);

    int numBetas = 9;
    complex betas[9] = {0};

    double abschiOld = 0, abschiNew = 0;
    int chiMinIdx = 0;

    do{
        //set betas
        betas[0] = beta0;
        betas[1] = cxadd(beta0,rInc); 
        betas[2] = cxsub(beta0,rInc);
        betas[3] = cxadd(beta0,iInc); 
        betas[4] = cxsub(beta0,iInc);
        betas[5] = cxadd(cxadd(beta0,rInc),iInc);
        betas[6] = cxsub(cxsub(beta0,rInc),iInc);
        betas[7] = cxadd(cxsub(beta0,rInc),iInc);
        betas[8] = cxsub(cxadd(beta0,rInc),iInc);
    
        double abschi[numBetas];
        int j=0;
        for (j=0; j<numBetas; j++)
        {
            abschi[j] = abschi_find(wavelength, thicknesses, indexesReal, indexesImag, 
                                    numLayers, real(betas[j]), imag(betas[j]));
        }
        chiMinIdx = argmin(abschi, numBetas);
        abschiOld = abschiNew;
        abschiNew = abschi[chiMinIdx];
        beta0 = betas[chiMinIdx];
    }while(abs(abschiOld - abschiNew)/abschiOld > beta_find_precision);

    *betaOut = real(beta0);
    *(betaOut+1) = imag(beta0);
}

