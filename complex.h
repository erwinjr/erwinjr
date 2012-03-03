/* complex.h            define type complex and arithmetic operations */
#ifndef COMPLEX_H_
#define COMPLEX_H_


typedef struct {double re; double im;} complex;

complex cmplx(double  a, double  b); /* create a complex number            */  
double  real (complex a);            /* get real part                      */
double  imag (complex a);            /* get imaginary part                 */
complex cxadd(complex a, complex b); /* add two complex numbers            */
complex cxadd4(complex a, complex b, complex c, complex d);
complex cxsub(complex a, complex b); /* subtract b from a, complex         */
complex cxmul(complex a, complex b); /* multiply two complex numbers       */
complex cmuld(complex a, double  b); /* multiply complex times double      */
complex cxdiv(complex a, complex b); /* divide complex by complex          */
complex cdivd(complex a, double  b); /* divide complex by double           */
double  cxabs(complex a);            /* return magnitude of complex munber */
double  cxarg(complex a);            /* return argument of complex number  */
complex cxneg(complex a);            /* negate a complex number            */
complex cxsqr(complex a);            /* complex square                     */
complex cxsqrt(complex a);           /* complex square root                */
complex cxsin(complex a);            /* complex sin                        */
complex cxcos(complex a);            /* complex cos                        */


complex cmplx(double  a, double  b)  /* make a complex number */
{
  complex c;
  c.re = a;
  c.im = b;
  return c;
}

double real(complex a) /* get real part */
{
  return a.re;
}

double imag(complex a) /* get imaginary part */
{
  return a.im;
}

complex cxadd(complex a, complex b)  /* add two complex numbers */
{
  complex c;
  c.re = a.re+b.re;
  c.im = a.im+b.im;
  return c;
}

complex cxadd4(complex a, complex b, complex c, complex d)
{
  complex e;
  e.re = a.re+b.re+c.re+d.re;
  e.im = a.im+b.im+c.im+d.im;
  return e;
}

complex cxsub(complex a, complex b)  /* subtract b from a, complex */
{
  complex c;
  c.re = a.re-b.re;
  c.im = a.im-b.im;
  return c;
}

complex cxmul(complex a, complex b)  /* multiply two complex numbers */
{
  complex c;
  c.re = a.re*b.re - a.im*b.im;
  c.im = a.re*b.im + a.im*b.re;
  return c;
}

complex cmuld(complex a, double b)  /* multiply complex times double */
{
  complex c;
  c.re = a.re*b;
  c.im = a.im*b;
  return c;
}

complex cxdiv(complex a, complex b)  /* divide complex by complex */
{
  complex c;
  double  r;
  r = b.re*b.re+b.im*b.im;
  c.re = (a.re*b.re+a.im*b.im)/r;
  c.im = (a.im*b.re-a.re*b.im)/r;
  return c;
}

complex cdivd(complex a, double b)  /* divide complex by double */
{
  complex c;
  c.re = a.re/b;
  c.im = a.im/b;
  return c;
}

double cxabs(complex a)             /* return magnitude of complex munber */
{
  return sqrt(a.re*a.re+a.im*a.im);
}

double cxarg(complex a)             /* return argument of complex munber */
{
  return atan2(a.im,a.re);
}

complex cxneg(complex a)             /* negate a complex */
{
  complex c;
  c.re = -a.re;
  c.im = -a.im;
  return c;
}

complex cxsqr(complex a)
{
    complex c;
    c.re = a.re*a.re - a.im*a.im;
    c.im = 2 * a.re * a.im;
    return c;
}

complex cxsqrt(complex a)
{
    complex c;
    c.re = sqrt((cxabs(a)+real(a))/2);
    if(imag(a)<0)
        c.im = -sqrt((cxabs(a)-real(a))/2);
    else
        c.im = sqrt((cxabs(a)-real(a))/2);
    return c;
}
    
complex cxsin(complex a)
{
    complex c;
    c.re = sin(a.re)*cosh(a.im);
    c.im = cos(a.re)*sinh(a.im);
    return c;
}

complex cxcos(complex a)
{
    complex c;
    c.re = cos(a.re)*cosh(a.im);
    c.im = -sin(a.re)*sinh(a.im);
    return c;
}


#endif /* COMPLEX_H_ */