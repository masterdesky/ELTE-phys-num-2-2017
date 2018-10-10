#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M_F 0.59736e25
#define M_H 7.349e22
#define xa 4.055e8
#define va 964
#define xp 3.633e8
#define vp 1076
#define G 0.667384e-10


double derivalt(double* x, int i)
{
double dx[4];
 if(i<2)
 {

     dx[i]=x[i+2];

 }
 if(i>=2)
 {

      double absz=pow((x[0]*x[0]+x[1]*x[1]),1.5);
      dx[i]=-G*M_F*x[i-2]/(absz);
 }


 return dx[i];
}


void energy (double vx, double vy, double x, double y, FILE* g, double t)
{
    double* k;
    double E;
    double absz=pow(x*x+y*y,0.5);
    k=&E;
    E=0.5*(vx*vx+vy*vy)*M_H-(G*M_H*M_F)/(absz);
    fprintf(g,"%f %f\n",t,*k);
}



void eulerstep(double t, double* yr, double* yn, double dt, double time , int N , FILE* f, FILE* g)
{
    int i;
    double* z;

    t=0;
    while(t < time)
    {

       for(i=0; i<N; ++i)
        {
        yn[i]=yr[i]+dt*derivalt( yr, i);
        }
    fprintf(f, "%f %f\n", yn[0] , yn[1]);
    t+=dt;
    energy(yn[2],yn[3], yn[0], yn[1], g,t);
    z = yr;
    yr = yn;
    yn = z;
    }
}





int main(int argc, char *argv[])
{



double r[4];
double rn[4];



r[0]=xa;
r[1]=0;
r[2]=0;
r[3]=va;
rn[0]=0;
rn[1]=0;
rn[2]=0;
rn[3]=0;

FILE* m=fopen("coord.txt","w");
FILE* e=fopen("energ.txt","w");


double time=30*24*60*60*10;
double dt=1000;
double t=0;


eulerstep(t,r,rn,dt,time,4,m,e);

fclose(m);
fclose (e);
return 0;
}
