#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// i edik derivált kiszámítása
double dy(double* t, double* y, int i, int N)
{
    double G=0.6674e-10;
    double M= 0.59736e25;
    double derivalt;

    if (i < 2)
    {
        derivalt=y[i+2];
    }
    else
    {
        derivalt=-G*M*y[i-2]/pow(y[0]*y[0] + y[1]*y[1],1.5);
    }

   return derivalt;
}

void energy(double* t, double* yn, FILE* e )
{
 double v;
 double m =0.7349e23;
 double G=0.6674e-10;
 double M= 0.59736e25;
 double absz=pow(yn[0]*yn[0]+yn[1]*yn[1],0.5);

 v=0.5 * m * (yn[2]*yn[2]+yn[3]*yn[3]) - (G*M*m)/(absz);;
   fprintf(e, "%f %f\n", *t , v);

}

void eulerstep(double* t, double* yr, double* yn, double dt, double time , int N , FILE* f, FILE* e)
{
    int i;
    double* c;
    *t=0;
    while(*t < time)
    {
        for(i=0; i<N; ++i)
        {
            yn[i]=yr[i]+dt*dy(t, yr, i, N);
        }
        fprintf(f, "%f %f\n", yn[0] , yn[1]);
        *t+=dt;
        energy(t,yn,e);

        yr = yn;
    }
}

int main()
{

double x=0;         double* X;    X=&x;
double y=0;         double* Y;    Y=&y;
double vx=0;        double* VX;   VX=&vx;
double vy=0;        double* VY;   VY=&vy;
double x_r=405500000;  double* X_R;  X_R=&x_r;
double y_r=0;       double* Y_R;  Y_R=&y_r;
double vx_r=0;      double* VX_R; VX_R=&vx_r;
double vy_r=964;    double* VY_R; VY_R=&vy_r;
double t=0;         double* T;    T=&t;

double yr[4];
double yn[4];

yn[0]=*X;
yn[1]=*Y;
yn[2]=*VX;
yn[3]=*VY;

yr[0]=*X_R;
yr[1]=*Y_R;
yr[2]=*VX_R;
yr[3]=*VY_R;


FILE* f=fopen("koord.txt","w");
FILE* e=fopen("energy.txt","w");

eulerstep(T, yr , yn ,700,20000000,4,f,e);


fclose(f);
fclose(e);


return 0;
}
