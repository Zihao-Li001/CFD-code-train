/*************
  1D inviscous compressible flow
  boundary condition ---    Dummy Cell
  flux calc          ---    AUSM Scheme
  Reconst            ---    MUSCL 
  Restrictor         ---    Van Albada
  discre T           ---    Runge-Kutta
 **************/

#include "stdio.h"
#include "conio.h"
#include "malloc.h"
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define h (1/400.0)  // mesh step
#define Nc 404  // total mesh number; h = 1/h +4
#define PI 3.1415927
#define It 1000
#define gama 1.4  // 气体比热比

double KAKA=0.0; // 限制器控制参数
double XS1,XS2; // 计算域的两个端点
double dt=2.5e-5; // 时间步长
double timesum; // 总的计算时间

// 函数声明
void output();
void SolveWtoU(double W[3],double U[3]); // W = conservation var 密度、动量密度、能量密度,
                                         // U = basic var  压力、速度和温度
void SolveUtoW(double W[3],double U[3]); // 基本变量和守恒变量之间的转换函数

struct cell
{
  int flag;
  double xc;
  double W[3],Wp[3]; //conser var
  double U[3];  //bas var
  double R[3];  // residuals
  double S; //entropy
};
struct cell cell[Nc];

// mesh generation && initialization
void initialsolve()
{
  int i;
  double x,xi,xe;
  XS1=-0.0;
  XS2=1.0;
  xi = XS1 - 2*h;
  xe = XS2 + 2*h;
  for(i=0;i<Nc;i++)
  {
    cell[i].xc = (xi+h/2) + i*h; // cell center location
    cell[i].flag = 0;
    x = cell[i].xc;
    if(x<0.5)
    {
      cell[i].U[0]=0.445; cell[i].U[1]=0.698;cell[i].U[2]=3.528;
    }
    else
    {
      cell[i].U[0]=0.5;cell[i].U[1]=0.0;cell[i].U[2]=0.571;
    }
    SolveUtoW(cell[i].W,cell[i].U);
  }
  printf("x=%f,rho=%f,u=%f,p=%f\n",
  cell[i].xc,cell[i].U[0],cell[i].U[1],cell[i].U[2]);
//  getchar();
};

// boundarycondition: Dummy Cell Method 
void boundarycondition()
{
  int i;
  double U[3];
  //赋值远场
  for(i=0;i<Nc;i++)
  {
    if(cell[i].flag == 1)
    {
      cell[i].U[0] = 0.445;
      cell[i].U[1] = 0.698;
      cell[i].U[2] = 3.528;
      SolveUtoW(cell[i].W,U);
    }
    else if(cell[i].flag == 2)
    {
      cell[i].U[0]=0.5;
      cell[i].U[1]=0.0;
      cell[i].U[2]=0.571;
      SolveUtoW(cell[i].W,U);
    }
  }
}

void SolveReconstruction(double Um1[],double UI[],double Up1[],double Up2[],double UL[],double UR[])
{
  int i;
  double epsilon,aR[3],aL[3],bL[3],bR[3],sL[3],sR[3];
   epsilon = 1.0e-5*h;
  for(i=0;i<3;i++)
  {
    aR[i] = Up2[i] - Up1[i]; //dealta + U(I+1)
    bR[i] = Up1[i] - UI[i];
    aL[i] = Up1[i] - UI[i];
    bL[i] = UI[i] - Um1[i];
  }
  for(i=1;i<3;i++)
  {
    sR[i] = (2.0*aR[i]*bR[i] + 1.0e-6)/(aR[i]*aR[i]+bR[i]*bR[i] + 1.0e-6);
    sL[i] = (2.0*aL[i]*bL[i] + 1.0e-6)/(aL[i]*aL[i]+bL[i]*bL[i] + 1.0e-6);
  }
  for(i=1;i<3;i++)
  {
    UL[i] = UI[i] + 0.25*sL[i]* ((1-KAKA*sL[i])*bL[i] + (1+KAKA*sL[i])*aL[i]);
    UR[i] = Up1[i] - 0.25*sR[i]* ((1-KAKA*sR[i])*aR[i] + (1+KAKA*sR[i])*bR[i]);
  }
}

// convert W to U  basic-->conserv
void SolveWtoU(double W[3],double U[3])
{
  U[0] = W[0]; //density
  U[1] = W[1]/W[0]; //velocity
  U[2] = (gama - 1.0)*(W[2]-0.5*U[0]*U[1]*U[1]); //e
}

// convert U to W  conserv-->basic
void SolveUtoW(double W[3],double U[3])
{
  W[0] = U[0];
  W[1] = U[1]*U[0];
  W[2] = U[2]/(gama - 1.0) + 0.5*U[0]*U[1]*U[1];
}


// AUSM scheme
void SolveAUSMFlux(double UL[],double UR[],double Fc[]) 
// left status, right status, flux
{
  double WL[3]={0.0},WR[3]={0.0};
  double deta,fMn;
  double ML,cl,MR,cr,Ml,Mr,Mn;
  double pl,pr,pn,PL,PR;
  deta = 0.25;
  SolveUtoW(WL,UL);
  SolveUtoW(WR,UR);
  
  PL = UL[2]; //calc pressure at left cell
  cl = sqrt(gama *PL/UL[0]);

  PR = UR[2]; //calc pressure at right cell
  cr = sqrt(gama *PR/UR[0]);

  ML = UL[1]/cl;
  MR = UR[1]/cr;

  //  calcu Ml && pl
  if(ML>=1.0)
  {
    Ml = ML;
    pl = PL;
  }
  else if(ML<-1.0)
  {
    Ml = 0.0;
    pl = 0.0;
  }
  else
  {
    Ml = 0.25*(ML+1)*(ML+1);
    pl = 0.25*PL*(ML+1)*(ML+1)*(2.0-ML);
  }

  // calcu Mr && pr
  if(MR<-1.0)
  {
    Mr = MR;
    pr = PR;
  }
  else if(MR>1.0)
  {
    Mr = 0.0;
    pr = 0.0;
  }
  else
  {
    Mr = 0.25*(MR-1)*(MR-1);
    pr = 0.25*PR*(MR-1)*(MR-1)*(2.0+MR);
  }

  Mn = Mr + Ml;
  pn = pr + pl;
  
  if(fabs(Mn) > deta) 
  {
    fMn = fabs(Mn);
  }
  else 
  {
    fMn = 0.5*(fabs(Mn)*fabs(Mn)+deta*deta)/deta;
  }

  Fc[0] = 0.5*Mn*(WL[0]*cl+WR[0]*cr)-0.5*fMn*(WR[0]*cr-WL[0]*cl);
  Fc[1] = 0.5*Mn*(WL[1]*cl+WR[1]*cr)-0.5*fMn*(WR[1]*cr-WL[1]*cl)+pn;
  Fc[2] = 0.5*Mn*((WL[2]+PL)*cl+(WR[2]+PR)*cr)-0.5*fMn*((WR[2]+PR)*cr-(WL[2]+PL)*cl);
}

void SolveResidual()
{
  int i,j;
  double UL[3],UR[3],Fp[3],Fm[3];
  for(i=0;i<Nc;i++)
  {
    if(cell[i].flag == 0)
    {
      //i+1/2
      SolveReconstruction(cell[i-1].U, cell[i].U, cell[i+1].U, cell[i+2].U, UL, UR);
      SolveAUSMFlux(UL,UR,Fp); // second order
      SolveAUSMFlux(cell[i].U,cell[i+1].U,Fp);  //first order
      
      //i-1/2
      SolveReconstruction(cell[i-2].U, cell[i-1].U, cell[i].U, cell[i+1].U, UL, UR);
      SolveAUSMFlux(UL,UR,Fm); // second order
      SolveAUSMFlux(cell[i-1].U,cell[i].U,Fm);  //first order      
      
      for(j=0;j<3;j++)
      {
        cell[i].R[j] = -(Fp[j]-Fm[j])/h;
      }
    }
  }
}

//update flow field
void SolveNextstep(double ar[],int ir)
{
  int i,j;
  for(i=0;i<=Nc;i++)
  {
    if(cell[i].flag == 0)
    {
      if(ir == 0)
      {
         for(j=0;j<3;j++)
            cell[i].Wp[j] = cell[i].W[j];
      }
      for(j=0;j<3;j++)
        cell[i].W[j] = cell[i].Wp[j] + ar[ir]*dt*cell[i].R[j];
      SolveWtoU(cell[i].W,cell[i].U);
    }
  }
}

//Runge-Kutta method 
void SolveRungeKutta()
{
  //
  double ar[4]={1/4.0, 1/3.0, 0.5, 1.0};
  
  int it,ir;
  timesum = 0.0;
  for(it=0;;it++)
  {
    for(ir=0;ir<4;ir++)
    {
      boundarycondition();
      SolveResidual();
      SolveNextstep(ar,ir);
    }
    timesum+=dt;
    printf("it=%d,timesum=%f\n",it,timesum);
    if(timesum>=0.16)
    {
      output();
      printf("Please enter any key to continue...");
      getchar();
      break;
    }
  }
}

void output()
{
     int i;
     FILE *fpd,*fpu,*fpp;

     // open density file
     if((fpd=fopen("Density.plt","w"))==NULL)
     {
    	printf("connot open infile");
    	return;

     }	    
    fprintf(fpd,"TITLE = \"TestCase\"n");
    fprintf(fpd,"VARIABLES = \"x\", \"density\"n");
    fprintf(fpd,"ZONE T=\"Only Zone\", I=%d, F=POINTN",Nc);

    // open pressure file
    if((fpp=fopen("Pressure.plt","w"))==NULL)
    {
        printf("connot open infile");
        return;
    }
    fprintf(fpp,"TITLE = \"TestCase\"n");
    fprintf(fpp,"VARIABLES = \"x\", \"pressure\"n");
    fprintf(fpp,"ZONE T=\"Only Zone\", I=%d, F=POINTN",Nc);

    // open velocity file
    if((fpu=fopen("Velocity.plt","w"))==NULL)
    {
        printf("connot open infile");
        return;
    }
    fprintf(fpu,"TITLE = \"TestCase\"n");
    fprintf(fpu,"VARIABLES = \"x\", \"velocity\"n");
    fprintf(fpu,"ZONE T=\"Only Zone\", I=%d, F=POINTN",Nc);

    for(i=0;i<Nc;i++)
    {
        fprintf(fpd,"%f %fn",cell[i].xc,cell[i].U[0]);
        fprintf(fpu,"%f %fn",cell[i].xc,cell[i].U[1]);
        fprintf(fpp,"%f %fn",cell[i].xc,cell[i].U[2]);
    
        printf("x=%f,d=%f,u=%f,p=%f\n",cell[i].xc,cell[i].U[0],cell[i].U[1],cell[i].U[2]);
        
    //    getchar();
    }
    fclose(fpd);

    fclose(fpp);

    fclose(fpu);
}

void SetMeshFlag()
{
  int i;
  for(i=0;i<Nc;i++)
  {
    if(cell[i].xc>=XS1 && cell[i].xc<=XS2)
      cell[i].flag=0;
    else if(cell[i].xc>XS2)
      cell[i].flag = 2;
    else
      cell[i].flag = 1;
    printf("cell[%d].flag=%d\n",i,cell[i].flag);
    getchar();
  }
}

void SolveEntropy()
{
  int i;
  for(i=0;i<Nc;i++)
  {
    cell[i].S=cell[i].U[2]/pow(cell[i].U[0],gama);
  }
}

int main(void)
{
  initialsolve();
  SetMeshFlag();
  SolveRungeKutta();
  output();
  return 0;
}
