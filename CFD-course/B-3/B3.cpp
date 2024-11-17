#include <iostream>
#include <math.h>

#define Nc 12  // cell number
#define It 2  // iteration number

int dt = 1;  // time step
double dx = 1.0/3.0;  // cell size in X direction
double dy = 1.0/2.0;  // cell size in Y direction

// cell structer


//        
//
struct cell
{   
    int flag;
    double U[2];
};
struct cell cell[Nc];
int Nx = 4;
int Cellx=4, Celly=3;
// int i,j=0 // Nc = Nx*j +i

void initialsolve()
{
    for(int i=0;i<Nc;i++)
    {
        if(i==5||i==6) 
            cell[i].flag = 1;  // point 6 or 7
        else 
            cell[i].flag = 0;  // boundary
        cell[i].U[0] = 0;   // u
        cell[i].U[1] = 0;   // v
    }
};

void SetMesh()
{

};

void SolveImplicit()
{
    //dudu/dxdx
    for (int time=0;time<It;time++)
    {
            for(int i=0;i<4;i++)
            {
                
                for(int j=0;j<3;j++)
                {
                    if(cell[Nx*j+i].flag != 0)
                    {
                        double u = cell[Nx*j+i].U[0];
                        double v = cell[Nx*j+i].U[1];
                        
                    }
                }

            }
    }

};



int main()
{
    initialsolve();
    SetMesh();
    SolveImplicit();
    return 0;
}