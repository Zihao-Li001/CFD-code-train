#include <iostream>
#include <math.h>

const int dt = 1;  // time step
const int N_step = 100;

typedef struct
{   
    int flag;
    int N_X=4;
    int N_Y=3;
    double l_x,l_y; // length of x, y
    double dx;  // cell size in X direction
    double dy;  // cell size in Y direction   
    double V[2]; //convection velocity
    double* u; // varibles [DOF]
    double** x; // position [DOF][2]
    double ALPHA = 1.0;
}CALC_POINTS;

typedef struct{
    double** mat;
    double*  rhs;
}LIN_SYS;

void initialize(CALC_POINTS* points, int n_x, int n_y, double l_x, double l_y)
{
    points->N_X = n_x;
    points->N_Y = n_y;
    points->l_x =  l_x;
    points->l_y =  l_y;
    points->dx  =   l_x / (n_x);
    points->dy  =   l_y / (n_y);
    points->ALPHA   =   1.0;
    points->V[0]    =   1.0;
    points->V[1]    =   1.0;

    int DOF = n_x * n_y;

    points->x   =   (double**)malloc(DOF * sizeof(double*));
    for (int i=0;i<DOF;i++)
    {
        points->x[i]    =   (double*)malloc(2 * sizeof(double));
    }

    // for(int i=0;i<Nc;i++)
    //     {
    //         if(i==5||i==6) 
    //             cell[i].flag = 1;  // point 6 or 7
    //         else 
    //             cell[i].flag = 0;  // boundary
    //         cell[i].U[0] = 0;   // u
    //         cell[i].U[1] = 0;   // v
    //     }
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