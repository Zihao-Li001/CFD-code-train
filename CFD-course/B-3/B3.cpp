#include <iostream>
#include <math.h>
void initialize(CALC_POINTS* points, int n_x, int n_y, double l_x, double l_y);
void SetMesh();
void Solve();

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
    double ALPHA = 1.0; // diffusion coefficience
}CALC_POINTS;

typedef struct{
    double** mat;
    double*  rhs;
}LIN_SYS;



int main() {
    CALC_POINTS points;
    int n_x, n_y;
    double l_x, l_y;
    initialize(&points, n_x, n_y, l_x, l_y);
    Solve();
    return 0;
}



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
    points->u = (double*)malloc(DOF*sizeof(double));
    for (int j = 0; j < n_y; j++)
    {
        for(int i =0; i < n_x; i++)
        {
            int idx = j* n_y + i;
            points->x[idx][0] = i * points->dx;
            points->x[idx][1] = j * points->dy;
            points->u[idx]  = 0.0;
        }
    }
};

void Solve()
{
    for (int it=0; it<N_step; it++)
    {
        for(int i=0;i<4;i++)
        {            
            for(int j=0;j<3;j++)
            {

            }
        }
    }
};