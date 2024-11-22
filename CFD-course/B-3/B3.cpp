#include <iostream>
#include <math.h>
void initialize(CALC_POINTS* points, int n_x, int n_y, double l_x, double l_y, double alpha, double V[2]);
void SetMesh();
void Solve(CALC_POINTS* point);

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
    double** V; //convection velocity
    double* u; // varibles [DOF]
    double** x; // position [DOF][2]
    double ALPHA = 1.0; // diffusion coefficience
}CALC_POINTS;

// typedef struct{
//     double** mat;
//     double*  rhs;
// }LIN_SYS;



int main() {
    CALC_POINTS points;
    int n_x, n_y;
    double l_x, l_y;
    double alpha;
    double V[2];
    initialize(&points, n_x, n_y, l_x, l_y, alpha, V);
    Solve(CALC_POINTS* point);
    return 0;
}



void initialize(CALC_POINTS* points, int n_x, int n_y, double l_x, double l_y,double alpha, double V[2])
{
    points->N_X = n_x;
    points->N_Y = n_y;
    points->l_x =  l_x;
    points->l_y =  l_y;
    points->dx  =   l_x / (n_x);
    points->dy  =   l_y / (n_y);
    points->ALPHA   =   1.0;
    // points->V[0]    =   1.0;
    // points->V[1]    =   1.0;
    int DOF = n_x * n_y;
    points->x=(double**)malloc(DOF * sizeof(double*));
    points->V=(double**)malloc(DOF * sizeof(double*));
    points->u=(double*)malloc(DOF*sizeof(double));

    for (int i=0;i<DOF;i++)
    {
        points->x[i]=(double*)malloc(2 * sizeof(double));
        points->V[i]=(double*)malloc(2 * sizeof(double));
    }
    for (int j = 0; j < n_y; j++)
    {
        for(int i =0; i < n_x; i++)
        {
            int idx = j* n_y + i;
            points->x[idx][0] = i * points->dx;
            points->x[idx][1] = j * points->dy;
            points->u[idx]  = 0.0;
            points->V[idx][0] = V[0];
            points->V[idx][1] = V[1];
        }
    }
};

void Solve(CALC_POINTS* points)
{
    for (int it=0; it<N_step; it++)
    {
        for(int j=0;j<points->N_X ;j++)
        {            
            for(int i=0;i<3;i++)
            {
                int idx = j * points->N_Y + i;
                if (i>0 && i<points->N_X)
                double conv_x = (points->u[])
            }
        }
    }
};