#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

const int dt = 0.0001;  // time step
const int t_step = 2;
const double SOURCE_TERM = 1.0;

typedef struct
{   
    int N_X;
    int N_Y;
    double l_x,l_y; // length of x, y
    double dx;  // cell size in X direction
    double dy;  // cell size in Y direction   
    // double** V; //convection velocity
    double* u; // varibles [DOF]
    // double** x; // position [DOF][2]
    double ALPHA = 1.0; // diffusion coefficience
}CALC_POINTS;

void initialize(CALC_POINTS* points, int n_x, int n_y, double l_x, double l_y, double alpha);
void Solve(int N_X, int N_Y, double dx, double dy, double v_x, double v_y,double alpha, double u[], double u_new[]);


int main() {
    CALC_POINTS points;
    
    // Parameters
    int n_x = 4;
    int n_y = 3;
    double l_x = 1.0;
    double l_y = 1.0;
    double alpha = 1.0;
    double V[2] = {1.0,1.0};

    double* u_new = (double*)malloc(n_x * n_y * sizeof(double));
    
    initialize(&points, n_x, n_y, l_x, l_y, alpha);
    for (int t=0;t<t_step;t++){
        Solve(points.N_X, points.N_Y, points.dx, points.dy, V[0], V[1], points.ALPHA, points.u, u_new);
        double* temp = points.u;
        points.u = u_new;
        u_new = temp;
    }

    cout << "Solution u at final time step:" << endl;
    for (int j = 0; j < points.N_Y; j++) {
        for (int i = 0; i < points.N_X; i++) {
            int idx = j * points.N_X + i;
            cout << points.u[idx] << " ";
        }
        cout << endl;
    }
    free(points.u);
    free(u_new);
    return 0;
}

void initialize(CALC_POINTS* points, int n_x, int n_y, double l_x, double l_y,double alpha)
{
    points->N_X = n_x;
    points->N_Y = n_y;
    points->l_x =  l_x;
    points->l_y =  l_y;
    points->dx  =   l_x / (n_x-1);
    points->dy  =   l_y / (n_y-1);
    points->ALPHA = alpha;
    int DOF = n_x * n_y;

    // points->x=(double**)malloc(DOF * sizeof(double*));
    // points->V=(double**)malloc(DOF * sizeof(double*));
    points->u=(double*)malloc(DOF*sizeof(double));
    // for (int i=0;i<DOF;i++)
    // {
    //     // points->x[i]=(double*)malloc(2 * sizeof(double));
    //     points->V[i]=(double*)malloc(2 * sizeof(double));
    // }

    for (int j = 0; j < n_y; j++)
    {
        for(int i =0; i < n_x; i++)
        {
            int idx = j* n_x + i;
            // points->x[idx][0] = i * points->dx;
            // points->x[idx][1] = j * points->dy;
            if (i==0 || i==n_x-1 || j==0 || j==n_y-1){
                points->u[idx] = 0.0;
            }
            else{
                points->u[idx] = 1.0;
            }
        }
    }

};

void Solve(int n_x, int n_y, double dx, double dy, double v_x, double v_y,double alpha, double u[], double u_new[])
{
    for(int j=1; j<n_y-1; j++)
    {            
        for(int i=1; i<n_x-1; i++)
        {
            int idx = j * n_x + i;
            int left= idx - 1;
            int right= idx + 1;
            int top = idx + n_x;
            int bottom = idx - n_x;
            
            double conv_x = (u[idx]-u[left]) / dx;
            double conv_y = (u[idx]-u[bottom]) / dy;
            double diff_x = (u[right]-2*u[idx]+u[left]) / (dx * dx);
            double diff_y = (u[top]-2*u[idx]+u[bottom]) / (dy * dy);
            u_new[idx]    =  u[idx] + dt*(-v_x * conv_x - v_y * conv_y+alpha*(diff_x+diff_y) + SOURCE_TERM);            
        }
    }
    for(int i=0; i<n_x ; i++){
        u_new[i] = 0.0;
        u_new[(n_y-1) * n_x + i] = 0.0;
    }
    for(int j=0; j<n_y; j++){
        u_new[j*n_x] = 0.0;
        u_new[j*n_x+n_x-1] = 0.0;
    }
};
