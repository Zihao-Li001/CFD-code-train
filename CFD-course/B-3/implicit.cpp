#include <stdio.h>
#include <stdlib.h>

const int n_step = 100;
const int n_x = 4;
const int n_y = 3;
const double l_x =1.0;
const double l_y=1.0;
const double alpha = 1.0;
double v[2]={1.0, 1.0};
const double dt = 1.0;
const double SOURCE_TERM = 1.0;

typedef struct
{
    int N_X;
    int N_Y;
    double l_x,l_y; // length of x, y
    double dx;  // cell size in X direction
    double dy;  // cell size in Y direction  
    double ALPHA; // diffusion coefficience   

    int* flag={0};
    double** V; //convection velocity [DOF][2]
    double* u; // varibles [DOF]
    double** x; // position [DOF][2]
}CALC_POINTS;

typedef struct{
    double** mat;
    double*  rhs;
}LIN_SYS;

//allocate memory
void memory_allocation(
        CALC_POINTS* points, 
        LIN_SYS* ls,
        int DOF)
{
    points->flag=(int*)malloc(DOF * sizeof(int));

    points->u=(double*)malloc(DOF * sizeof(double));
    
    points->x=(double**)malloc(DOF * sizeof(double*));
    points->V=(double**)malloc(DOF * sizeof(double*));
    for (int i=0;i<DOF;i++){
        points->x[i]=(double*)malloc(2 * sizeof(double));
        points->V[i]=(double*)malloc(2 * sizeof(double));
    }

    ls->rhs=(double*)malloc(DOF * sizeof(double));

    ls->mat=(double**)malloc(DOF * sizeof(double*));
    for (int i=0;i<DOF;i++){
        ls->mat[i]=(double*)malloc(2 * sizeof(double));
    }
}

void memory_free(
        CALC_POINTS* points, 
        LIN_SYS* ls,
        int DOF)
{
    free(points->flag);

    free(points->u);
    
    for (int i=0;i<DOF;i++){
        free(points->x[i]);
        free(points->V[i]);
    }
    free(points->x);
    free(points->V);

    free(ls->rhs);

    for (int i=0;i<DOF;i++){
        free(ls->mat[i]);
    }
    free(ls->mat);
}

int index(
    int i,
    int j,
    int n_x
)
{
    return(n_x*j + i);
}

void set_calc_points(
        CALC_POINTS* points, 
        int n_x, 
        int n_y, 
        double l_x, 
        double l_y)
{
    points->N_X = n_x;
    points->N_Y = n_y;
    points->l_x =  l_x;
    points->l_y =  l_y;
    points->dx  =   l_x / (n_x-1);
    points->dy  =   l_y / (n_y-1);
    for (int j = 0; j < n_y; j++)
    {
        for(int i =0; i < n_x; i++)
        {
            int idx = index(i,j,n_x);
            points->x[idx][0] = i * points->dx;
            points->x[idx][1] = j * points->dy;
        }
    }
};

void set_matrix(
        CALC_POINTS* points, 
        LIN_SYS* ls,
        double dt
        ){
    int n_x = points->N_X;
    int n_y = points->N_Y;
    double dx = points->dx;
    double dy = points->dy;
    double alpha = points->ALPHA;
    double conv_x[n_x*n_y];
    double conv_y[n_x*n_y];
    for (int j = 1; j < n_y-1; j++)
    {
        for(int i =1; i < n_x-1; i++)
        {
            int idx = index(i,j,n_x);
            conv_x[idx] = points->V[idx][0];
            conv_y[idx] = points->V[idx][1];
        }
    }
    
    for(int j=1; j<n_y-1; j++){
        for (int i=1; i<n_x-1; i++){
            int k = index(i,j,n_x);
            ls->mat[k][k-n_x] = -( conv_y[k]/dy + alpha / (dy*dy) ) * dt;
            ls->mat[k][k-1  ] = -( conv_y[k]/dy + alpha / (dy*dy) ) * dt;
            ls->mat[k][k    ] = 1 + ( conv_y[k]/dy + alpha / (dy*dy) + conv_y[k]/dy + alpha / (dy*dy) ) * dt;
            ls->mat[k][k+1  ] = - alpha * dt / (dx*dx);
            ls->mat[k][k+n_x] = - alpha * dt / (dy*dy);

            ls->rhs[k] = points->u[k] + SOURCE_TERM * dt;
        }
    }
}

// initiallization 
void initialize(
        CALC_POINTS* points, 
        int n_x, 
        int n_y, 
        double alpha, 
        double v[2]
        )
{
    points->ALPHA   =   alpha;

    for (int j = 1; j < n_y-1; j++)
    {
        for(int i =1; i < n_x-1; i++)
        {
            int idx = index(i,j,n_x);
            points->u[idx]    = 1.0;
            points->V[idx][0] = v[0];
            points->V[idx][1] = v[1];
            points->flag[idx] = 1;
        }
    }
};

int main(
        int argc,
        char* argv[]
) {
    CALC_POINTS points;
    LIN_SYS     ls;
    memory_allocation(&points,&ls,n_x*n_y);
    set_calc_points(&points, n_x, n_y, l_x, l_y);
    initialize(&points, n_x, n_y, alpha, v);
    for(int n=0; n < n_step; n++){
        set_matrix(&points, &ls, dt);
    };
    memory_free(&points,&ls,n_x*n_y);
    return 0;
}