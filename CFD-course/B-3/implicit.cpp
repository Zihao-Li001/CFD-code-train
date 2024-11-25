#include <iostream>
#include <math.h>

const int dt = 1;  // time step
const int N_step = 100;

typedef struct
{   
    int* flag;
    int N_X;
    int N_Y;
    double l_x,l_y; // length of x, y
    double dx;  // cell size in X direction
    double dy;  // cell size in Y direction  
    double ALPHA; // diffusion coefficience   

    double** V; //convection velocity
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

void initialize(
        CALC_POINTS* points, 
        int n_x, 
        int n_y, 
        double l_x, 
        double l_y, 
        double alpha, 
        double V[2]);

int main(
        int argc,
        char* argv[]
) {
    CALC_POINTS points;
    LIN_SYS     ls;
    int n_x, n_y;
    double l_x, l_y;
    double alpha;
    double V[2];
    memory_allocation(&points,&ls,n_x*n_y);
    initialize(&points, n_x, n_y, l_x, l_y, alpha, V);


    memory_free(&points,&ls,n_x*n_y);
    return 0;
}

// initiallization of cell
void initialize(CALC_POINTS* points, int n_x, int n_y, double l_x, double l_y,double alpha, double v[2])
{
    points->N_X = n_x;
    points->N_Y = n_y;
    points->l_x =  l_x;
    points->l_y =  l_y;
    points->dx  =   l_x / (n_x);
    points->dy  =   l_y / (n_y);
    points->ALPHA   =   alpha;

    for (int j = 0; j < n_y; j++)
    {
        for(int i =0; i < n_x; i++)
        {
            int idx = j* n_y + i;
            points->x[idx][0] = i * points->dx;
            points->x[idx][1] = j * points->dy;
            points->u[idx]  = 0.0;
            points->V[idx][0] = v[0];
            points->V[idx][1] = v[1];
            // bottom boundary
            if(i==0)
                points->flag[idx] = 0;
            // top boundary
            else if(i==(n_x-1))
                points->flag[idx] = 0;
        }
        // left boundary
        if(j==0)
            points->flag = 0;
        //right boundary
        else if(j==(n_y-1))
            points->flag = 0;
    }
};