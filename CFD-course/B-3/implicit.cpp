#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int n_step = 100;  //number of iteration
const int n_x = 31;      //point number of x
const int n_y = 31;      //              of y
const double l_x = 1.0;  // lenth of x
const double l_y = 1.0;   //          y
const double alpha = 1.0; // diffusion coefficience
const double v[2] = {1.0, 1.0};     // convection velocity
const double dt = 0.1;      // time step
const double SOURCE_TERM = 1.0; //source term
static const int max_iteration = 10000; //maxmium iteration
static const double EPSILON = 1e-7;

//  define cell 
typedef struct
{
    int N_X;
    int N_Y;
    double l_x,l_y; // length of x, y
    double dx;  // cell size in X direction
    double dy;  // cell size in Y direction  
    double ALPHA; // diffusion coefficience   

    double V[2]; //convection velocity [DOF][2]
    double* u; // varibles [DOF]
    double** x; // position [DOF][2]
}CALC_POINTS;

//define determinant
typedef struct{
    double** mat;
    double*  rhs;
}LIN_SYS;

//allocate memory
void memory_allocation(
        CALC_POINTS* points, 
        LIN_SYS* ls,
        int DOF){
    points->u=(double*)malloc(DOF * sizeof(double));
    
    points->x=(double**)malloc(DOF * sizeof(double*));
    // points->V=(double**)malloc(DOF * sizeof(double*));
    for (int i=0;i<DOF;i++){
        points->x[i]=(double*)malloc(2 * sizeof(double));
        // points->V[i]=(double*)malloc(2 * sizeof(double));
    }

    ls->rhs=(double*)malloc(DOF * sizeof(double));

    ls->mat=(double**)malloc(DOF * sizeof(double*));
    for (int i=0;i<DOF;i++){
        ls->mat[i]=(double*)malloc(DOF * sizeof(double));
        if (ls->mat[i] == NULL) {
        fprintf(stderr, "Memory allocation failed for mat[%d]\n",i);
        exit(EXIT_FAILURE);
    }
    }
}

// free memory
void memory_free(
        CALC_POINTS* points, 
        LIN_SYS* ls,
        int DOF){
    free(points->u);
    
    // free memory layer by layer is IMPORTANT!
    for (int i=0;i<DOF;i++){
        free(points->x[i]);
        // free(points->V[i]);
    }
    free(points->x);
    // free(points->V);

    free(ls->rhs);

    for (int i=0;i<DOF;i++){
        free(ls->mat[i]);
    }
    free(ls->mat);
}

int index(
    int i,
    int j,
    int n_x){
    return(n_x*j + i);
}

void set_calc_points(
        CALC_POINTS* points, 
        int n_x, 
        int n_y, 
        double l_x, 
        double l_y){
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

// initiallization 
void initialize(
        CALC_POINTS* points, 
        int n_x, 
        int n_y, 
        double alpha){
    points->ALPHA   =   alpha;
    points->V[0] = 1.0;
    points->V[1] = 1.0;

    for (int j = 1; j < n_y-1; j++)
    {
        for(int i =1; i < n_x-1; i++)
        {
            int idx = index(i,j,n_x);
            points->u[idx]    = 0.0;
        }
    }
};

void set_matrix(
        CALC_POINTS* points, 
        LIN_SYS* ls,
        double dt  ){
    int n_x = points->N_X;
    int n_y = points->N_Y;
    double dx = points->dx;
    double dy = points->dy;
    double alpha = points->ALPHA;
    double conv_x = points->V[0];
    double conv_y = points->V[1];

    // set matrix coefficients
    for(int j=1; j<n_y-1; j++){
        for (int i=1; i<n_x-1; i++){
            int k = index(i,j,n_x);
            ls->mat[k][k-n_x] = -( conv_y/dy + alpha / (dy*dy) ) * dt;
            ls->mat[k][k-1  ] = -( conv_x/dx + alpha / (dx*dx) ) * dt;
            ls->mat[k][k    ] = 1 + ( conv_x/dx + 2.0 * alpha / (dx*dx) + conv_y/dy + 2.0 * alpha / (dy*dy) ) * dt;
            ls->mat[k][k+1  ] = - alpha * dt / (dx*dx);
            ls->mat[k][k+n_x] = - alpha * dt / (dy*dy);
            ls->rhs[k] = points->u[k] + SOURCE_TERM * dt;
            // printf("%f %f %f %f %f = %f\n",ls->mat[k][k-n_x],ls->mat[k][k-1  ],ls->mat[k][k    ],ls->mat[k][k+1  ],ls->mat[k][k+n_x],ls->rhs[k]);
        }
    }
}

void set_bc(
        CALC_POINTS* points,
        LIN_SYS* ls){
    int n_x = points->N_X;
    int n_y = points->N_Y;

    for( int j=0; j<n_y; j++ ){
        for( int i=0; i<n_x; i++ ){
            int k = index(i,j,n_x);

            if( i==0 || j==0 || i==n_x-1 || j==n_y-1){
                for( int l=0; l<n_x*n_y; l++){
                    ls->mat[k][l] = 0.0;
                    ls->mat[l][k] = 0.0;
                }
                ls->mat[k][k] = 1.0;
                ls->rhs[k] = 0.0;
            }
        }
    }
}

void mat_solve(double* sol, 
        LIN_SYS* ls,
        int max_iteration,
        int length){
    for (int count = 1; count < max_iteration; count++){
        for( int i=0; i<length; i++ ){
            sol[i] = ls->rhs[i];
            for( int j=0; j<length; j++){
                if( i != j){
                    sol[i] -= ls->mat[i][j]*sol[j]; 
                }
            }
            sol[i] /= ls->mat[i][i];         
        }

        double residual = 0.0;
        for( int i=0; i<length; i++){
            double r_i = -ls->rhs[i];
            for( int j=0; j<length; j++){
                r_i += ls->mat[i][j]*sol[j];
            }
            residual += r_i * r_i;
        }
        residual = sqrt(residual);
        printf("GS_loop %d: %e\n", count, residual);
        if( residual < EPSILON ) return;
    }
}

void write_vtk(
        CALC_POINTS* points,
        int file_num){
    FILE* fp;
    char filename[10000];
    snprintf(filename, 10000, "./data/out_%06d.vtk", file_num);
    fp = fopen( filename, "w" );

    if (!fp) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    int num = points->N_X * points->N_Y;
    fprintf(fp, "POINTS %d float\n", num);
    for( int i=0; i<num; i++){
        fprintf(fp, "%e %e %e\n",points->x[i][0],points->x[i][1], 0.0);
    }

    fprintf(fp, "CELLS %d %d\n", num, 2*num);
    for(int i=0; i<num; i++){
        fprintf(fp, "%d %d\n", 1, i);
    }

    fprintf(fp, "CELL_TYPES %d\n", num);
    for( int i=0; i<num; i++){
        fprintf(fp, "%d\n", 1);
    }

    fprintf(fp, "POINT_DATA %d\n", num);
    fprintf(fp, "SCALARS u float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");   
    for(int i=0; i<num; i++){
        fprintf(fp, "%e\n", points->u[i]);
    } 
    fclose(fp);
}

int main(
        int argc,
        char* argv[]) {
    CALC_POINTS points;
    LIN_SYS     ls;

    memory_allocation(&points,&ls,n_x*n_y);

    set_calc_points(&points, n_x, n_y, l_x, l_y);

    initialize(&points, n_x, n_y, alpha);

    for(int n=0; n < n_step; n++){
        set_matrix(&points, &ls, dt);
        set_bc(&points, &ls);
        mat_solve(points.u,&ls,max_iteration,n_x*n_y);
        write_vtk(&points, n);
    };

    memory_free(&points,&ls,n_x*n_y);
    return 0;
}