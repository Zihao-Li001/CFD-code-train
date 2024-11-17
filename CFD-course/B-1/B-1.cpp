#include <iostream>
#include <math.h>

int main()
{
    // theoretical solve
    double theo;
    theo = exp (1)*cos(1)+exp(1)*sin(1);
    printf("theoretical solve is %f.\n",theo);
    std::cin.get();
    
    // define \delta x
    double h[3] ={0.1,0.05,0.0025}; 
    double u[3];
    double dudx1,dudx2;
    double Error1[3],Error2[3];
////////////////////////////////
//      x-1    x    x+1       //
//       0     1     2        //
////////////////////////////////
    double x = 1.0;

    //1st upwind scheme
    //double upwind[3]={0.0};
    for(int i=0;i<3;i++)
    {
        u[0] = exp(x-h[i])*sin(x-h[i]);
        u[1] = exp(x)*sin(x);
        dudx1 = (u[1]-u[0])/h[i];
        // Calculate Error
        Error1[i] = abs(
            (dudx1-theo)
        );
        printf("Δx = %f, 1st-order upwind error is %.8f\n",h[i],Error1[i]);
    }
    std::cin.get();

    //2nd central scheme
    for(int i=0;i<3;i++)
    {
        u[0] = exp(x-h[i])*sin(x-h[i]);
        // u[1] = exp(x)*sin(x);
        u[2] = exp(x+h[i])*sin(x+h[i]);
        dudx2 = (u[2]-u[0])/(2.0*h[i]);
        // Calculate Error
        Error2[i] = abs(
            (dudx2-theo)
        );
        printf("Δx = %f, 2nd-order central error is %.8f\n",h[i],Error2[i]);

    }
    std::cin.get();

    return 0;
}