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
    double dudx;
    double Error[3];
////////////////////////////////
//       x    x+1    x+2       //
//       0     1      2        //
////////////////////////////////
    double x = 1.0;

    //3-rd order upwind scheme
    for(int i=0;i<3;i++)
    {
        u[0] = exp(x)*sin(x);
        u[1] = exp(x+h[i])*sin(x+h[i]);
        u[2] = exp(x+2*h[i])*sin(x+2*h[i]);
        dudx = (-3*u[0]+4*u[1]-u[2])/(2*h[i]);
        // Calculate Error
        Error[i] = abs(
            (dudx-theo)
        );
        printf("Δx = %f, 1st-order upwind error is %.8f\n",h[i],Error[i]);
    }
    std::cin.get();

    return 0;
}