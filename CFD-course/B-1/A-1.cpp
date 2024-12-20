#include <iostream>
#include <math.h>

int main()
{
    // theoretical solve
    double theo;
    theo = exp(1)*cos(1)+exp(1)*sin(1);
    double x1,x2;
    x1 = cos(1)+sin(1);
    x2 = exp(1) * x1;
    printf("theoretical solve is %10f.\n",theo);
    std::cin.get();
    
    // define \delta x
    double h[3] ={0.1, 0.05,0.025}; 
    double u[6];
    double dudx1,dudx2,dudx3;
    double Error[3];
////////////////////////////////
//      x-1    x    x+1       //
//       0     1     2        //
////////////////////////////////
    double x = 1.0;

    // // 1st upwind scheme
    // double upwind[3]={0.0};
    // for(int i=0;i<3;i++)
    // {
    //     u[0] = exp(x-h[i])*sin(x-h[i]);
    //     u[1] = exp(x)*sin(x);
    //     dudx1 = (u[1]-u[0])/h[i];
    //     // Calculate Error
    //     Error[i] = abs(
    //         (dudx1-theo)
    //     );
    //     printf("Δx = %f, 1st-order upwind error is %.10f\n",h[i],Error[i]);
    // }
    // std::cin.get();

    //2nd central scheme
    for(int i=0;i<3;i++)
    {
        u[0] = exp(x-h[i])*sin(x-h[i]);
        // u[1] = exp(x)*sin(x);
        u[2] = exp(x+h[i])*sin(x+h[i]);
        dudx2 = (u[2]-u[0])/(2.0*h[i]);
        // Calculate Error
        Error[i] = abs(
            (dudx2-theo)
        );
        printf("Δx = %f, 2nd-order central error is %.10f\n",h[i],Error[i]);

    }
    std::cin.get();

    //5st central scheme
    for(int i=0;i<3;i++)
    {
        u[0] = exp(x+h[i])*sin(x+h[i]);
        u[1] = exp(x+2*h[i])*sin(x+2*h[i]);
        u[2] = exp(x+3*h[i])*sin(x+3*h[i]);
        u[3] = exp(x-h[i])*sin(x-h[i]);
        u[4] = exp(x-2*h[i])*sin(x-2*h[i]);
        u[5] = exp(x-3*h[i])*sin(x-3*h[i]);        
        double a1, a2, a3;
        a1 = 3.0 / 4.0;
        a2 = -3.0 / 20.0 ;
        a3 = 1.0 / 60.0;
        dudx3 = ( a1*(u[0]-u[3]) + a2*(u[1]-u[4]) + a3*(u[2]-u[5]) ) / h[i];
        // Calculate Error
        Error[i] = abs(
            (dudx3-theo)
        );
        printf("Δx = %f, 5th-order central error is %.10f\n",h[i],Error[i]);
    }
    std::cin.get();

    return 0;
}
