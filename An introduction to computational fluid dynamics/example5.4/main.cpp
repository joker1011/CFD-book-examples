#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#define Num 5 //cell number

struct coefficient
{
    double p = 0;
    double e = 0;
    double w = 0;
    double ww = 0;
    double sp = 0;
};

void coefficient_calculate(double *Flux, double *D, double phiA, double phiB, coefficient *a)
{
    //ap
    a[0].p = 3*D[0] + D[1] + (7.0 / 8.0) * Flux[1];
    for (int i = 1; i < Num-1; i++){
        a[i].p = D[i] + D[i + 1] - (3.0 / 8.0) * Flux[i] + (3.0 / 4.0) * Flux[i + 1];
    }
    a[Num - 1].p = D[Num - 1] + 3 * D[Num] - (3.0 / 8.0) * Flux[Num - 1];

    //ae
    a[0].e = (1.0 / 3.0) * D[0] + D[1] - (3.0 / 8.0) * Flux[1];
    for (int i = 1; i < Num - 1; i++){
        a[i].e = D[i + 1] - (3.0 / 8.0) * Flux[i+1] ;
    }
    a[Num - 1].e = 0;

    //aw
    a[0].w = 0;
    a[1].w = D[1] + (7.0 / 8.0) * Flux[1] + (1.0 / 8.0) * Flux[2];
    for (int i = 2; i < Num - 1; i++){
        a[i].w = D[i] + (3.0 / 4.0) * Flux[i] + (1.0 / 8.0) * Flux[i + 1];
    }
    a[Num - 1].w = D[Num - 1] + (1.0 / 3.0) * D[Num] + (3.0 / 4.0) * Flux[Num - 1];

    //aww
    a[0].ww = 0;
    a[1].ww = 0;
    for (int i = 2; i <= Num - 1; i++){
        a[i].ww = (-1.0 / 8.0) * Flux[2];
    }

    //sp
    a[0].sp = ((8.0 / 3.0) * D[0] + Flux[0] + (1.0 / 4.0) * Flux[1]) * phiA;
    a[1].sp = (-1.0 / 4.0) * Flux[1] * phiA;
    a[Num - 1].sp = ((8.0 / 3.0) * D[Num] - Flux[Num]) * phiB;
}
void matrix_initialize(double (*A)[Num], coefficient* a)
{
    A[0][0] = a[0].p;
    A[0][1] = -a[0].e;// "-"!!!  ae is in the right side of the equation
    A[1][0] = -a[1].w;
    A[1][1] = a[1].p;
    A[1][2] = -a[1].e;
    for (int i = 2; i < Num - 1; i++){
        A[i][i - 2] = -a[i].ww;
        A[i][i - 1] = -a[i].w;
        A[i][i] = a[i].p;
        A[i][i+1] = -a[i].e;
    }
    A[Num-1][Num - 3] = -a[Num-1].ww;
    A[Num-1][Num - 2] = -a[Num-1].w;
    A[Num-1][Num-1] = a[Num-1].p;
}
void exact_solution()
{
    double phi_exact[Num] = { 0 };
    double deltaX = 1.0 / Num;
    double x = 0;
    for (int i = 0; i < Num; i++)
    {
        x = (0.5 + i) * deltaX;
        phi_exact[i] = 1 - (exp(2 * x) - 1) / (exp(2) - 1);
        std::cout << "phi " << i << " = " << phi_exact[i]<<std::endl;
    }
}
bool ifconverge(double* phi_new, double *phi,double rtol)
{
    for (int i = 0; i < Num; i++){
        if (fabs(phi_new[i] - phi[i]) >= rtol* phi_new[i]){
            return false;
        }
    }
    return true;
}
void GS_SOR(double(*A)[Num], coefficient* a, double* phi)
{
    //iteration parameter
    int ITERATION_LIMIT = 100;
    int i, j, step = 0;
    double relaxtion = 1;
    double rtol = 1e-5;
    //create new phi vector
    double phi_new[Num] = { 0 };
    //two summation in GS_SOR
    double sum1 = 0, sum2= 0;
    //convergence criteria
    bool converge = 0;

    //initialize new phi vector
    for (int k = 0; k < Num; k++){
        phi_new[k] = phi[k];
    }

    //GS_SOR method
    while (!converge){
        if (step >= ITERATION_LIMIT){
            return;
        }
        std::cout << "Iteration " << step << ": ";
        for (i = 0; i < Num; i++){
            std::cout << std::setw(8) << phi_new[i] << "	";
            for (j = 0; j < i ; j++){
                sum1 = sum1 + A[i][j] * phi_new[j];
            }
            for (j = i+1; j < Num; j++){
                sum2 = sum2 + A[i][j] * phi[j];
            }
            phi_new[i] = (1.0 - relaxtion) * phi[i] + (relaxtion / A[i][i]) * (a[i].sp - sum1 - sum2);
            sum1 = 0;
            sum2 = 0;
        }
        std::cout << std::endl;
        step++;
        converge = ifconverge(phi_new, phi, rtol);
        //upadte phi!!!!
        for (i = 0; i < Num; i++){
            phi[i] = phi_new[i];
        }
    }
}
int main(int argc,char *argv[])
{
    //cell center value
    double phi[Num] = { 0 };
    //boundary conditions
    double phiA = 1;
    double phiB = 0;

    //constant information
    double Flux[Num + 1] = { 0 };
    std::fill_n(Flux,Num+1, 0.2);
    double D[Num + 1] = { 0 };
    std::fill_n(D, Num+1, 0.5);
    //equation cofficients
    coefficient a[Num];

    //matrixA
    double A[Num][Num] = { 0 };
    coefficient_calculate(Flux, D, phiA, phiB, a);

    //output cofficients
    std::cout.precision(4);
    std::cout<<std::left;
    std::cout << "----------------------cofficient---------------------------" << std::endl;
    std::cout << std::setw(8) << "cell" << std::setw(8) << "aww" << std::setw(8) << "aw"
              << std::setw(8) << "ap"   << std::setw(8) << "ae"  << std::setw(8) << "sp" << std::endl;
    for (int i = 0; i < Num; i++){
        std::cout << std::setw(8) << i
                  << std::setw(8) << a[i].ww << std::setw(8) << a[i].w
                  << std::setw(8) << a[i].p  << std::setw(8) << a[i].e << std::setw(8) << a[i].sp << std::endl;
    }
    matrix_initialize(A, a);
    std::cout <<std::endl<< "------------------linear equation -------------------------"<<std::endl;
    for (int i = 0; i < Num; i++){
        std::cout <<"[";
        for (int j = 0; j < Num; j++){
            std::cout << std::setw(8)<< A[i][j];
        }
        std::cout << "] " ;
        std::cout << "[" << "phi" << i << "] ";
        std::cout << "[" <<  std::setw(5) << a[i].sp << "]" << std::endl;

    }
    std::cout << std::endl << "------------------iterative solution----------------------" << std::endl;

    clock_t start, end;
    start = clock();		//time start
    GS_SOR(A, a, phi);
    end = clock();		//time end
    double endtime = (double)(end - start) / CLOCKS_PER_SEC;
    std::cout << "Total time:" << endtime * 1000 << "ms" << std::endl;
    std::cout << std::endl << "------------------exact solution----------------------" << std::endl;
    exact_solution();
    return 0;
}