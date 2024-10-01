#include "wigxjpf.h"
#include <iostream>
#include <cmath>
#include <iomanip>

double d_JMM(double M, double Mprime, double J, double beta);
double KlepsTo3j(double j1, double m1, double j2, double m2, double J, double M);


int main(){

    double M, Mprime, J, beta; // beta insert in degrees 

    std::cout << "enter M: "; std::cin >> M;
    std::cout << std::endl;
    std::cout << "enter Mprime: "; std::cin >> Mprime;
    std::cout << std::endl;
    std::cout << "enter J: "; std::cin >> J;
    std::cout << std::endl;
    std::cout << "enter beta (in degrees): "; std::cin >> beta;
    std::cout << std::endl;

    double printt = d_JMM(M, Mprime, J, beta);
    printf("Value of Wigner-d-function is: %f", printt);

    return 0;
}

double d_JMM(double M, double Mprime, double J, double b){
    double beta = b*3.14159/180; //converting degrees to radians
    double result = 0;
    double j1 = (J + Mprime)/2;
    double j2 = (J - Mprime)/2;

    for(double m1 = -j1; m1 < (j1 + 1); m1++){
        for(double m2 = -j2; m2 < (j2 + 1); m2++){
            if ((m1 + m2) != M){
                //std::cout << "Error: m1+m2 is not equal M" << std::endl;
                continue;
            }
            else{
                result += pow(-1, j2 + m2)*KlepsTo3j(j1, m1, j2, m2, J, M)*((pow(cos(beta/2), j1+j2+m1-m2)*pow(sin(beta/2), j1+j2-m1+m2))/pow(tgamma(j1+m1+1)*tgamma(j1-m1+1)*tgamma(j2+m2+1)*tgamma(j2-m2+1), 0.5));
            }
        }
    }

    return result*(pow((tgamma(j1 + j2 + J + 1 + 1)*tgamma(j1 + j2 - J + 1))/(2*J + 1), 0.5));
}

double KlepsTo3j(double j1, double m1, double j2, double m2, double J, double M){

    wig_table_init(2*100, 3);
    wig_temp_init(2*100);

    double val3j;
    val3j = wig3jj(int(2*j1), int(2*j2), int(2*J),
                   int(2*m1), int(2*m2), -int(2*M));

    wig_temp_free();
    wig_table_free();

    return (pow(-1, (-j1+j2-M)))*sqrt((2*J)+1)*val3j;
}