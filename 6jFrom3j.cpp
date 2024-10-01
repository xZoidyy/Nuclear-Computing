#include "wigxjpf.h"
#include <iostream>
#include <cmath>
#include <iomanip>

double x3jTo6j(double j1, double j2, double j3, double j12, double j23, double J);

int main(){

    double val6j;
    double prec = pow(1, -10);

    wig_table_init(2*100, 3);
    wig_temp_init(2*100);

    double j1, j2, j3, j12, j23, J;
    
    std::cout << "Insert values of j1, j2, j12, j3, J and j23" << std::endl;
    std::cin >> j1;
    std::cin >> j2;
    std::cin >> j12;
    std::cin >> j3;
    std::cin >> J;
    std::cin >> j23;
    std::cout << std::endl;

    val6j = wig6jj(2*j1, 2*j2, 2*j12,
                   2*j3, 2*J, 2*j23);

    double a = val6j;
    double b = x3jTo6j(j1, j2, j3, j12, j23, J); //because my function is sooooo presice
    std::cout << std::fixed << std::setprecision(10) << a << std::endl;
    std::cout << std::fixed << std::setprecision(10) << b << std::endl;

    if (abs(a-b) < prec){
        std::cout << "The relation between 6j-symbols to 3j-symbols works!" << std::endl;
        std::cout << a << " = " << b << std::endl;
    }
    else{
        std::cout << "ERROR: The relation between 6j-symbols to 3j-symbols does not works!" << std::endl;
        std::cout << a << " != " << b << std::endl;
    }

    wig_temp_free();
    wig_table_free();

    return 0;
}

double x3jTo6j(double j1, double j2, double j3, double j12, double j23, double J){
    double val;

    wig_table_init(2*100, 3);
    wig_temp_init(2*100);

    for (double m1 = -j1; m1 < (j1+1); m1++){
        for (double m2 = -j2; m2 < (j2+1); m2++){
            for (double m3 = -j3; m3 < (j3+1); m3++){
                for (double m12 = -j12; m12 < (j12+1); m12++){
                    for (double m23 = -j23; m23 < (j23+1); m23++){
                        for (double M = -J; M < (J+1); M++){
                            val += (pow(-1, (j3 + J + j23 - m3 - M - m23)))*wig3jj(int(2*j1), int(2*j2), int(2*j12),
                                                                                   int(2*m1), int(2*m2), int(2*m12))*
                                                                            wig3jj(int(2*j1), int(2*J), int(2*j23),
                                                                                   int(2*m1), -int(2*M), int(2*m23))*
                                                                            wig3jj(int(2*j3), int(2*j2), int(2*j23),
                                                                                   int(2*m3), int(2*m2), -int(2*m23))*
                                                                            wig3jj(int(2*j3), int(2*J), int(2*j12),
                                                                                  -int(2*m3), int(2*M), int(2*m12));
                        }
                    }
                }
            }
        }
    }

    wig_temp_free();
    wig_table_free();

    return val;
}