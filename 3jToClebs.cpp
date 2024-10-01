#include "wigxjpf.h"
#include <iostream>
#include <cmath>
#include <vector>

double KlepsTo3j(double j1, double m1, double j2, double m2, double J, double M);

int main(){

    wig_table_init(2*100, 3);
    wig_temp_init(2*100);

    double j1, m1, j2, m2, J, M;
    
    //insert values
    std::cout << "Insert values of j1, m1, j2, m2, J and M" << std::endl;
    std::cin >> j1;
    std::cin >> m1;
    std::cin >> j2;
    std::cin >> m2;
    std::cin >> J;
    std::cin >> M;

    std::cout << "Clebs is equal = " << KlepsTo3j(j1, m1, j2, m2, J, M) << std::endl;

    wig_temp_free();
    wig_table_free();

    return 0;
}

double KlepsTo3j(double j1, double m1, double j2, double m2, double J, double M){

    double val3j;
    val3j = wig3jj(int(2*j1), int(2*j2), int(2*J),
                   int(2*m1), int(2*m2), -int(2*M));

    return (pow(-1, (-j1+j2-M)))*sqrt((2*J)+1)*val3j;
}
