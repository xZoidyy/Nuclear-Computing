#include "wigxjpf.h"
#include <iostream>
#include <cmath>
#include <vector>

double KlepsTo3j(double j1, double m1, double j2, double m2, double J, double M);

int main(){

    wig_table_init(2*100, 3);
    wig_temp_init(2*100);

    double j1, j2;
    
    std::cout << "Insert values of j1 and j2" << std::endl;
    std::cin >> j1;
    std::cin >> j2;
    std::cout << std::endl;

    // create MATRIX
    std::vector<std::vector<double> > matrix(m, std::vector<double>(n));

    // ORTONORMALITY TEST
    bool q = true;

    for (double J = abs(j1-j2); J < (j1+j2+1); J++){
        for (double Jprime = abs(j1-j2); Jprime < (j1+j2+1); Jprime++){
            for (double M = -J; M < (J+1); M++){
                for (double Mprime = -Jprime; Mprime < (Jprime+1); Mprime++){
                    double tmp = 0;
                    for (double m1 = -j1; m1 < (j1+1); m1++){
                        for (double m2 = -j2; m2 < (j2+1); m2++){
                            tmp += KlepsTo3j(j1, m1, j2, m2, J, M)*KlepsTo3j(j1, m1, j2, m2, Jprime, Mprime);
                        }
                    }
                    if (tmp == 1 && (J != Jprime || M != Mprime)){
                        std::cout << "Error: Relace orthonormality nefunguji!" << std::endl;
                        q = false;
                    }
                }
            }
        }
    }
    if (q){
        std::cout << "Relace orthonormality funguji!" << std::endl;
    }

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
