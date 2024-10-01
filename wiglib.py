#include "./header/read_input.h"

int main(int argc, char const *argv[])
{   
    
    inputs input_parameters = read_input(); // class INPUTS with input parameters
    std::cout << input_parameters.d << std::endl;

    return 0;
}
