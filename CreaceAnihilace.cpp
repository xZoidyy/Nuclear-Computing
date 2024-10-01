#include <iostream>
#include <cmath>
#include <vector>

class hilbert {
    private:
        std::vector<bool> state;
        int phase;

    public:
        std::pair<std::vector<bool>, int> Creation(const std::vector<bool>& vec, int where) {

            std::vector<bool> newBits = vec;
            int number_of_l = 0;
            int Phase;

            // make sure that my vectors of bits is of correct lenght
            if (vec.size() >= (where+1)){}
            else{
                for(int i = 0; i < ((where+1) - vec.size()); i++){
                    newBits.push_back(0);
                }
            }

            // Number of l in a bit string before (where) bit
            for(int j = 0; j < where; j++){
                if(newBits[j] == 1){
                    number_of_l++;
                }
            }

            if(newBits[where] == 0){
                newBits[where] = 1;
                Phase = int(pow(-1, number_of_l));
            }

            else{
                newBits[where] = 0;
                Phase = 0;
            }

            state = newBits;
            phase = Phase;

            return std::make_pair(state, phase);
        }

        std::pair<std::vector<bool>, int> Anihilation(const std::vector<bool>& vec, int where) {

            std::vector<bool> newBits = vec;
            int number_of_l = 0;
            int Phase;

            // make sure that my vectors of bits is of correct lenght
            if (vec.size() >= (where+1)){}
            else{
                for(int i = 0; i < ((where+1) - vec.size()); i++){
                    newBits.push_back(0);
                }
            }

            // Number of l in a bit string before (where) bit
            for(int j = 0; j < where; j++){
                if(newBits[j] == 1){
                    number_of_l++;
                }
            }

            if(newBits[where] == 1){
                newBits[where] = 0;
                Phase = int(pow(-1, number_of_l));
            }

            else{
                newBits[where] = 1;
                Phase = 0;
            }

            state = newBits;
            phase = Phase;

            return std::make_pair(state, phase);
        }
};

// print my state in hilbert space
void printState(const std::vector<bool>& vec){
    for (int i = 0; i < vec.size(); i++){
        std::cout << int(vec[i]) << "  ";
    }
    std::cout << std::endl;
}


int main(){

    hilbert Hilbert;

    std::vector<bool> test = {0, 0, 1, 1, 1, 0, 1, 1}; // 8 terms 

    std::vector<bool> newState;
    int phase;

    std::tie(newState, phase) = Hilbert.Anihilation(test, 3); // position in state from 0!
    
    printState(newState);
    std::cout << "Phase = " << phase << std::endl;

    return 0;
}