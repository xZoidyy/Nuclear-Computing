#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

// basis template
class basis{
    public:
        int n;
        int l;
        int j; // 2*j
        int m; // 2*m
        int t_z; // 2* t_z

        // Constructor
        basis(int a, int b, int c, int d, int e) : n(a), l(b), j(c), m(d), t_z(e) {}

        // Define the equality operator for checking for duplicates 
        bool operator==(const basis& other) const {
            return n == other.n && l == other.l && j == other.j && m == other.m && t_z == other.t_z;
        }
    };

bool areVectorsEqual(const std::vector<basis>& vec1, const std::vector<basis>& vec2); // compare my basis vector to not be same at the end
void print_vec(std::vector<basis> vec); // print vector basis
void print_vec_vec(std::vector<std::vector<basis>> vec_vec); // print vector of vector basis which where chossen to have total M

void generate_combinations_recursive(const std::vector<basis>& input_vector, int num_to_choose, int start_idx, std::vector<basis>& current_combination, std::vector<std::vector<basis>>& all_combinations);
std::vector<std::vector<basis>> generate_combinations(const std::vector<basis>& input_vector, int num_to_choose);

// How many nucleons will be in my shell (protons and neutrons)
std::vector<int> number_of_choosing (int A, int Z, int n);

// Define a function which creates a vector of basis
std::vector<basis> create_basis(int n_max, int nucleon);

// Define function with takes desired M of nucleus with A nucleons and Z protons, generetes it in m-scheme
std::vector<std::vector<basis>> m_scheme_basis_forM (std::vector<std::vector<basis>> vec_p, std::vector<std::vector<basis>> vec_n, int M, int Nprotons, int Nneutrons);

bool correctShell(int A, int Z, int n);

int main(){

    // parameters of nucleus and totoal projection
    int A = 6; // Mass number
    int Z = 4;  // Proton number
    int M = -4; // My total 2*M!! what I want 

    int n = 1; // main quantum number n

    if (correctShell(A, Z, n)){

        //create vector of basis 
        std::vector<basis> my_basis_p = create_basis(n, -1); // -1 for protons
        std::vector<basis> my_basis_n = create_basis(n, 1); // 1 for neutrons

        // print my basis
        std::cout << "My neutron basis: "<< std::endl;
        print_vec(my_basis_n);
        std::cout << std::endl;
        std::cout << "My proton basis: "<< std::endl;
        print_vec(my_basis_p);
        std::cout << std::endl;
        std::cout << "#################################################################################" << std::endl;
        std::cout << "#################################################################################" << std::endl;
        std::cout << std::endl;

        std::vector<int> num_to_choose = number_of_choosing (A, Z, n); // Calculate How many nucleons are in my shell

        // make an m_sceme vector for desirable value of M of nucleus (A,Z)
        std::vector<std::vector<basis>> vec_p = generate_combinations(my_basis_p, num_to_choose[0]); // generate all combinations of states in shell for protons
        std::vector<std::vector<basis>> vec_n = generate_combinations(my_basis_n, num_to_choose[1]); // generate all combinations of states in shell for neutrons
        std::vector<std::vector<basis>> m_scheme = m_scheme_basis_forM(vec_p, vec_n, M, num_to_choose[0], num_to_choose[1]); // combinations of :basis_p and basis_n, M projection, Number of protons and neutrons in shell
        print_vec_vec(m_scheme); // printing combination of basis coresponding to M
        std::cout << "dim = " << m_scheme.size() << std::endl;

        return 0;
    }

    else{
        std::cout << "Error: In my chosen shell level n are not any nucleons (empty shell), therefore dim = 0." << std::endl;
        return 1;
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool correctShell(int A, int Z, int n){ // checks if shell level has a meaning or its just empty
    int N = A - Z;
    int number = (n + 1)*(n + 2);
    int total = 0;

    for (int k = 0; k <= n; k++){
        total += (k + 1)*(k + 2);
    }

    if ((Z <= (total - number)) && (N <= (total - number))){
        return false;
    }
    else{
        return true;
    }
}

bool areVectorsEqual(const std::vector<basis>& vec1, const std::vector<basis>& vec2) {
    // Check if vectors are of the same size
    if (vec1.size() != vec2.size()) {
        return false;
    }
    // Compare vectors element-wise
    for (size_t i = 0; i < vec1.size(); ++i) {
        if (!(vec1[i] == vec2[i])) {
            return false;
        }
    }
    return true;
}

std::vector<basis> create_basis(int n_max, int nucleon){ // nucleon: -1 for proton and 1 for neutron

    int number_of_states = (n_max + 1)*(n_max + 2);

    std::vector<basis> vec_basis;

    for(int l = n_max; l >= 0; l -= 2){ // l = n, n-2, n-4, ..., 1 or 0

        int j_up = int(2*(l+0.5));
        int j_down = int(2*fabs(l-0.5));

        if (j_up == j_down){
            for (int j_1 = j_down; j_1 > -j_down -1; j_1 -= 2){
                basis obj(n_max, l, j_down, j_1, nucleon);
                vec_basis.push_back(obj);
            }
        }
        else{
            for (int j_1 = j_down; j_1 > -j_down -1; j_1 -= 2){
                basis obj(n_max, l, j_down, j_1, nucleon);
                vec_basis.push_back(obj);
            }
            for (int j_2 = j_up; j_2 > -j_up -1; j_2 -= 2){
                basis obj(n_max, l, j_up, j_2, nucleon);
                vec_basis.push_back(obj);
            }
        }
    }
    return vec_basis;
}

void generate_combinations_recursive(const std::vector<basis>& input_vector, int num_to_choose, int start_idx, std::vector<basis>& current_combination, std::vector<std::vector<basis>>& all_combinations) {
    if (num_to_choose == 0) {
        all_combinations.push_back(current_combination);
        return;
    }

    // If choosing all variables, use the entire input_vector as one combination
    if (num_to_choose == input_vector.size()) {
        all_combinations.push_back(input_vector);
        return;
    }

    for (int i = start_idx; i < input_vector.size(); ++i) {
        current_combination.push_back(input_vector[i]);
        generate_combinations_recursive(input_vector, num_to_choose - 1, i + 1, current_combination, all_combinations);
        current_combination.pop_back();
    }
}

std::vector<std::vector<basis>> generate_combinations(const std::vector<basis>& input_vector, int num_to_choose) {
    std::vector<std::vector<basis>> all_combinations;
    std::vector<basis> current_combination;
    generate_combinations_recursive(input_vector, num_to_choose, 0, current_combination, all_combinations);
    return all_combinations;
}

std::vector<int> number_of_choosing (int A, int Z, int n){
    
    std::vector<int> vector;

    int N = A - Z; // number of neutrons

    int Nneutrons = N;
    int Nprotons = Z;

    for(int k = 0; k > -1; k++){
        int Number_in_levels = (k + 1)*(k + 2);
        if (Nneutrons < Number_in_levels){
            break;
        }
        else{
            Nneutrons = Nneutrons - Number_in_levels;
        }
    }
    for(int k = 0; k > -1; k++){
        int Number_in_levels = (k + 1)*(k + 2);
        if (Nprotons < Number_in_levels){
            break;
        }
        else{
            Nprotons = Nprotons - Number_in_levels;
        }
    }

    if (Nprotons == 0){ //maximum occupation of shell
        vector.push_back((n + 1)*(n + 2));
    }
    else{
        vector.push_back(Nprotons);
    }

    if (Nneutrons == 0){ //maximum occupation of shell
        vector.push_back((n + 1)*(n + 2));
    }
    else{
        vector.push_back(Nneutrons);
    }

    return vector;
}

std::vector<std::vector<basis>> m_scheme_basis_forM (std::vector<std::vector<basis>> vec_p, std::vector<std::vector<basis>> vec_n, int M, const int Nprotons, const int Nneutrons){

    std::vector<std::vector<basis>> m_scheme_vector;

    for (int i = 0; i < vec_n.size(); i++){ // going thru all combinations of vector basis for protons where I already chose how many neutrons are there
        for (int j = 0; j < vec_p.size(); j++){ // same for protons as neutrons
            std::vector<basis> tmp_vec;
            std::vector<basis> tmp_vec_n;
            std::vector<basis> tmp_vec_p;
            int my_M = 0;

            for (int k1 = 0; k1 < Nneutrons; k1++){
                my_M = my_M + (vec_n[i])[k1].m;
                tmp_vec_n.push_back((vec_n[i])[k1]);  
            }

            for (int k2 = 0; k2 < Nprotons; k2++){
                my_M = my_M + (vec_p[j])[k2].m;
                tmp_vec_p.push_back((vec_p[j])[k2]);
            }

            // sort this for better comparison of values of basis vectors
            std::sort(tmp_vec_n.begin(), tmp_vec_n.end(), [](const basis& a, const basis& b) {
                if (a.l != b.l) {
                    return a.l < b.l;
                }
                if (a.j != b.j) {
                    return a.j < b.j;
                }
                return a.m < b.m;
            });

            for (int h1 = 0; h1 < Nneutrons; h1++){
                tmp_vec.push_back(tmp_vec_n[h1]); // put in one vector where will be neutrons and protons simulataniouselly
            }

            // sort this for better comparison of values of basis vectors
            std::sort(tmp_vec_p.begin(), tmp_vec_p.end(), [](const basis& a, const basis& b) {
                if (a.l != b.l) {
                    return a.l < b.l;
                }
                if (a.j != b.j) {
                    return a.j < b.j;
                }
                return a.m < b.m;
            });

            for (int h2 = 0; h2 < Nprotons; h2++){
                tmp_vec.push_back(tmp_vec_p[h2]);
            }

            if (my_M == (M)){ // chech if my total M is what I want

                if (m_scheme_vector.size() < 1){
                    m_scheme_vector.push_back(tmp_vec);
                }
                else{
                    bool verita = false;
                    for (int k = 0; k < m_scheme_vector.size(); k++){
                        if (areVectorsEqual(tmp_vec, m_scheme_vector[k])){
                            verita = true;
                        }
                    }
                    if (!verita){
                        m_scheme_vector.push_back(tmp_vec);
                    }
                }
            }
        }
    }
    
    return m_scheme_vector;
}

void print_vec(std::vector<basis> vec){
    for (int i = 0; i < vec.size(); i++){
        if (vec[i].m > 0){
            std::cout << "n = " << vec[i].n << "\t" << "l = " << vec[i].l << "\t" << "2*j = " << vec[i].j << "  " << "2*m =  " << vec[i].m << "  " << "2*t_z = " << vec[i].t_z << std::endl;
        }
        else{
            std::cout << "n = " << vec[i].n << "\t" << "l = " << vec[i].l << "\t" << "2*j = " << vec[i].j << "  " << "2*m = " << vec[i].m << "  " << "2*t_z = " << vec[i].t_z << std::endl;
        }
    }
}

void print_vec_vec(std::vector<std::vector<basis>> vec_vec){
    for(int i = 0; i < vec_vec.size(); i++){
        for(int j = 0; j < vec_vec[i].size(); j++){
            if ((vec_vec[i])[j].m > 0){
                std::cout << "n = " << vec_vec[i][j].n << "\t" << "l = " << vec_vec[i][j].l << "\t" << "2*j = " << vec_vec[i][j].j << "  " << "2*m =  " << vec_vec[i][j].m << "  " << "2*t_z = " << vec_vec[i][j].t_z << std::endl;
            }
            else{
                std::cout << "n = " << vec_vec[i][j].n << "\t" << "l = " << vec_vec[i][j].l << "\t" << "2*j = " << vec_vec[i][j].j << "  " << "2*m = " << vec_vec[i][j].m << "  " << "2*t_z = " << vec_vec[i][j].t_z << std::endl;
            }
        }
        std::cout << std::endl;
    }
}