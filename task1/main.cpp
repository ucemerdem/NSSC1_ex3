#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

using namespace std;

class CCS{
    private:
        vector<double> vals;
        vector<unsigned int> row_idx;
        vector<unsigned int> col_ptr;
        unsigned int n;

    public:
        /// @brief Constructor of the CCS class taking a matrix as input and converting it to the CCS format.
        CCS(vector<vector<double>> matrix){
            n = matrix.size();
            unsigned int nnz = 0;
            col_ptr.push_back(0); // Initialize the first column pointer to 0
            for (unsigned int j = 0; j < n; j++){
                for (unsigned int i = 0; i < n; i++){
                    if (matrix[i][j] != 0){
                    vals.push_back(matrix[i][j]);
                    row_idx.push_back(i);
                    nnz++;
                }
            }
            col_ptr.push_back(nnz);
            }

            // TODO: adapt for symmetric matrices
        }

        /// @brief overload the << operator for nice output of a CCS object
        friend ostream& operator<< (ostream &os, const CCS &ccs){
            os << "Values: ";
            for (unsigned int i = 0; i < ccs.vals.size(); i++){
                os << ccs.vals[i] << " ";
            }
            os << endl;

            os << "Row indices: ";
            for (unsigned int i = 0; i < ccs.row_idx.size(); i++){
                os << ccs.row_idx[i] << " ";
            }
            os << endl;

            os << "Column pointers: ";
            for (unsigned int i = 0; i < ccs.col_ptr.size(); i++){
                os << ccs.col_ptr[i] << " ";
            }
            os << endl;

            return os;
        }
};

/// @brief This function parses a matrix stored in the Matrix Market format (.mtx).
/// @param filename The name of the file containing the matrix in the Matrix Market format.
/// @return A matrix in the form of a vector of vectors.
vector<vector<double>> parseMTX(string filename){
    // TODO: implement
}

/// @brief This function generates a random matrix of size n x n with a given percentage of nonzero entries.
/// @remark This function is for testing purposes only.
/// @param n The number of rows and columns of the matrix.
/// @param perc_non_zeros The percentage of nonzero entries in the matrix (a value between 0 and 1).
/// @return A matrix of size n x n with a given percentage of nonzero entries.
vector<vector<double>> generateMatrix(unsigned int n, double perc_non_zeros=0.1){
    vector<vector<double>> matrix(n, vector<double>(n, 0));
    vector<double> nonzeros(int(ceil(perc_non_zeros * n * n)));

    // initialize random seed
    srand(chrono::steady_clock::now().time_since_epoch().count());
    // initialize the vector of nonzeros with random values
    for(unsigned int i = 0; i < nonzeros.size(); i++){
        nonzeros[i] = (rand() % 10 + 1) / 10.0;
    }

    unsigned int i, j;

    for(double val : nonzeros){
        // generate random indices and check whether the generated pair of indices has already been used before to place a nonzero entry in the matrix
        while(true){
            i = int(rand() % n);
            j = int(rand() % n);
            if(matrix[i][j] == 0){
                matrix[i][j] = val;
                break;
            }
        }
    }
    return matrix;
}

int main(){
    vector<vector<double>> matrix = generateMatrix(5, 0.1);
    CCS ccs(matrix);

    for (unsigned int i = 0; i < matrix.size(); i++){
        for (unsigned int j = 0; j < matrix[i].size(); j++){
            cout << setw(3) << matrix[i][j] << " ";
        }
        cout << endl;
    }

    cout << ccs << endl;

    return 0;
}