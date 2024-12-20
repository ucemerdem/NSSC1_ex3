#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

using namespace std;

class CCS_symm{
    private:
        vector<double> vals;
        vector<unsigned int> row_idx;
        vector<unsigned int> col_ptr;
        unsigned int n;

    public:
        /// @brief Constructor of the CCS_symm class taking a symmetric matrix as input and converting it to the CCS_symm format.
        CCS_symm(vector<vector<double>> matrix){
            n = matrix.size();
            unsigned int nnz = 0;
            col_ptr.push_back(0); // Initialize the first column pointer to 0
            for (unsigned int j = 0; j < n; j++){
                for (unsigned int i = j; i < n; i++){
                    if (matrix[i][j] != 0){
                        vals.push_back(matrix[i][j]);
                        row_idx.push_back(i);
                        nnz++;
                }
            }
            col_ptr.push_back(nnz);
            }
        }

        vector<double> getVals(){
            return vals;
        }
        vector<unsigned int> getRowIdx(){
            return row_idx;
        }
        vector<unsigned int> getColPtr(){
            return col_ptr;
        }

        /// @brief overload the << operator for nice output of a CCS_symm object
        friend ostream& operator<< (ostream &os, const CCS_symm &ccs){
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

/// @brief overload the << operator for nice output of matrices
ostream& operator<< (ostream &os, const vector<vector<double>> &matrix){
    for (unsigned int i = 0; i < matrix.size(); i++){
        for (unsigned int j = 0; j < matrix[i].size(); j++){
            os << setw(4) << matrix[i][j] << " ";
        }
        os << endl << endl;
    }
    return os;
}

/// @brief This function parses a matrix stored in the Matrix Market format (.mtx).
/// @param filename The name of the file containing the matrix in the Matrix Market format.
/// @return A matrix in the form of a vector of vectors.
vector<vector<double>> parseMTX(const string& filename){
    ifstream reader(filename);
    string line;
    unsigned int rows, cols, entries;

    getline(reader, line); // skip the first line

    // extract the number of rows, columns, and the number of entries from the second line
    reader >> rows >> cols >> entries;

    vector<vector<double>> matrix(rows, vector<double>(cols, 0)); // 'rows' and 'cols' are interchangeable since the matrix is symmetric

    // continue with the third line until the end of the file to extract the non-zero entries
    unsigned int i, j;
    double val;
    for (int k = 0; k < entries; k++){
        reader >> i >> j >> val;
        matrix[i-1][j-1] = val;
        matrix[j-1][i-1] = val; // make the matrix symmetric
    }

    reader.close();

    return matrix;
}

/// @brief This function generates a random symmetric matrix of size n x n with a given percentage of nonzero entries.
/// @remark This function is for testing purposes only.
/// @param n The number of rows and columns of the matrix.
/// @param perc_non_zeros The percentage of nonzero entries in the matrix (a value between 0 and 1).
/// @return A symmetric matrix of size n x n with a given percentage of nonzero entries.
vector<vector<double>> generateSymmetricMatrix(unsigned int n, double perc_non_zeros=0.1){
    vector<vector<double>> matrix(n, vector<double>(n, 0));
    vector<double> nonzeros(int(ceil(perc_non_zeros * n * n)));

    // initialize random seed
    srand(chrono::steady_clock::now().time_since_epoch().count());
    // initialize the vector of nonzeros with random values
    for(unsigned int i = 0; i < nonzeros.size(); i++){
        nonzeros[i] = (rand() % 10 + 1) / 10.0;
    }

    unsigned int i, j;

    for(int k = 0; k < nonzeros.size(); k++){
        // generate random indices and check whether the generated pair of indices has already been used before to place a nonzero entry in the matrix
        while(true){
            i = int(rand() % n);
            j = int(rand() % n);
            if(matrix[i][j] == 0){
                matrix[i][j] = nonzeros[k];
                // if the previous line did not place an entry on the main diagonal line, place the same entry mirrored along the main diagonal (to make the matrix symmetric)
                if (i != j){
                    matrix[j][i] = nonzeros[k];
                    // skip one iteration to avoid placing too many entries in the matrix
                    k++;
                }
                break;
            }
        }
    }
    return matrix;
}

int main(){
    vector<vector<double>> matrix_test = generateSymmetricMatrix(5, 0.2);
    CCS_symm ccs_test(matrix_test);

    cout << matrix_test << endl;
    cout << ccs_test << endl;

    vector<vector<double>> matrix_mtx = parseMTX("bcsstk13.mtx");

    // cout << matrix_mtx << endl;

    CCS_symm ccs_mtx(matrix_mtx);

    cout << "MTX parsing and conversion to CCS done:" << endl;
    cout << "Number of values: " << ccs_mtx.getVals().size() << endl;
    cout << "Number of row indices: " << ccs_mtx.getRowIdx().size() << endl;
    cout << "Numbe rof column pointers: " << ccs_mtx.getColPtr().size() << endl;

    return 0;
}