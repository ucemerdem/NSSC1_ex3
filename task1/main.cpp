#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <chrono>

using namespace std;

namespace program_options {

struct Options {
  string name;
  unsigned int iters;
  void print() const {
    printf("name: %s\n", name.c_str());
    printf("iters: %u\n", iters);
  }
};

auto parse(int argc, char *argv[]) {
  if (argc != 3)
    throw runtime_error("unexpected number of arguments");
  Options opts;
  opts.name = argv[1];
  if (sscanf(argv[2], "%u", &opts.iters) != 1 || opts.iters <= 0)
    throw runtime_error("invalid parameter for iters");
  return opts;
}

} // namespace program_options

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
        unsigned int getSize(){
            return n;
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

/// @brief overload the << operator for nice output of vectors
ostream& operator<< (ostream &os, const vector<double> &vec){
    for (unsigned int i = 0; i < vec.size(); i++){
        os << vec[i] << " ";
    }
    os << endl;
    return os;
}

/// @brief overload the - operator for vector subtraction
vector<double> operator- (vector<double> a, vector<double> b){
    vector<double> result(a.size());
    for (unsigned int i = 0; i < a.size(); i++){
        result[i] = a[i] - b[i];
    }
    return result;
}

/// @brief overload the * operator for vector multiplication (inner product)
double operator* (vector<double> a, vector<double> b){
    double result = 0.0;
    for (unsigned int i = 0; i < a.size(); i++){
        result += a[i] * b[i];
    }
    return result;
}

/// @brief overload the * operator for scalar-vector multiplication
vector<double> operator* (double a, vector<double> b){
    vector<double> result(b.size());
    for (unsigned int i = 0; i < b.size(); i++){
        result[i] = a * b[i];
    }
    return result;
}

/// @brief overload the + operator for vector addition
vector<double> operator+ (const vector<double>& a, const vector<double>& b){
    vector<double> result(a.size());
    for (unsigned int i = 0; i < a.size(); i++){
        result[i] = a[i] + b[i];
    }
    return result;
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
    vector<vector<double>> matrix(n, vector<double>(n, 0.0));
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

/// @brief This function performs a matrix-vector multiplication using the CCS_symm format.
/// @param ccs The CCS_symm object representing the matrix.
/// @param x The vector to be multiplied with the matrix.
/// @return The resulting vector of the matrix-vector multiplication.
vector<double> MVmultCCS_symm(CCS_symm& ccs, vector<double>& x){
    // initialize the result vector with zeros
    vector<double> result(ccs.getSize(), 0.0);

    // iterate over all columns
    for (unsigned int j = 0; j < ccs.getSize(); j++){
        // retrieve the indices for vals and row_idx for the current column
        // and perform the matrix-vector multiplication
        for (unsigned int i = ccs.getColPtr()[j]; i < ccs.getColPtr()[j+1]; i++){
            result[ccs.getRowIdx()[i]] += ccs.getVals()[i] * x[j];
            // check not to double-count the diagonal
            if (ccs.getRowIdx()[i] != j){
                // transposed multiplication
                result[j] += ccs.getVals()[i] * x[ccs.getRowIdx()[i]];
            }
        }
    }

    return result;
}

/// @brief This function calculates the norm of a vector with respect to the A-norm of a matrix stored in CCS format.
/// @param vec The vector for which the norm is to be calculated.
/// @param ccs The matrix based on which the norm is to be calculated, stored in CCS format.
/// @return The norm of the vector with respect to the A-norm of the matrix.
double getVectorANorm(vector<double>& vec, CCS_symm& ccs){
    vector<double> result = MVmultCCS_symm(ccs, vec);
    return sqrt(result * result);
}

/// @brief Calculate the solution x to Ax=b via the non-preconditioned conjugate gradient method, using a symmetric matrix stored in CCS format
/// @param ccs The symmetric, positive definite system matrix stored in CCS format
/// @param x0 The initial guess
/// @param b The right-hand side of the matrix equation Ax=b
/// @param max_iters The maximum number of iterations
/// @return The solution vector
vector<double> executeCG(CCS_symm& ccs, vector<double>& x0, vector<double>& b, unsigned int max_iters){
    // r is the residual vector (i.e. the difference between the right-hand side and the current approximation)
    // p is the search direction for the next iteration
    // x is the current approximation
    // Ap is the matrix-vector product of A and p used in the CG method
    // alpha and beta are the coefficients used in the CG method (i.e. the step sizes)
    // r_old_sq is the dot-product of the residual vector from the previous iteration
    // r_norm and r_norm_old are norms used as stopping criteria
    // EPS is the machine epsilon used as a threshold for the stopping criterion
    vector<double> r = b - MVmultCCS_symm(ccs, x0);
    vector<double> p = r;
    vector<double> x = x0;
    vector<double> Ap(ccs.getSize(), 0.0);
    double alpha, beta, r_old = r * r, r_norm, r_norm_old = sqrt(r * r);
    double EPS = 1e-6;
    unsigned int iterations = 0;
    
    // These are needed for plotting the convergence of the CG method via the Python script
    vector<double> r0 = b - MVmultCCS_symm(ccs, x0);
    vector<double> expected_solution = vector<double>(ccs.getSize(), 1.0);
    vector<double> x_diff;

    while(iterations < max_iters && iterations < ccs.getSize()){
        Ap = MVmultCCS_symm(ccs, p);
        alpha = (r * r) / (p * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        r_norm = sqrt(r * r);
        // FIXME: reconsider whether the stopping criterion is correct
        if (r_norm < EPS || abs(r_norm - r_norm_old) < EPS){
            break;
        }
        beta = (r * r) / r_old;
        p = r + beta * p;
        r_norm_old = r_norm;
        r_old = r * r;
        iterations++;

        x_diff = expected_solution - x;
        cout << r_norm / sqrt(r0 * r0) << " ";
        cout << getVectorANorm(x_diff, ccs) << endl;
    }

    // cout << "Number of iterations: " << iterations << endl;
    
    // cout << "Norm of last residual divided by norm of initial residual: " << r_norm / sqrt(r0 * r0) << endl;
    
    // cout << "Error on A-norm: " << getVectorANorm(x_diff, ccs) << endl;

    return x;
}

int main(int argc, char *argv[]) try{
    auto opts = program_options::parse(argc, argv);

    // do a test run of the CCS_symm class, using a randomly generated symmetric matrix
    // vector<vector<double>> matrix_test = generateSymmetricMatrix(6, 0.25);
    // CCS_symm ccs_test(matrix_test);

    // cout << "TEST MATRIX: " << endl;
    // cout << matrix_test << endl;
    // cout << ccs_test << endl;

    // obtain the matrix from file BCSSTK13 in the Matrix Market format
    vector<vector<double>> matrix_mtx = parseMTX(opts.name);
    CCS_symm ccs_mtx(matrix_mtx);

    // cout << "MTX parsing and conversion to CCS done:" << endl;
    // cout << "Number of values: " << ccs_mtx.getVals().size() << endl;
    // cout << "Number of row indices: " << ccs_mtx.getRowIdx().size() << endl;
    // cout << "Number of column pointers: " << ccs_mtx.getColPtr().size() << endl;

    // TEST MATRIX
    // prescribe the right-hand side b to the given solution x=[1,...,1]
    // vector<double> b;
    // vector<double> x(ccs_test.getSize(), 1.0);

    // b = MVmultCCS_symm(ccs_test, x);

    // cout << endl << "Multiplying the above test matrix with a vector of ones ..." << endl;
    // cout << "Result: " << b << endl;

    // excecute CG method using the previously calculated b and the initial guess x0=[0,...,0]
    // cout << "Executing the CG method for the test matrix ..." << endl;
    
    // vector<double> x0(ccs_test.getSize(), 0.0);
    // vector<double> result = executeCG(ccs_test, x0, b, opts.iters);

    // cout << "Result: " << result << endl;


    // MTX MATRIX
    // prescribe the right-hand side b to the given solution x=[1,...,1]
    vector<double> x(ccs_mtx.getSize(), 1.0);
    vector<double> b = MVmultCCS_symm(ccs_mtx, x);

    // cout << endl << "Multiplying BCSSTK13.mtx matrix with a vector of ones ..." << endl;

    // excecute CG method using the previously calculated b and the initial guess x0=[0,...,0]
    // cout << "Executing the CG method for BCSSTK13 ..." << endl;
    
    vector<double> x0(ccs_mtx.getSize(), 0.0);
    vector<double> result = executeCG(ccs_mtx, x0, b, opts.iters);

    // cout << "Result: " << result << endl;

  return EXIT_SUCCESS;
} catch (exception &e) {
  cout << e.what() << endl;
  return EXIT_FAILURE;
}
