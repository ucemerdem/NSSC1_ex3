#include <eigen3/Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <vector>

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
    unsigned int iterations = 0;
    
    // These are needed for plotting the convergence of the CG method via the Python script
    vector<double> r0 = b - MVmultCCS_symm(ccs, x0);
    vector<double> expected_solution = vector<double>(ccs.getSize(), 1.0);

    while(iterations < max_iters && iterations < ccs.getSize()){
        Ap = MVmultCCS_symm(ccs, p);
        alpha = (r * r) / (p * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        r_norm = sqrt(r * r);

        beta = (r * r) / r_old;
        p = r + beta * p;
        r_norm_old = r_norm;
        r_old = r * r;
        iterations++;

        cout << r_norm / sqrt(r0 * r0) << endl;
    }

    // cout << "Number of iterations: " << iterations << endl;
    
    // cout << "Norm of last residual divided by norm of initial residual: " << r_norm / sqrt(r0 * r0) << endl;

    return x;
}

int main(int argc, char *argv[]) try{
    auto opts = program_options::parse(argc, argv);
    
    vector<vector<double>> matrix_mtx = parseMTX("bcsstk11.mtx");
    CCS_symm ccs_mtx(matrix_mtx);

    // own implementation of CG algorithm from task 1
    // prescribe the right-hand side b to the given solution x=[1,...,1]
    vector<double> x(ccs_mtx.getSize(), 1.0);
    vector<double> rhs = MVmultCCS_symm(ccs_mtx, x);

    // excecute CG method using the previously calculated b and the initial guess x0=[0,...,0]
    vector<double> x0(ccs_mtx.getSize(), 0.0);
    vector<double> result_1 = executeCG(ccs_mtx, x0, rhs, opts.iters);
    
    
    // Create a sparse matrix in Eigen format by reading from the CCS_symm object
    Eigen::SparseMatrix<double> A(ccs_mtx.getSize(), ccs_mtx.getSize());
    for (unsigned int j = 0; j < ccs_mtx.getSize(); j++) {
        for (unsigned int i = ccs_mtx.getColPtr()[j]; i < ccs_mtx.getColPtr()[j+1]; i++) {
            A.insert(ccs_mtx.getRowIdx()[i], j) = ccs_mtx.getVals()[i];
            if (ccs_mtx.getRowIdx()[i] != j) {
                A.insert(j, ccs_mtx.getRowIdx()[i]) = ccs_mtx.getVals()[i];
            }
        }
    }

    // prescribe the right-hand side b to the given solution x=[1,...,1]
    Eigen::VectorXd b(ccs_mtx.getSize());
    for (unsigned int i = 0; i < ccs_mtx.getSize(); i++) {
        b[i] = 1.0;
    }
    b = A * b;

    // Eigen’s Conjugate Gradients using DiagonalPreconditioner
    // Output the residual norm for each iteration (to be used for plotting by the Python script)
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg_dp;
    for(unsigned int i = 0; i < opts.iters; i++) {
        cg_dp.setMaxIterations(i);
        cg_dp.compute(A);
        Eigen::VectorXd x_dp = cg_dp.solve(b);
        cout << cg_dp.error() << endl;
    }

    // Eigen’s Conjugate Gradients using IncompleteCholesky
    // Output the residual norm for each iteration (to be used for plotting by the Python script)
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> cg_ic;
    for(unsigned int i = 0; i < opts.iters; i++) {
        cg_ic.setMaxIterations(i);
        cg_ic.compute(A);
        Eigen::VectorXd x_ic = cg_ic.solve(b);
        cout << cg_ic.error() << endl;
    }
    
    return EXIT_SUCCESS;
} catch (exception &e) {
  cout << e.what() << endl;
  return EXIT_FAILURE;
}
