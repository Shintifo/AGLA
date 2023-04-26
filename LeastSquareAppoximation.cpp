// Timofey Brayko
// t.brayko@innopolis.university
#include <iostream>
#include <cmath>

using namespace std;

class IdentityMatrix;

class ColumnVector {
private:
    double *data;
    int size;
public:
    ColumnVector(int n) {
        data = new double[n];
        size = n;
    }

    ColumnVector() {
        size = 0;
    }

    //Put a value at certain position
    void put(double value, int index) {
        data[index] = value;
    }

    void operator=(ColumnVector v) {
        if (size == 0) {
            size = v.size;
            data = new double[size];
        }
        for (int i = 0; i < size; ++i) {
            data[i] = v.data[i];
        }
    }

    /**
     * Swaps two rows.
     * @param first index of first row.
     * @param second index of second row.
     */
    void swapRows(int first, int second) {
        double temp = data[first];
        data[first] = data[second];
        data[second] = temp;
    }

    /**
     * Subtract rows
     * @param multiple multiplier of subtrahend.
     * @param fRow index of first row.
     * @param sRow index of second row.
     */
    void subtract(double multiple, int fRow, int sRow) {
        if (data[sRow] != 0) {
            double res = (data[fRow] + multiple * data[sRow]);
            data[fRow] = res;
        }
    }

    //Divide the row by particular value
    void division(double value, int index) {
        if (data[index] != 0 && value != 0) {
            double res = data[index] / value;
            data[index] = res;
        }
    }

    //Returns value at certain position
    double valueAt(int index) {
        return data[index];
    }

    int getSize() {
        return size;
    }

    ColumnVector operator+(ColumnVector &matrix) {
        ColumnVector result(size);
        for (int i = 0; i < size; ++i) {
            result.data[i] = this->data[i] + matrix.data[i];
        }
        return result;
    }

    ColumnVector operator-(ColumnVector &matrix) {
        ColumnVector result(size);
        for (int i = 0; i < size; ++i) {
            result.data[i] = this->data[i] - matrix.data[i];
        }
        return result;
    }

    double norm() {
        double sum = 0.0;
        for (int i = 0; i < size; i++) {
            sum += data[i] * data[i];
        }
        return sqrt(sum);
    }
};

istream &operator>>(istream &in, ColumnVector &vector) {
    double value;
    for (int i = 0; i < vector.getSize(); ++i) {
        cin >> value;
        vector.put(value, i);
    }
    return in;
}

ostream &operator<<(ostream &out, ColumnVector &vector) {
    out.precision(4);
    for (int j = 0; j < vector.getSize(); ++j) {
        out << fixed << vector.valueAt(j) << endl;
    }
    return out;
}

class Matrix {
private:
    int rows;
    int columns;
    double **matrix{};
    string name;

public:
    Matrix(int n, int m, string name) {
        matrix = new double *[n];
        for (int i = 0; i < n; ++i) {
            matrix[i] = new double[m];
        }
        rows = n;
        columns = m;
        this->name = name;
    }

    Matrix() {
        rows = 0;
        columns = 0;
    }

    string getName() {
        return name;
    }

    void setName(string name) {
        this->name = name;
    }


    virtual void operator=(const Matrix &a) {
        this->rows = a.rows;
        this->columns = a.columns;
        matrix = new double *[rows];
        for (int i = 0; i < rows; ++i) {
            matrix[i] = new double[columns];
        }
        for (int i = 0; i < this->rows; ++i) {
            for (int j = 0; j < this->columns; ++j) {
                matrix[i][j] = a.matrix[i][j];
            }
        }
    }

    /**
     * Performs addition of two matrices
     * @param a Second matrix
     * @return result of summation
     */
    virtual Matrix operator+(Matrix &a) const {
        if (this->rows == a.rows && this->columns == a.columns) {
            Matrix newMatrix(rows, columns, name);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < columns; ++j) {
                    newMatrix.matrix[i][j] = matrix[i][j] + a.matrix[i][j];
                }
            }
            return newMatrix;
        } else {
            return {0, 0, 0};
        }
    }

    /**
     * Subtract two matrices
     * @param a Second matrix
     * @return result of difference
     */
    virtual Matrix operator-(Matrix a) const {
        if (this->rows == a.rows && this->columns == a.columns) {
            Matrix newMatrix(rows, columns, name);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < columns; ++j) {
                    newMatrix.matrix[i][j] = matrix[i][j] - a.matrix[i][j];
                }
            }
            return newMatrix;
        } else {
            return {0, 0, 0};
        }
    }

    /**
     * Multiply two matrices
     * @param a Second Matrix
     * @return Result of multiplication two matrices
     */
    virtual Matrix operator*(const Matrix a) const {
        if (columns == a.rows) {
            Matrix newMatrix(rows, a.columns, name);
            for (int k = 0; k < rows; ++k) {
                for (int i = 0; i < a.columns; ++i) {
                    double sum = 0;
                    for (int j = 0; j < columns; ++j) {
                        sum += matrix[k][j] * a.matrix[j][i];
                    }
                    newMatrix.matrix[k][i] = sum;
                }
            }
            return newMatrix;
        } else {
            return {0, 0, name};
        }
    }

    virtual ColumnVector operator*(ColumnVector b) const {
        ColumnVector vector(this->rows);
        for (int i = 0; i < this->rows; ++i) {
            double sum = 0;
            for (int j = 0; j < columns; ++j) {
                sum += matrix[i][j] * b.valueAt(j);
            }
            vector.put(sum, i);
        }
        return vector;
    }

    /**
     * Transposes matrix
     * @return transposed matrix
     */
    Matrix T() {
        Matrix transposed(columns, rows, "transpose");
        for (int i = 0; i < columns; ++i) {
            for (int j = 0; j < rows; ++j) {
                transposed.matrix[i][j] = matrix[j][i];
            }
        }
        return transposed;
    }

    /**
     * Shows number of rows
     * @return
     */
    virtual int getRows() const {
        return rows;
    }

    /**
     * Returns number of columns
     * @return
     */
    virtual int getColumns() const {
        return columns;
    }

    /**
     * Returns a value inside matrix at certain position
     * @param x Row number
     * @param y Column number
     * @return Value at position (x,y)
     */
    virtual double getAtPosition(int x, int y) {
        return matrix[x][y];
    }

    /**
     * Puts value into matrix at certain position
     * @param value Value of double type
     * @param x Row number
     * @param y Column number
     */
    virtual void putValueAt(double value, int x, int y) {
        matrix[x][y] = value;
    }
};

istream &operator>>(istream &in, Matrix &matrix) {
    double value;
    for (int i = 0; i < matrix.getRows(); ++i) {
        for (int j = 0; j < matrix.getColumns(); ++j) {
            in >> value;
            matrix.putValueAt(value, i, j);
        }
    }
    return in;
}

ostream &operator<<(ostream &out, Matrix &m) {
    out.precision(4);
    if (m.getRows() != 0) {
        for (int i = 0; i < m.getRows(); ++i) {
            for (int j = 0; j < m.getColumns(); ++j) {
                out << fixed << m.getAtPosition(i, j);
                if (j != m.getColumns() - 1) {
                    out << " ";
                }
            }
            out << endl;
        }
    } else {
        out << "Error: the dimensional problem occurred" << endl;
    }
    return out;
}

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int n, string name) : Matrix(n, n, name) {}

    SquareMatrix() : Matrix() {}

    SquareMatrix operator+(SquareMatrix &matrix) {
        Matrix res = Matrix::operator+(matrix);
        const SquareMatrix &result = static_cast<const SquareMatrix &> (res);
        return result;
    }

    SquareMatrix operator-(SquareMatrix &matrix) {
        Matrix res = Matrix::operator-(matrix);
        const SquareMatrix &result = static_cast<const SquareMatrix &> (res);
        return result;
    }

    SquareMatrix operator*(SquareMatrix &matrix) {
        Matrix res = Matrix::operator*(matrix);
        const SquareMatrix &result = static_cast<const SquareMatrix &> (res);
        return result;
    }

    ColumnVector operator*(ColumnVector &vector) {
        ColumnVector res = Matrix::operator*(vector);
        return res;
    }

    SquareMatrix T() {
        const Matrix &bR = Matrix::T();
        const SquareMatrix &res = static_cast<const SquareMatrix &>(bR);
        return res;
    }

    IdentityMatrix inverseMatrix();

    void ZeroLowerTriangle(IdentityMatrix &augmented, int &counterOfOperations);

    void ZeroUpperTriangle(IdentityMatrix &augmented, int &counterOfOperations);

    void DiagonalNormalization(IdentityMatrix &augmented);

    ColumnVector solveEquation(ColumnVector &vector) {
        return RREF(vector);
    }

    ColumnVector RREF(ColumnVector &vector);

    void DiagonalNormalization(ColumnVector &vector);

    void ZeroLowerTriangle(ColumnVector &vector, int &counterOfOperations);

    void ZeroUpperTriangle(ColumnVector &vector, int &counterOfOperations);

};

class IdentityMatrix : public SquareMatrix {
public :
    IdentityMatrix(int n,string name) : SquareMatrix(n,name) {
        for (int i = 0; i < this->getRows(); ++i) {
            for (int j = 0; j < this->getColumns(); ++j) {
                if (i == j) this->putValueAt(1, i, i);
                else this->putValueAt(0, i, j);
            }
        }
    }

    Matrix operator+(Matrix &matrix) {
        Matrix res = Matrix::operator+(matrix);
        return res;
    }

    Matrix operator-(Matrix &matrix) {
        Matrix res = Matrix::operator-(matrix);
        return res;
    }

    SquareMatrix operator*(SquareMatrix &matrix) {
        return matrix;
    }

    virtual void operator=(const Matrix &a) {
        Matrix::operator=(a);
    }
};

class EliminationMatrix : public IdentityMatrix {
public:
    EliminationMatrix(int row, int column, Matrix &matrix) : IdentityMatrix(matrix.getRows(), matrix.getName()) {
        double value = (-1) * matrix.getAtPosition(row, column) / matrix.getAtPosition(column, column);
        this->putValueAt(value, row, column);
    }

    SquareMatrix operator*(SquareMatrix &matrix) {
        Matrix res = Matrix::operator*(matrix);
        const SquareMatrix &result = static_cast<const SquareMatrix &> (res);
        return result;
    }

};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(int n, int firstRow, int secondRow, string name) : IdentityMatrix(n, name) {
        this->putValueAt(0, firstRow, firstRow);
        this->putValueAt(0, secondRow, secondRow);
        this->putValueAt(1, firstRow, secondRow);
        this->putValueAt(1, secondRow, firstRow);
    }

    Matrix operator*(Matrix &matrix) {
        Matrix res = Matrix::operator*(matrix);
        return res;
    }

    SquareMatrix operator*(SquareMatrix &matrix) {
        Matrix res = Matrix::operator*(matrix);
        const SquareMatrix &result = static_cast<const SquareMatrix &> (res);
        return result;
    }
};

void fillMatrixVector(Matrix &A, int n, Matrix &input, ColumnVector &b) {
    for (int i = 0; i < A.getRows(); ++i) {
        A.putValueAt(1, i, 0);
        b.put(input.getAtPosition(i, 1), i);
    }
    for (int i = 1; i <= n; ++i) {
        for (int j = 0; j < A.getRows(); ++j) {
            double val = pow(input.getAtPosition(j, 0), i);
            A.putValueAt(val, j, i);
        }
    }
}

int main() {
    int m, n;
    cin >> m;
    Matrix input(m, 2,"input");
    cin >> input;
    cin >> n;
    Matrix A(m, n + 1, "A");
    ColumnVector b(m);
    fillMatrixVector(A, n, input, b);

    cout << "A:" << endl << A;
    Matrix A_T = A.T();
    Matrix temp = A_T * A;
    SquareMatrix A_TA = static_cast<SquareMatrix &> (temp);
    A_TA.setName("A_T*A");
    cout << "A_T*A:" << endl << A_TA;
    SquareMatrix inverse = A_TA.inverseMatrix();
    inverse.setName("inverse");
    cout << "(A_T*A)^-1:" << endl << inverse;

    ColumnVector A_T_b = A_T * b;
    cout << "A_T*b:" << endl << A_T_b;
    ColumnVector ans = inverse * A_T_b;
    cout << "x~:" << endl << ans;

    return 0;
}

ColumnVector SquareMatrix::RREF(ColumnVector &vector) {
    int counterOfOperations = 1;
    cout << "step #0:" << endl;
    cout << *this << vector;

    ZeroLowerTriangle(vector, counterOfOperations);

    ZeroUpperTriangle(vector, counterOfOperations);

    cout << "Diagonal normalization:" << endl;
    DiagonalNormalization(vector);

    cout << *this << vector;
    return vector;
}


void printA(SquareMatrix& matrix, IdentityMatrix& au) {
    for (int i = 0; i < matrix.getRows(); ++i) {
        for (int j = 0; j < matrix.getColumns(); ++j) {
            cout << matrix.getAtPosition(i,j) << " ";
        }
        cout << "    ";
        for (int j = 0; j < matrix.getColumns(); ++j) {
            cout << au.getAtPosition(i,j) << " ";
        }
        cout << endl;
    }
    cout << endl;
}

IdentityMatrix SquareMatrix::inverseMatrix() {
    int counterOfOperations = 1;
    IdentityMatrix augmented(this->getRows(), "Aug");

    ZeroLowerTriangle(augmented, counterOfOperations);
//    printA(*this, augmented);

    ZeroUpperTriangle(augmented, counterOfOperations);
//    printA(*this, augmented);
    DiagonalNormalization(augmented);
//    printA(*this, augmented);
    return augmented;
}

void SquareMatrix::ZeroLowerTriangle(IdentityMatrix &augmented, int &counterOfOperations) {
    for (int i = 0; i < this->getColumns() - 1; ++i) {
        double maximum = abs(this->getAtPosition(i, i));
        int maxIndex = i;
        //Finding max pivot
        for (int j = i + 1; j < this->getRows(); ++j) {
            if (abs(this->getAtPosition(j, i)) > maximum) {
                maximum = abs(this->getAtPosition(j, i));
                maxIndex = j;
            }
        }
        //Permutation
        if (maxIndex != i) {
            PermutationMatrix P(this->getRows(), i, maxIndex, "P");
            *this = P * *this;
            augmented = P * augmented;
        }
        //Elimination
        for (int j = i + 1; j < this->getRows(); ++j) {
            if (this->getAtPosition(j, i) != 0) {
                EliminationMatrix E(j, i, *this);
                *this = E * *this;
                augmented = E * augmented;
            }
        }
//        printA(*this, augmented);
    }
}

void SquareMatrix::ZeroLowerTriangle(ColumnVector &vector, int &counterOfOperations) {
    for (int i = 0; i < this->getColumns() - 1; ++i) {
        double maximum = abs(this->getAtPosition(i, i));
        int maxIndex = i;
        //Finding max pivot
        for (int j = i + 1; j < this->getRows(); ++j) {
            if (abs(this->getAtPosition(j, i)) > maximum) {
                maximum = abs(this->getAtPosition(j, i));
                maxIndex = j;
            }
        }
        //Permutation
        if (maxIndex != i) {
            PermutationMatrix P(this->getRows(), i, maxIndex, "E");
            *this = P * *this;
            vector.swapRows(i, maxIndex);
        }
        //Elimination
        for (int j = i + 1; j < this->getRows(); ++j) {
            if (this->getAtPosition(j, i) != 0) {
                EliminationMatrix E(j, i, *this);
                vector.subtract(E.getAtPosition(j, i), j, i);
                *this = E * *this;
            }
        }
    }
}

void SquareMatrix::ZeroUpperTriangle(IdentityMatrix &augmented, int &counterOfOperations) {
    for (int i = this->getColumns() - 1; i >= 0; --i) {
        for (int j = i - 1; j >= 0; --j) {
            EliminationMatrix E(j, i, *this);
            *this = E * *this;
            augmented = E * augmented;
        }
    }

}

void SquareMatrix::ZeroUpperTriangle(ColumnVector &vector, int &counterOfOperations) {
    for (int i = this->getColumns() - 1; i >= 0; --i) {
        for (int j = i - 1; j >= 0; --j) {
            if (this->getAtPosition(j, i) != 0 && this->getAtPosition(i, i) != 0) {
                EliminationMatrix E(j, i, *this);
                vector.subtract(E.getAtPosition(j, i), j, i);
                *this = E * *this;
            }
        }
    }
}

void SquareMatrix::DiagonalNormalization(IdentityMatrix &augmented) {
    for (int i = 0; i < this->getRows(); ++i) {
        if (this->getAtPosition(i, i) != 0) {
            for (int j = 0; j < augmented.getColumns(); ++j) {
                double newVal = augmented.getAtPosition(i, j) / this->getAtPosition(i, i);
                augmented.putValueAt(newVal, i, j);
            }
            this->putValueAt(1, i, i);
        }
    }
}

void SquareMatrix::DiagonalNormalization(ColumnVector &vector) {
    for (int i = 0; i < this->getRows(); ++i) {
        if (this->getAtPosition(i, i) != 0) {
            vector.division(this->getAtPosition(i, i), i);
            this->putValueAt(1, i, i);
        }
    }
}
