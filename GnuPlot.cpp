ColumnVector LeastSquares(int n, int m, Matrix& input) {
    Matrix A(m, n + 1, "A");
    ColumnVector b(m);
    fillMatrixVector(A, n, input, b);

    Matrix A_T = A.T();
    Matrix temp = A_T * A;
    SquareMatrix A_TA = static_cast<SquareMatrix &> (temp);
    A_TA.setName("A_T*A");
    SquareMatrix inverse = A_TA.inverseMatrix();
    inverse.setName("inverse");

    ColumnVector A_T_b = A_T * b;
    ColumnVector ans = inverse * A_T_b;
    return ans;
}


int main() {
    FILE* file = _popen("C:\\gnuplot\\bin\\gnuplot -persist","w");

    if (file != nullptr) {
        int m, n;
        cin >> m;
        Matrix input(m, 2,"input");
        cin >> input;
        cin >> n;

        double minX = 1000000, minY = 1000000, maxY = -1000000, maxX = -1000000;

        ColumnVector ans = LeastSquares(n,m,input);

        for (int i = 0; i < input.getRows(); ++i) {
            minX = min(minX,input.getAtPosition(i,0));
            maxX = max(maxX, input.getAtPosition(i,0));
            minY = min(minY,input.getAtPosition(i,1));
            maxY = max(maxY, input.getAtPosition(i,1));
        }

        fprintf(file, "set xrange [%f:%f]\n", minX-10, maxX+10);
        fprintf(file, "set yrange [%f:%f]\n", minY-10, maxY+10);
        fprintf(file, "set samples %d\n", 1000000);

        string formula;
        formula += to_string(ans.valueAt(0));
        for (int i = 1; i < ans.getSize(); ++i) {
            if (ans.valueAt(i) >= 0) {
                formula += '+';
            }
            formula += to_string(ans.valueAt(i));
            formula += '*';
            for (int j = 0; j < i; ++j) {
                formula += 'x';
                if (j!=i-1) {
                    formula += '*';
                }
            }
        }

        fprintf(file, "plot '-' with points title '', %s with lines title ''\n", formula.c_str());
        for (int i = 0; i < input.getRows(); ++i) {
            fprintf(file, "%lf %lf\n", input.getAtPosition(i,0) , input.getAtPosition(i,1));
        }
        cout << formula;
        fprintf(file, "plot %s with lines\n", formula.c_str());

        pclose(file);
    }
    return 0;
}
