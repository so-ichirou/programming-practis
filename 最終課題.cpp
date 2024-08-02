#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <chrono>
#include <string>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>

#define GNUPLOT_PATH "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\""/*Gnuplotを用いるため、使うときはここのexeを変えること*/


using namespace std;
using cplx = complex<double>;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Wing 
{
public:

    int point;
    float y_max;
    float x_max;
    float a;
    string path;
    int nodes;
    vector<cplx> z_foil;
    vector<double> l;
    vector<vector<cplx>> S1;
    vector<vector<cplx>> S2;
    vector<vector<double>> A;
    vector<double> B;
    float B_;


    Wing(int temp_point, float temp_y_max, float temp_x_max, float temp_a, float temp_B, string temp_name) 
    {
        point = temp_point;
        y_max = temp_y_max;
        x_max = temp_x_max;
        a = temp_a;
        B_ = temp_B;
        path = temp_name + ".dat";
        generate_wings();
        nodes = Check_FileSize(path) - 1;
        z_foil = vector<cplx>(nodes + 1);
    }

    /*翼型を生成してファイルに格納する関数*/
    void generate_wings()
    {
        float b, c, d;
        float A, B, C, D;

        // y1 = a * x1.^3 + b * x1.^2 + c * x1 + d //
        b = -2 * a * x_max - y_max / pow(x_max, 2);
        c = -3 * a * pow(x_max, 2) - 2 * b * x_max;
        d = 0;

        // y2 = A * x2.^3 + B_ * x2.^2 + C * x2 + D //
       
        A = (-B_ * (pow(x_max, 2) - 2 * x_max + 1) - y_max) / (2 * pow(x_max, 3) - 3 * pow(x_max, 2) + 1);
        C = -3 * A * pow(x_max, 2) - 2 * B_ * x_max;
        D = -A - B_ - C;

        float dx = (float)1 / point;

        ofstream outfile(path, ios_base::out);
        vector<float> x, y;
        vector<float>::iterator it = x.begin();         //iterator
        for (int i = 0; i <= point; i++)
        {
            x.push_back(1 - dx * i);
        }
        for (auto it = x.begin(); it != x.end(); ++it)
        {
            if (*it < x_max)

            {
                y.push_back(a * pow(*it, 3) + b * pow(*it, 2) + c * *it + d);
            }
            else
            {
                y.push_back(A * pow(*it, 3) + B_ * pow(*it, 2) + C * *it + D);
            }

        }
        for (int i = 0; i <= point; i++)
        {

            outfile << x[i] << " " << y[i] << endl;

        }
        for (int i = point - 1; i >= 0; i--)
        {

            outfile << x[i] << " " << -y[i] << endl;

        }
        outfile.close();
    }

    /*一行目から最終行まで読み込んで行数を数える関数*/
    int Check_FileSize(const string& file_name)
    {
        ifstream file(file_name);
        if (!file.is_open())
        {
            cerr << "Not Found " << file_name << "," << endl;
            exit(EXIT_FAILURE);
        }
        int nodes = 0;
        string line;
        while (getline(file, line))
        {
            nodes++;
        }
        file.close();
        return nodes;
    }

    /*ファイルを読み込み、複素座標に置き換える*/
    void Input_Airfoil()
    {
        ifstream infile(path);
        if (!infile)
        {
            cerr << "Not Found " << path << "." << std::endl;
            exit(EXIT_FAILURE);
        }

        vector<vector<double>> val;
        string line;
        while (getline(infile, line))
        {
            istringstream iss(line);
            vector<double> row;
            double num;
            while (iss >> num)
            {
                row.push_back(num);
            }
            val.push_back(row);
        }

        infile.close();

        for (int n = 0; n < nodes + 1; ++n)
        {
            double x = val[n][0];
            double y = val[n][1];
            z_foil[n] = cplx(x, y);
        }
    }

    /*S1,S2行列の計算を行う*/
    void Make_S()
    {
        vector<cplx> z_ref(nodes);
        cplx I(0, 1);
        for (int i = 0; i < nodes; i++)
        {
            l[i] = std::abs(z_foil[i + 1] - z_foil[i]);
            z_ref[i] = (z_foil[i] + z_foil[i + 1]) / 2.0;
        }

        for (int k = 0; k < nodes; k++)
        {
            for (int j = 0; j < nodes; j++)
            {
                S1[k][j] = ((I * l[j]) / (2.0 * M_PI * (z_foil[j + 1] - z_foil[j]))) *
                    (-1.0 + (z_foil[j + 1] - z_ref[k]) / (z_foil[j + 1] - z_foil[j]) *
                        std::log((z_foil[j + 1] - z_ref[k]) / (z_foil[j] - z_ref[k])));
                S2[k][j] = ((I * l[j]) / (2.0 * M_PI * (z_foil[j + 1] - z_foil[j]))) *
                    (1.0 - (z_foil[j] - z_ref[k]) / (z_foil[j + 1] - z_foil[j]) *
                        std::log((z_foil[j + 1] - z_ref[k]) / (z_foil[j] - z_ref[k])));
            }
        }
    }

    /*A行列を求める*/
    void Make_A()
    {
        vector<cplx> z_n(nodes);
        cplx I(0, 1);

        // Make z_n
        for (int i = 0; i < nodes; i++) 
        {
            z_n[i] = ((z_foil[i + 1] - z_foil[i]) / l[i]) * (1.0 / I); // eq(4.29)
        }

        // Make A eq(4.40)
        for (int k = 0; k < nodes; k++) 
        {
            for (int j = 0; j <= nodes; j++)
            {
                if (j == 0) 
                {
                    A[k][j] = std::real(S1[k][j] * z_n[k]);
                }
                else if (j == nodes) {
                    A[k][j] = std::real(S2[k][j - 1] * z_n[k]);
                }
                else {
                    A[k][j] = std::real((S2[k][j - 1] + S1[k][j]) * z_n[k]);
                }
            }
        }
        A[nodes][0] = 1.0;
        A[nodes][nodes] = 1.0; // Kutta condition
    }

    /*B行列を求める*/
    void Make_B(cplx U_inf)
    {
        vector<cplx> z_n(nodes);
        cplx I(0, 1);

        // Make z_n
        for (int i = 0; i < nodes; i++) 
        {
            z_n[i] = ((z_foil[i + 1] - z_foil[i]) / l[i]) * (1.0 / I); // eq(4.29)
        }

        // Make B eq(4.40d)
        for (int k = 0; k < nodes; k++)
        {
            B[k] = -std::real(U_inf * z_n[k]);
        }
        B[nodes] = 0; // Kutta Condition
    };

    /*ガウスの消去法*/
    void Gaussian_Elimination(int nodes, vector<vector<double>>& A,
        vector<double>& B)
    {
        double p, pivot;
        int i, j, ii, jj;

        for (i = 0; i <= nodes; i++) 
        {
            pivot = A[i][i];
            for (j = i; j <= nodes; j++)
            {
                A[i][j] = A[i][j] / pivot;
            }
            B[i] = B[i] / pivot;

            for (ii = i + 1; ii <= nodes; ii++)
            {
                p = A[ii][i];
                for (jj = i; jj <= nodes; jj++)
                {
                    A[ii][jj] = A[ii][jj] - p * A[i][jj];
                }
                B[ii] = B[ii] - p * B[i];
            }
        }

        for (i = nodes; i > 0; i--)
        {
            for (ii = i - 1; ii >= 0; ii--)
            {
                p = A[ii][i];
                for (j = i; j <= nodes; j++) 
                {
                    A[ii][j] = A[ii][j] - p * A[i][j];
                }
                B[ii] = B[ii] - p * B[i];
            }
        }
    }

    /*行列計算*/
    void Calculate_Matrix(cplx U_inf)
    {
        l = vector<double>(nodes);

        S1 = vector<vector<cplx>>(nodes, vector<cplx>(nodes));
        S2 = vector<vector<cplx>>(nodes, vector<cplx>(nodes));
        Make_S();
        A = vector<vector<double>>(nodes + 1, vector<double>(nodes + 1));
        Make_A();
        B = vector<double>(nodes + 1);
        Make_B(U_inf);
    }

};




void output_cp(const char* filename, vector<double>& cp)
{
    ofstream outfile(filename, ios::out);

    for (const auto& value : cp)
    {
        outfile << value << endl;
    }
    outfile.close();
}

void Plot_Airfoil(int nodes, vector<cplx>& z_foil, vector<double>& cp)
{
    ofstream pipe("airfoil_plot.plt", ios_base::out);
    pipe << "set terminal pngcairo enhanced\n";
    pipe << "set output 'airfoil.png'\n";
    pipe << "set xrange [-0.1:1.1]\n";
    pipe << "set yrange [-0.3:0.3]\n";
    pipe << "plot '-' using 1:2 title 'Airfoil' with lines\n";

    for (int i = 0; i <= nodes; i++) 
    {
        pipe << real(z_foil[i]) << " " << imag(z_foil[i]) << endl;
    }

    pipe << "e\n";
    pipe << "set terminal pngcairo enhanced\n";
    pipe << "set output 'cp.png'\n";
    pipe << "set xrange [-0.1:1.1]\n";
    pipe << "set yrange [1:-1]\n";
    pipe << "plot '-' using 1:2 title 'Cp' with lines\n";

    for (int i = 0; i <= nodes; i++) 
    {
        pipe << real(z_foil[i]) << " " << cp[i] << endl;
    }

    pipe << "e\n";
    pipe.close();
    system(GNUPLOT_PATH " airfoil_plot.plt");
}

int main(int argc, char* argv[])
{
    srand(time(nullptr));
    Wing wing = Wing(100, 0.1, 0.5, 0.1, 0.3, "wing1");
    wing.generate_wings();
    wing.Input_Airfoil();


    double alpha = 0.0;
    double amin = 0.0;
    double amax = 10.0;
    double da = 2.0;

    cplx I(0, 1);
    cplx U_inf = 1.0 * exp(I * alpha * M_PI / 180.0);
    wing.Calculate_Matrix(U_inf);

    int nodes = wing.nodes;
    vector<cplx> z_foil = wing.z_foil;

    vector<double> cp(nodes + 1);
    for (int i = 0; i <= nodes; i++) 
    {
        cp[i] = 1 - pow(wing.B[i], 2);
    }

    output_cp("cp.dat", cp);
    Plot_Airfoil(nodes, z_foil, cp);

    return 0;
}
