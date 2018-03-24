#include <iostream>
#include <cmath>

const double eps = 0.0001;

using namespace std;

//функции
double f1(double x, double y){ return 2 * x * x + 3 * y * y - 6 * y - 4; }
double f2(double x, double y){ return x * x - 3 * y * y - 4 * x - 2; }

//частные производные ф-ий
double f1_dx(double x){ return 4 * x; }
double f1_dy(double y){ return 6 * y - 6; }

double f2_dx(double x){ return 4 * x - 4; }
double f2_dy(double y){ return 6 * y; }

double sum(double *x)
{
    double tmp = 0;
    for(int i = 0; i < sizeof(x)/sizeof(x[0]);i++)
        tmp += x[i];
    return tmp;
}

int main()
{
    setlocale(LC_ALL, "Russian");

    double X[2] = {0.5, -0.5};
    double XX[2];
    double W[2][2];
    double WW[2][2];
    double F[2];
    double det;
    int k = 0;
    do{
        k++;
  
        for(int i = 0; i < 2; i++)
            XX[i] = X[i];

        //матрица Якоби
        W[0][0] = f1_dx(X[0]);
        W[0][1] = f1_dy(X[1]);
        W[1][0] = f2_dx(X[0]);
        W[1][1] = f2_dy(X[1]);

        //матрица функций f1 и f2
        F[0] = f1(X[0], X[1]);
        F[1] = f2(X[0], X[1]);
        
        det = W[0][0] * W[1][1] - W[0][1] * W[1][0];//определитель матрицы(Якобиан)
        if(det == 0){
            cout << "система имеет бесконечное множество решений(Якобиан = 0)" << endl;
            return 0;
        }
        
        //поиск обратной матрицы Якоби
        WW[0][0] = W[1][1] / det;
        WW[1][0] = W[1][0] * (-1) / det;
        WW[0][1] = W[0][1] * (-1) / det;
        WW[1][1] = W[0][0] / det;
        
        double dx[2];//dx = W^(-1) * F(X)

        //перемножение матриц(поиск dx)
        for(int j = 0; j < 2; j++){
            dx[j] = 0;    
            for(int i = 0; i < 2; i++){
                dx[j] += WW[j][i] * F[i];
            }
        }    
        for(int i = 0; i < 2; i++)
            X[i] = XX[i] - dx[i];//X^k = X^(k-1) - W^(-1) * F(X)
           

    }while(fabs( sum(X) - sum(XX) ) > eps && k < 100);

    cout << "кол-во итераций: " << k << endl; 
    cout << "корни: " << X[0] << " " << X[1] << endl;
    
    return 0;
}
