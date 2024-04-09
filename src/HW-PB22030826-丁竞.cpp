#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
const int N = 101; // 矩阵大小
const double epsilon = 1e-4; // 精度要求
 int n;
double er;
double a;
double h;
// 列主元Gauss消元法
void gaussElimination(double A[N][N], double b[N], double x_g[N]) 
{
    // 消元过程
    for (int i = 0; i < N - 1; ++i) {
        // 选取主元
        int pivot = i;
        for (int j = i + 1; j < N; ++j) {
            if (abs(A[j][i]) > abs(A[pivot][i]))
                pivot = j;
        }
        // 交换行
        if (pivot != i) 
        {
            swap(b[i], b[pivot]);
            for (int j = i; j < N; ++j)
                swap(A[i][j], A[pivot][j]);
        }
        // 消元
        for (int j = i + 1; j < N; ++j) {
            double ratio = A[j][i] / A[i][i];
            for (int k = i; k < N; ++k)
                A[j][k] -= ratio * A[i][k];
            b[j] -= ratio * b[i];
        }
    }

    // 回代求解
    for (int i = N - 1; i >= 0; --i) {
        x_g[i] = b[i] / A[i][i];
        for (int j = i + 1; j < N; ++j)
            x_g[i] -= A[i][j] * x_g[j] / A[i][i];
    }
}
// Gauss-Seidel迭代法
void gaussSeidel(double A[N][N],double b[N],double x_gs[N]) 
{
    double tolerance = 1e-6;
    int maxIter  = 100;
    for (int k = 0; k < maxIter; ++k) 
    {
        double maxError = 0.0;
        for (int i = 0; i < N; ++i) 
        {
            double sum1 = 0.0, sum2 = 0.0;
            for (int j = 0; j < i; ++j)
                sum1 += A[i][j] * x_gs[j];
            for (int j = i + 1; j < N; ++j)
                sum2 += A[i][j] * x_gs[j];
            double newX = (b[i] - sum1 - sum2) / A[i][i];
            maxError = max(maxError, abs(newX - x_gs[i]));
            x_gs[i] = newX;
        }
        if (maxError < tolerance)
            break;
    }
}

void initial_A(double (&A)[N][N])
{
     
    for ( auto i = 0; i < N; i++)
     {
        for (auto j = 0; j < N; j++)
        {
            if(i==j)
            A[i][j]=-(2*er+h);
            if(i-j==1)
            A[i][j]=er+h;
            if(i-j==-1)
            A[i][j]=er;
        }
        
     }

}
void initial_b(double (&b)[N])
{
    b[N]={0};
 for (auto i = 0; i < N-1; i++)
      {
        b[i]=a*h*h;
      }
      b[N-1]=a*h*h-er-h;

}

//计算向量范数差值
double compare(double (&a)[N],double (&b)[N])
{
    double sum=0;
  for (auto i = 0; i < N; i++)
  {
    sum+=pow(a[i]-b[i],2);
  }
  return sqrt(sum);
}

double max_error(double (&a)[N],double (&b)[N])
{
    double error=0;
    for (auto i = 0; i <N; i++)
    {
        error=max(error,abs(a[i]-b[i]));
    }
    return error;
}

int main() 
{
   
    cout<<"Please enter the value of n"<<endl;
    cin>>n;
    
    cout<<"PLease enter rhe value of epsilon"<<endl;
    cin>>er;

    cout<<"PLease enter rhe value of a"<<endl;
    cin>>a;

    h=0.01;

    double A[N][N]={0}; 
    double b[N]{0}; 
    double x_g[N]={0};
    double x_gs[N]={0}; 

    initial_A(A);
    initial_b(b);
    
    // 调用列主元Gauss消元法求解线性方程组
    gaussElimination(A, b, x_g);
   
   //再次初始化
    initial_A(A);
    initial_b(b);

    // 调用Gauss-Seidel迭代方法求解线性方程组
    gaussSeidel(A, b, x_gs);
 
    // 输出解向量x
    cout << "Gauss:x:" << endl;
    for (int i = 0; i < N; i++)
    //cout << x_g[i] << endl;
    printf("%.4f  ",x_g[i]);
    printf("\n");
     cout << "Gauss-Seidel:x:" << endl;
    for (int i = 0; i < N; i++)
    //cout <<fixed<<setprecision(4)<<x_gs[i] << endl;
     printf("%.4f  ",x_gs[i]);
     printf("\n");
    double ans[N]={0};
    for (auto i = 0; i < N; i++)
    {
        ans[i]=(1-a)*(1-exp(-h*i/er))/(1-exp(-1/er))+a*i*h;
    }
    cout<<"the right value is:"<<endl;
     for (int i = 0; i < N; i++)
    printf("%.4f  ",ans[i]);
    printf("\n");
    
    cout<<"the error of Gauss is:"<<compare(x_g,ans)<<endl;
    cout<<"the max_error of Gauss is:"<<max_error(x_g,ans)<<endl;
    cout<<"the error of Gauss-Seidel is:"<<compare(x_gs,ans)<<endl;
    cout<<"the max_error of Gauss-Seidel is:"<<max_error(x_gs,ans)<<endl;
    
   system("pause");
   
    return 0;
}
