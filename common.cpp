#include <iostream>
#include <windows.h>
#include <stdlib.h>
using namespace std;
const long long int N = 10000;
float m[N][N];
void Create_Matrix(int n)
{
    for(int i=0;i<n;i++)
        {
            m[i][i]=1.0;
            for(int j=i+1;j<n;j++)
                {
                    m[i][j]=rand();
                }
        }
    for(int k=0;k<n;k++){
            for(int i=k+1;i<n;i++){
                    for(int j=0;j<n;j++){
                        m[i][j]+=m[k][j];
                    }
            }
    }

}
void Gaussian_elimination(int n){
    for(int k=0;k<n;k++){
        for(int j=k+1;j<n;j++){
            m[k][j]=m[k][j]/m[k][k];//division
        }
        m[k][k]=1.0;
        for(int i=k+1;i<n;i++){
            for(int t=k+1;t<n;t++){
                m[i][t]=m[i][t]-m[i][k]*m[k][t];//subtraction
            }
            m[i][k]=0;
        }
    }
}
int main()
{
    int n;
    cin>>n;
    Create_Matrix(n);
    long long head,tail,freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq );
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    Gaussian_elimination(n);
    QueryPerformanceCounter((LARGE_INTEGER *)&tail );
    cout << "Col:" << (tail - head) * 1000.0 / freq << "ms" << endl ;
    return 0;
}
