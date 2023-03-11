#include <iostream>
#include <windows.h>
#include <stdlib.h>

using namespace std;

const int N = 10240;

double b[N][N] , a[N],sum[N] ;

void CreateMatrix(int n)
{
    for(int i=0;i<n;i++)
        {
            a[i]=i;
            for(int j=0;j<n;j++)
                {
                    b[i][j]=i+j;
                }
        }
}

void Trivial_algorithm(int n)
{
    for(int i=0;i<n;i++)
        {
            sum[i]=0.0;
            for(int j=0;j<n;j++)
            {
                sum[i]+=b[j][i]*a[j];
            }
        }
}
int main()
{
    int n;
    cin>>n;
    CreateMatrix(n);
    long long head,tail,freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq );
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int t=0;t<2000;t++)
        {
            Trivial_algorithm(n);
        }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail );
    cout << "Col:" << (tail - head) * 1000.0 / freq << "ms" << endl ;
    cout<<"result:"<<endl;
    for(int i=0;i<n;i++)
        {
            cout<<sum[i]<<" ";
        }
    return 0;
}
