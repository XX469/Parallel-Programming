#include <iostream>
#include <windows.h>
#include <stdlib.h>

using namespace std;

const int N = 9999999;
double a[N];
double sum=0;
void CreateArray(int n)
{
    for(int i=0;i<n;i++)
        {
            a[i]=i;
        }
}
void unroll(int n)
{
    sum=0;
    double sum1=0;
    double sum2=0;
    double sum3=0;
    double sum4=0;
    double sum5=0;
    double sum6=0;
    double sum7=0;
    double sum8=0;
    for(int i=0;i<n;i+=8)
        {
            sum1+=a[i];
            sum2+=a[i+1];
            sum3+=a[i+2];
            sum4+=a[i+3];
            sum5+=a[i+4];
            sum6+=a[i+5];
            sum7+=a[i+6];
            sum8+=a[i+7];
        }
    sum=sum1+sum2+sum3+sum4+sum5+sum6+sum7+sum8;
}
int main()
{
    int n;
    cin>>n;
    CreateArray(n);
    long long head,tail,freq;
    QueryPerformanceFrequency((LARGE_INTEGER *)&freq );
    QueryPerformanceCounter((LARGE_INTEGER *)&head);
    for(int t=0;t<1000;t++)
        {
          unroll(n);
        }
    QueryPerformanceCounter((LARGE_INTEGER *)&tail );
    cout << "Col:" << (tail - head) * 1000.0 / freq << "ms" << endl ;
    cout<<"result:"<<sum<<endl;
    return 0;
}
