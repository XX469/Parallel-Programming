#include <iostream>
#include <windows.h>
#include <stdlib.h>
#include <xmmintrin.h> //SSE
#include <emmintrin.h> //SSE2
#include <pmmintrin.h> //SSE3
#include <tmmintrin.h> //SSSE3
#include <smmintrin.h> //SSE4.1
#include <nmmintrin.h> //SSSE4.2
#include <immintrin.h> //AVX¡¢AVX2
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
    __m256 vt,va;
    __m256 vaik,vakj,vaij,vx;
    for(int k=0;k<n;k++){
        int j;
        vt = _mm256_set1_ps(m[k][k]);//dupTo4Float(m[k,k])
        for(j=k+1;j+8<=n;j+=8){
            va = _mm256_loadu_ps(&m[k][j]);//load4FloatFrom(&m[k,j])
            va = _mm256_div_ps(va,vt);//va/vt
            _mm256_storeu_ps(&m[k][j],va);//store4FloatTo(&m[k,j],va)
        }
        for(int t=j ;t<n;t++){//The remaining part
            m[k][t]=m[k][t]/m[k][k];
        }
        m[k][k]=1.0;
        for(int i=k+1;i<n;i++){
            vaik=_mm256_set1_ps(m[i][k]);//dupToVector4(m[i,k])
            int j;
            for(j=k+1;j+8<=n;j+=8){
                vakj= _mm256_loadu_ps(&m[k][j]);//load4FloatFrom(&m[k,j])
                vaij= _mm256_loadu_ps(&m[i][j]);//load4FloatFrom(&m[i,j])
                vx=_mm256_mul_ps(vakj,vaik);//mul vakj*vaik
                vaij=_mm256_sub_ps(vaij,vx);//sub vaij-vx
                _mm256_storeu_ps(&m[i][j],vaij);//store4FloatTo(&m[i,j],vaij)
            }
            for(int t=j;t<n;t++){//The remaining part
                m[i][t]=m[i][t]-m[i][k]*m[k][t];
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
