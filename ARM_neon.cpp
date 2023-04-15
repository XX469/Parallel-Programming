#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <arm_neon.h>
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
void Gaussian_elimination(int n) {
	float32x4_t vt, va, vaik, vakj, vaij, vx;
	for (int k = 0; k < n; k++) {
		int j;
		vt = vdupq_n_f32(m[k][k]);
		for (j = k + 1; j + 4 <= n; j += 4) {
			va = vld1q_f32(&m[k][j]);
			va = vdivq_f32(va,vt);
			vst1q_f32(&m[k][j], va);
		}
		for (int t = j; t < n; t++) {
			m[k][t] = m[k][t] / m[k][k];
		}
		m[k][k] = 1.0;
		for (int i = k + 1; i < n; i++) {
			vaik = vdupq_n_f32(m[i][k]);
			int j;
			for (j = k + 1; j + 4 <= n; j += 4) {
				vakj = vld1q_f32(&m[k][j]);
				vaij = vld1q_f32(&m[i][j]);
				vx = vmulq_f32(vakj, vaik);
				vaij = vsubq_f32(vaij, vx);
				vst1q_f32(&m[i][j], vaij);
			}
			for (int t = j; t < n; t++) {
				m[i][t] = m[i][t] - m[i][k] * m[k][t];
			}
			m[i][k] = 0;
		}
	}
}
int main()
{
  struct timespec sts,ets;
  int n[7]={10,50,100,200,500,1000,2000};
  for(int i=0;i<7;i++){
        Create_Matrix(n[i]);
  	timespec_get(&sts, TIME_UTC);
        Gaussian_elimination(n[i]);
	  timespec_get(&ets, TIME_UTC);
	  time_t dsec=ets.tv_sec-sts.tv_sec;
	  long dnsec=ets.tv_nsec-sts.tv_nsec;
          if (dnsec<0){
		    dsec--;
		    dnsec+=1000000000ll;
	  }
	cout<<dsec<<"."<<dnsec<<endl;
  }
	return 0;
}

