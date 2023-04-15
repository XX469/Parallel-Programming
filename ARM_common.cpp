#include <iostream>
#include <time.h>
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
			m[k][j]=m[k][j]/m[k][k];
		}
		m[k][k]=1.0;
		for(int i=k+1;i<n;i++){
			for(int t=k+1;t<n;t++){
				m[i][t]=m[i][t]-m[i][k]*m[k][t];
			}
			m[i][k]=0;
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


