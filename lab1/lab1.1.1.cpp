#include <iostream>
#include<sys/time.h>
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
	float timeuse;
	struct timeval starttime,endtime;
	gettimeofday(&starttime,NULL);
	Trivial_algorithm(n);
	gettimeofday(&endtime,NULL);
	timeuse = 1000000*(endtime.tv_sec - starttime.tv_sec)+(endtime.tv_usec - starttime.tv_usec);
	timeuse/=1000;
	cout << "timeuse:"<<timeuse<<endl;
	cout<<"result:"<<endl;
	for(int i=0;i<n;i++)
	{
		//cout<<sum[i]<<" ";
	}
	return 0;
}
