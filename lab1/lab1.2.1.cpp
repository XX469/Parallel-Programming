#include <iostream>
#include<sys/time.h>
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
void Trivial_algorithm(int n)
{
	sum=0;
	for(int i=0;i<n;i++)
	{
		sum+=a[i];
	}
}
int main()
{
	int n;
	cin>>n;
	CreateArray(n);
	float timeuse;
	struct timeval starttime,endtime;
	gettimeofday(&starttime,NULL);
	Trivial_algorithm(n);
	gettimeofday(&endtime,NULL);
	timeuse = 1000000*(endtime.tv_sec - starttime.tv_sec)+(endtime.tv_usec - starttime.tv_usec);
	timeuse /= 1000;
	cout << "timeuse:"<<timeuse<<endl;
	cout<<"result:"<<sum<<endl;
	return 0;
}

