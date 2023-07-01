#include<iostream>
#include <cstdio>
#include <fstream>
#include<string>
#include<vector>
#include<sstream>
#include<math.h>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include<cstdlib>
#include<xmmintrin.h>
#include<immintrin.h>
#include <mmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include<pthread.h>
#include<semaphore.h>
#include "omp.h"
using namespace std;

int col=1011;//��������
int rowx=539;//��Ԫ������
int rowb=263;//����Ԫ����


int **X;//��Ԫ�Ӿ���
int **B;//����Ԫ�о���
int **Xe;//ʵ�ʲ�������Ԫ�Ӿ���
int **Be;//ʵ�ʲ����ı���Ԫ�о���
int *Xe_index;//��Ԫ����Ԫ��λ��-�Ӿ����Ӧ����

string pathx="D:\\data\\1011_539_263\\X.txt";//��Ԫ���ļ�·��
string pathb="D:\\data\\1011_539_263\\B.txt";//����Ԫ���ļ�·��

void ini_matrix()//�����ʼ��Ϊ0
{
    X=new int*[rowx];
    Xe_index=new int[rowx];
    for(int i=0;i<rowx;i++){
        X[i]=new int[col];
    }
    for (int i = 0; i < rowx; i++) {
        for (int j = 0; j < col; j++) {
            X[i][j] = 0;
        }
    }

    Xe=new int*[rowx];
    for(int i=0;i<rowx;i++){
        Xe[i]=new int[col];
    }
    for (int i = 0; i < rowx; i++) {
        for (int j = 0; j < col; j++) {
            Xe[i][j] = 0;
        }
    }

    B=new int*[rowb];
    for(int i=0;i<rowb;i++){
        B[i]=new int[col];
    }
    for (int i = 0; i < rowb; i++) {
        for (int j = 0; j < col; j++) {
            B[i][j] = 0;
        }
    }

    Be=new int*[rowb];
    for(int i=0;i<rowb;i++){
        Be[i]=new int[col];
    }
    for (int i = 0; i < rowb; i++) {
        for (int j = 0; j < col; j++) {
            Be[i][j] = 0;
        }
    }
}

void string_to_num(string s,int r,int c,int **a)//���ַ���ת��Ϊ���ֲ��������ж�Ӧλ����1
{
    string temp;
    int num;
    stringstream ss(s);//�ַ�����
    while(ss>>temp){
        stringstream t;
        t<<temp;
        t>>num;
        a[r][c-num-1]=1;
    }
}

void readfile(){//��ȡ�ļ�
    ifstream x;//��Ԫ��
    ifstream b;//����Ԫ��
    x.open(pathx,ios::in);
    if(!x.is_open()){
        cout<<"��Ԫ���ļ���ȡʧ��"<<endl;
        return ;
    }
    vector<string> xiaoyuanzi;
    string temp;
    while (getline(x, temp))
        xiaoyuanzi.push_back(temp);
    x.close();
    for (int i = 0; i < xiaoyuanzi.size(); i++)
        string_to_num(xiaoyuanzi[i], i, col, X);

    b.open(pathb,ios::in);
    if(!b.is_open()){
        cout<<"����Ԫ���ļ���ȡʧ��"<<endl;
        return ;
    }
    vector<string> bei;
    string temp2;
    while (getline(b, temp2))
        bei.push_back(temp2);
    b.close();
    for (int i = 0; i < bei.size(); i++)
        string_to_num(bei[i], i, col, B);
}

void restart()//��������Ҫִ�еľ���ָ��ɳ�ʼ״̬
{
    for (int i = 0; i < rowx; i++) {
        for (int j = 0; j < col; j++) {
            //cout<<Xe[i][j]<<" ";
            Xe[i][j] = X[i][j];
        }
        //cout<<endl;
    }
    for (int i = 0; i < rowb; i++) {
        for (int j = 0; j < col; j++) {
            Be[i][j] = B[i][j];
        }
    }
}

int find_first(int *a,int l,int start=0)//�ҵ���һ��1���ڵ�λ��
{
    for(int i=start;i<l;i++){
        if(a[i]==1)
            return i;
    }
    return -1;
}

int b_is_in_x(int *a,int l)//�ҵ���Ԫ���ж�Ӧ����
{
    for(int i=0;i<rowb;i++){
        if(find_first(Xe[i],l)==find_first(a,l)){
            return i;
        }
    }
    return -1;
}
//===========================�����㷨==============================
void common()//�����㷨
{
    int xnum=rowx-rowb;
    int count_xiao=0;
    int count_zhuan=0;
    long long head,tail,freq;
    long double sumtime=0;
    long double sumtime2=0;
    for(int i=0;i<rowb;i++){
        while(find_first(Be[i],col)!=-1){
            QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
            QueryPerformanceCounter((LARGE_INTEGER*)&head);
            int t=b_is_in_x(Be[i],col);//�ҵ���Ԫ�ӵ�Ŀ����
            QueryPerformanceCounter((LARGE_INTEGER*)&tail);
            sumtime2+=(tail - head) * 1000.0 / freq;
            //cout<<t<<" ";
            if(t!=-1){
                count_xiao++;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(int j=0;j<col;j++){
                    Be[i][j]=Be[i][j]^Xe[t][j];//������
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime+=(tail - head) * 1000.0 / freq;
            }
            else{
                count_zhuan++;
                for(int j=0;j<col;j++){
                    Xe[xnum][j]=Be[i][j];//ת������Ԫ��
                }
                xnum++;
                break;
            }
        }
    }
    cout<<"��Ԫ����:"<<count_xiao<<endl;
    cout<<"��Ԫ��ת��������"<<count_zhuan<<endl;
    cout<<"��Ԫʱ�䣺"<<sumtime<<"ms"<<endl;
    cout<<"����Ԫ��ʱ�䣺"<<sumtime2<<"ms"<<endl;
}
//=============================1��SSE==================================
void elimination1()
{
    int xnum=rowx-rowb;
    int count_xiao=0;
    int count_zhuan=0;
    long long head,tail,freq;
    long double sumtime=0;
    for(int i=0;i<rowb;i++){
        while(find_first(Be[i],col)!=-1){
            int t=b_is_in_x(Be[i],col);//�ҵ���Ԫ�ӵ�Ŀ����
            if(t!=-1){
                int j=0;
                count_xiao++;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(j=0;j+4<=col;j+=4){//������
                    __m128i v1 = _mm_loadu_si128((__m128i*)(Be[i]+j));
                    __m128i v2 = _mm_loadu_si128((__m128i*)(Xe[t]+j));
                    v1 = _mm_xor_si128(v1, v2);
                    _mm_storeu_si128((__m128i*)(Be[i]+j), v1);
                }
                for(;j<col;j++){
                    Be[i][j]=Be[i][j]^Xe[t][j];
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime+=(tail - head) * 1000.0 / freq;
            }
            else{
                count_zhuan++;
                int j=0;
                for(j=0;j+4<=col;j+=4){
                    __m128i v1 = _mm_loadu_si128((__m128i*)(Be[i]+j));
                    _mm_storeu_si128((__m128i*)(Xe[xnum]+j), v1);
                }
                for(;j<col;j++){
                    Xe[xnum][j]=Be[i][j];
                }
                xnum++;
                break;
            }
        }
    }
    cout<<"��Ԫ������"<<count_xiao<<endl;
    cout<<"��Ԫ��ת��������"<<count_zhuan<<endl;
    cout<<"��Ԫʱ�䣺"<<sumtime<<"ms"<<endl;
}
//=====================================================================

//=============================2��AVX==================================
void elimination2()
{
    int xnum=rowx-rowb;
    long long head,tail,freq;
    long double sumtime=0;
    for(int i=0;i<rowb;i++){
        while(find_first(Be[i],col)!=-1){
            int t=b_is_in_x(Be[i],col);//�ҵ���Ԫ�ӵ�Ŀ����
            if(t!=-1){
                int j=0;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(j=0;j+8<=col;j+=8){//������
                    __m256i v1 = _mm256_loadu_si256((__m256i*)(Be[i]+j));
                    __m256i v2 = _mm256_loadu_si256((__m256i*)(Xe[t]+j));
                    v1 = _mm256_xor_si256(v1, v2);
                    _mm256_storeu_si256((__m256i*)(Be[i]+j), v1);
                }
                for(;j<col;j++){
                    Be[i][j]=Be[i][j]^Xe[t][j];
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime+=(tail - head) * 1000.0 / freq;
            }
            else{
                for(int j=0;j<col;j++){
                    Xe[xnum][j]=Be[i][j];//ת������Ԫ��
                }
                xnum++;
                break;
            }
        }
    }
     cout<<"��Ԫʱ�䣺"<<sumtime<<"ms"<<endl;
}
//=====================================================================

//===========================3:����Ԫ���Ż�==============================
void ini_xeindex()
{
    for(int i=0;i<rowx;i++){
        Xe_index[i]=-1;
    }
    for(int i=0;i<rowx;i++){
        for(int j=0;j<col;j++){
            if(Xe[i][j]==1){
                Xe_index[j]=i;//��Ԫ��Ϊj����Ԫ������Ԫ�ӵ�i��
                break;
            }
        }
    }
}

void print_index(){
    for(int i=0;i<rowx;i++){
        cout<<i<<":"<<Xe_index[i]<<" ";
    }
    cout<<endl;
}

int find_first2(int *a,int l,int &first,int start=0)//�ҵ���һ��1���ڵ�λ��
{
    for(int i=start;i<l;i++){
        if(a[i]==1){
            first=i;
            return i;
        }
    }
    return -1;
}

void elimination3()
{
    int xnum=rowx-rowb;
    int count_xiao=0;
    int count_zhuan=0;
    long long head,tail,freq;
    long double sumtime=0;//��Ԫʱ��
    long double sumtime2=0;//����Ԫ��ʱ��
    long double sumtime3=0;//ת������Ԫ��ʱ��
    int first;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    ini_xeindex();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    sumtime2+=(tail - head) * 1000.0 / freq;
    //print_index();
    for(int i=0;i<rowb;i++){
        while(find_first2(Be[i],col,first)!=-1){
            int t=Xe_index[first];
            //cout<<t<<" ";
            if(t!=-1){
                count_xiao++;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(int j=first;j<col;j++){
                    Be[i][j]=Be[i][j]^Xe[t][j];//������
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime+=(tail - head) * 1000.0 / freq;
            }
            else{
                count_zhuan++;
                for(int j=first;j<col;j++){
                    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                    QueryPerformanceCounter((LARGE_INTEGER*)&head);
                    Xe[xnum][j]=Be[i][j];//ת������Ԫ��
                    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                    sumtime3+=(tail - head) * 1000.0 / freq;
                }
                Xe_index[first]=xnum;
                xnum++;
                break;
            }
        }
    }
    cout<<"��Ԫ������"<<count_xiao<<endl;
    cout<<"��Ԫ��ת��������"<<count_zhuan<<endl;
    cout<<"��Ԫʱ�䣺"<<sumtime<<"ms"<<endl;
    cout<<"����Ԫ��ʱ�䣺"<<sumtime2<<"ms"<<endl;
    cout<<"����Ԫ��ת������Ԫ��ʱ�䣺"<<sumtime3<<"ms"<<endl;
}
//===============================================================

//=======================4:AVXȫ�Ż�==============================
void elimination4()
{
    int xnum=rowx-rowb;
    int count_xiao=0;
    int count_zhuan=0;
    long long head,tail,freq;
    long double sumtime=0;//��Ԫʱ��
    long double sumtime2=0;//����Ԫ��ʱ��
    long double sumtime3=0;//ת������Ԫ��ʱ��
    int first;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    ini_xeindex();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    sumtime2+=(tail - head) * 1000.0 / freq;
    //print_index();
    for(int i=0;i<rowb;i++){
        while(find_first2(Be[i],col,first)!=-1){
            int t=Xe_index[first];
            //cout<<t<<" ";
            if(t!=-1){
                count_xiao++;
                int j=0;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(j=0;j+8<=col;j+=8){//������
                    __m256i v1=_mm256_loadu_si256((__m256i*)(Be[i]+j));
                    __m256i v2=_mm256_loadu_si256((__m256i*)(Xe[t]+j));
                    v1 = _mm256_xor_si256(v1, v2);
                    _mm256_storeu_si256((__m256i*)(Be[i]+j), v1);
                }
                for(;j<col;j++){
                    Be[i][j]=Be[i][j]^Xe[t][j];
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime+=(tail - head) * 1000.0 / freq;
            }
            else{
                count_zhuan++;
                int j=0;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(j=0;j+8<=col;j+=8){
                    __m256i v1 = _mm256_loadu_si256((__m256i*)(Be[i]+j));
                    _mm256_storeu_si256((__m256i*)(Xe[xnum]+j), v1);
                }
                for(;j<col;j++){
                    Xe[xnum][j]=Be[i][j];
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime3+=(tail - head) * 1000.0 / freq;
                Xe_index[first]=xnum;
                xnum++;
                break;
            }
        }
    }
    cout<<"��Ԫ������"<<count_xiao<<endl;
    cout<<"��Ԫ��ת��������"<<count_zhuan<<endl;
    cout<<"��Ԫʱ�䣺"<<sumtime<<"ms"<<endl;
    cout<<"����Ԫ��ʱ�䣺"<<sumtime2<<"ms"<<endl;
    cout<<"����Ԫ��ת������Ԫ��ʱ�䣺"<<sumtime3<<"ms"<<endl;
}
//=================================================================

//=======================5:SSEȫ�Ż�==============================
void elimination5()
{
    int xnum=rowx-rowb;
    int count_xiao=0;
    int count_zhuan=0;
    long long head,tail,freq;
    long double sumtime=0;//��Ԫʱ��
    long double sumtime2=0;//����Ԫ��ʱ��
    long double sumtime3=0;//ת������Ԫ��ʱ��
    int first;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    ini_xeindex();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    sumtime2+=(tail - head) * 1000.0 / freq;
    //print_index();
    for(int i=0;i<rowb;i++){
        while(find_first2(Be[i],col,first)!=-1){
            int t=Xe_index[first];
            //cout<<t<<" ";
            if(t!=-1){
                count_xiao++;
                int j=0;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(j=0;j+4<=col;j+=4){//������
                    __m128i v1=_mm_loadu_si128((__m128i*)(Be[i]+j));
                    __m128i v2=_mm_loadu_si128((__m128i*)(Xe[t]+j));
                    v1 = _mm_xor_si128(v1, v2);
                    _mm_storeu_si128((__m128i*)(Be[i]+j), v1);
                }
                for(;j<col;j++){
                    Be[i][j]=Be[i][j]^Xe[t][j];
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime+=(tail - head) * 1000.0 / freq;
            }
            else{
                count_zhuan++;
                int j=0;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(j=0;j+4<=col;j+=4){
                    __m128i v1 = _mm_loadu_si128((__m128i*)(Be[i]+j));
                    _mm_storeu_si128((__m128i*)(Xe[xnum]+j), v1);
                }
                for(;j<col;j++){
                    Xe[xnum][j]=Be[i][j];
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime3+=(tail - head) * 1000.0 / freq;
                Xe_index[first]=xnum;
                xnum++;
                break;
            }
        }
    }
    cout<<"��Ԫ������"<<count_xiao<<endl;
    cout<<"��Ԫ��ת��������"<<count_zhuan<<endl;
    cout<<"��Ԫʱ�䣺"<<sumtime<<"ms"<<endl;
    cout<<"����Ԫ��ʱ�䣺"<<sumtime2<<"ms"<<endl;
    cout<<"����Ԫ��ת������Ԫ��ʱ�䣺"<<sumtime3<<"ms"<<endl;
}
//=================================================================

//=============================6:pthread=======================================
const int NUM_THREADS = 4;
pthread_barrier_t barrier_getBe;//��ȡ���ֱ���Ԫ�м���
pthread_barrier_t barrier_Elimination;//���߳���Ԫ
pthread_barrier_t barrier_Transform;//ת������Ԫ��

int *Be_status;//ÿ������Ԫ�е�ǰ״̬��0��Ϊ�� 1����Ҫ��Ԫ -1����Ҫת������Ԫ�� -2:ת����ɣ�
int **Be_Elimination;//��ǰ��Ҫ������Ԫ�ı���Ԫ�м�����Ԫ�ؼ���
int **Be_transform;//��ǰ��Ҫת������Ԫ�ӵı���Ԫ�м�����Ԫ�ؼ���

int finished_Elimination;//��ǰ��Ԫ��ɵı���Ԫ����
int need_Elimination;//������Ҫ��Ԫ�ı���Ԫ����
int need_transform;//������Ҫת���ı���Ԫ����
int xnum=rowx-rowb;//��Ԫ������

typedef struct {
	int t_id;
}threadParm_t;

void initiate()//��ʼ��
{
    finished_Elimination=0;
    need_Elimination=0;
    need_transform=0;

    Be_status=new int[rowb];
    for(int i=0;i<rowb;i++){
        Be_status[i]=1;
    }

    Be_Elimination=new int*[rowb];
    for(int i=0;i<rowb;i++){
        Be_Elimination[i]=new int[2];
    }

    Be_transform=new int*[rowb];
    for(int i=0;i<rowb;i++){
        Be_transform[i]=new int[2];
    }
}

void clean() //deleteָ��
{
    delete[] Be_status;

    for (int i = 0; i < rowb; i++) {
        delete[] Be_Elimination[i];
    }
    delete[] Be_Elimination;

    for (int i = 0; i < rowb; i++) {
        delete[] Be_transform[i];
    }
    delete[] Be_transform;
}


void *threadFunc(void *param)
{
    threadParm_t *p = (threadParm_t*)param;
    int t_id=p->t_id;
    //cout<<t_id<<" ";
    while(finished_Elimination<rowb){
        //cout<<finished_Elimination<<" ";
        need_Elimination=0;
        need_transform=0;
        //0�Ž��̻�ȡ������Ԫ�ı���Ԫ���к�
        if(t_id==0){
            for(int i=0;i<rowb;i++){
                if(Be_status[i]==1){
                    Be_Elimination[need_Elimination][0]=i;
                    need_Elimination++;
                }
                else if(Be_status[i]==-1){
                    Be_transform[need_transform][0]=i;
                    Be_transform[need_transform][1]=find_first(Be[i],col);
                    need_transform++;
                    Be_status[i]=-2;
                }
            }
        }
        //ͬ��
        pthread_barrier_wait(&barrier_getBe);
        //ÿ���̷߳��������н�����Ԫ������״̬����
        for(int i=t_id;i<need_Elimination;i+=NUM_THREADS){
            int curr=Be_Elimination[i][0];//��ǰ��������Ԫ���к�
            Be_Elimination[i][1]=find_first(Be[curr],col);
            int first=Be_Elimination[i][1];//��Ԫ��
            int t=Xe_index[first];
            if(t!=-1){//��Ԫ�����и��У�ֱ����Ԫ
                for(int j=first;j<col;j++){
                    Be[curr][j]=Be[curr][j]^Xe[t][j];//������
                }
                first=find_first(Be[curr],col,first);
                //cout<<curr<<":"<<first<<" ";
                if(first==-1){//��Ԫ����б�ɿ�
                    Be_status[curr]=0;
                    //cout<<"2 ";
                    finished_Elimination++;
                }
            }
            else{//��Ԫ����û�и��У���������һ��ת������Ԫ��
                Be_status[curr]=-1;
                //cout<<"1 ";
                finished_Elimination++;
            }
        }
        //ͬ��
        pthread_barrier_wait(&barrier_Elimination);
        //0���߳̽���Ԫ��ת��
        if(t_id==0){
            for(int i=0;i<need_transform;i++){
                int curr=Be_transform[i][0];//��ǰ�����к�
                int first=Be_transform[i][1];//��Ԫ��
                if(Xe_index[first]==-1){
                    for(int j=0;j<col;j++){
                        Xe[xnum][j]=Be[curr][j];
                    }
                    Xe_index[first]=xnum;
                    xnum++;
                }
                else{//�����Ѿ��ж�Ӧ��Ԫ�ص���ת������Ԫ�ӣ���һ��������Ԫ
                    Be_status[curr]=1;
                    finished_Elimination--;
                }
            }
        }
        pthread_barrier_wait(&barrier_Transform);
    }
    pthread_exit(NULL);
}

void elimination6()
{
    int xnum=rowx-rowb;
    initiate();
    int first;
    ini_xeindex();
    //��ʼ��barrier
    pthread_barrier_init(&barrier_getBe,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Transform,NULL,NUM_THREADS);
    //�����߳�
    pthread_t handles[NUM_THREADS];
    threadParm_t param[NUM_THREADS];
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        param[t_id].t_id=t_id;
        pthread_create(&handles[t_id],NULL,threadFunc,(void*)&param[t_id]);
    }
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        pthread_join(handles[t_id],0);
    }
    pthread_barrier_destroy(&barrier_getBe);
    pthread_barrier_destroy(&barrier_Elimination);
    pthread_barrier_destroy(&barrier_Transform);
    clean();
}
//=====================================================================

//============================7:pthread+AVX=============================
void *threadFunc2(void *param)
{
    threadParm_t *p = (threadParm_t*)param;
    int t_id=p->t_id;
    //cout<<t_id<<" ";
    while(finished_Elimination<rowb){
        //cout<<finished_Elimination<<" ";
        need_Elimination=0;
        need_transform=0;
        //0�Ž��̻�ȡ������Ԫ�ı���Ԫ���к�
        if(t_id==0){
            for(int i=0;i<rowb;i++){
                if(Be_status[i]==1){
                    Be_Elimination[need_Elimination][0]=i;
                    need_Elimination++;
                }
                else if(Be_status[i]==-1){
                    Be_transform[need_transform][0]=i;
                    Be_transform[need_transform][1]=find_first(Be[i],col);
                    need_transform++;
                    Be_status[i]=-2;
                }
            }
        }
        //ͬ��
        pthread_barrier_wait(&barrier_getBe);
        //ÿ���̷߳��������н�����Ԫ������״̬����
        for(int i=t_id;i<need_Elimination;i+=NUM_THREADS){
            int curr=Be_Elimination[i][0];//��ǰ��������Ԫ���к�
            Be_Elimination[i][1]=find_first(Be[curr],col);
            int first=Be_Elimination[i][1];//��Ԫ��
            int t=Xe_index[first];
            if(t!=-1){//��Ԫ�����и��У�ֱ����Ԫ
                int j;
                for(j=first;j+8<=col;j+=8){//������
                    __m256i v1=_mm256_loadu_si256((__m256i*)(Be[curr]+j));
                    __m256i v2=_mm256_loadu_si256((__m256i*)(Xe[t]+j));
                    v1 = _mm256_xor_si256(v1, v2);
                    _mm256_storeu_si256((__m256i*)(Be[curr]+j), v1);
                }
                for(;j<col;j++){
                    Be[curr][j]=Be[curr][j]^Xe[t][j];
                }
                first=find_first(Be[curr],col,first);
                //cout<<curr<<":"<<first<<" ";
                if(first==-1){//��Ԫ����б�ɿ�
                    Be_status[curr]=0;
                    //cout<<"2 ";
                    finished_Elimination++;
                }
            }
            else{//��Ԫ����û�и��У���������һ��ת������Ԫ��
                Be_status[curr]=-1;
                //cout<<"1 ";
                finished_Elimination++;
            }
        }
        //ͬ��
        pthread_barrier_wait(&barrier_Elimination);
        //0���߳̽���Ԫ��ת��
        if(t_id==0){
            for (int i = 0; i < need_transform; i++) {
                int curr = Be_transform[i][0];//��ǰ�����к�
                int first = Be_transform[i][1];//��Ԫ��
                if (Xe_index[first] == -1) {
                    int j;
                    for (j = first; j + 8 <= col; j += 8) {
                        __m256i v1 = _mm256_loadu_si256((__m256i*)(Be[curr] + j));
                        _mm256_storeu_si256((__m256i*)(Xe[xnum] + j), v1);
                    }
                    for (; j < col; j++) {
                        Xe[xnum][j] = Be[curr][j];
                    }
                    Xe_index[first] = xnum;
                    xnum++;
                }
                else {//�����Ѿ��ж�Ӧ��Ԫ�ص���ת������Ԫ�ӣ���һ��������Ԫ
                    Be_status[curr] = 1;
                    finished_Elimination--;
                }
            }
        }
        pthread_barrier_wait(&barrier_Transform);
    }
    pthread_exit(NULL);
}

void elimination7()
{
    int xnum=rowx-rowb;
    initiate();
    int first;
    ini_xeindex();
    //��ʼ��barrier
    pthread_barrier_init(&barrier_getBe,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Transform,NULL,NUM_THREADS);
    //�����߳�
    pthread_t handles[NUM_THREADS];
    threadParm_t param[NUM_THREADS];
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        param[t_id].t_id=t_id;
        pthread_create(&handles[t_id],NULL,threadFunc2,(void*)&param[t_id]);
    }
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        pthread_join(handles[t_id],0);
    }
    pthread_barrier_destroy(&barrier_getBe);
    pthread_barrier_destroy(&barrier_Elimination);
    pthread_barrier_destroy(&barrier_Transform);
    clean();
}
//======================================================================

//=======================8:OpenMP+AVX=======================================
void elimination8()
{
    int xnum=rowx-rowb;
    initiate();
    int first;
    ini_xeindex();
    #pragma omp parallel num_threads(NUM_THREADS)
    while(finished_Elimination<rowb){
        cout<<finished_Elimination<<" ";
        need_Elimination=0;
        need_transform=0;
        #pragma omp single
        for(int i=0;i<rowb;i++){
            if(Be_status[i]==1){
                Be_Elimination[need_Elimination][0]=i;
                need_Elimination++;
            }
            else if(Be_status[i]==-1){
                Be_transform[need_transform][0]=i;
                Be_transform[need_transform][1]=find_first(Be[i],col);
                need_transform++;
                Be_status[i]=-2;
            }
        }
        #pragma omp for
        for(int i=0;i<need_Elimination;i++){
            int curr=Be_Elimination[i][0];//��ǰ��������Ԫ���к�
            Be_Elimination[i][1]=find_first(Be[curr],col);
            int first=Be_Elimination[i][1];//��Ԫ��
            int t=Xe_index[first];
            if(t!=-1){//��Ԫ�����и��У�ֱ����Ԫ
                int j;
                for(j=first;j+8<=col;j+=8){//������
                    __m256i v1=_mm256_loadu_si256((__m256i*)(Be[curr]+j));
                    __m256i v2=_mm256_loadu_si256((__m256i*)(Xe[t]+j));
                    v1 = _mm256_xor_si256(v1, v2);
                    _mm256_storeu_si256((__m256i*)(Be[curr]+j), v1);
                }
                for(;j<col;j++){
                    Be[curr][j]=Be[curr][j]^Xe[t][j];
                }
                first=find_first(Be[curr],col,first);
                //cout<<curr<<":"<<first<<" ";
                if(first==-1){//��Ԫ����б�ɿ�
                    Be_status[curr]=0;
                    //cout<<"2 ";
                    finished_Elimination++;
                }
            }
            else{//��Ԫ����û�и��У���������һ��ת������Ԫ��
                Be_status[curr]=-1;
                //cout<<"1 ";
                finished_Elimination++;
            }
        }
        #pragma omp single
        for (int i = 0; i < need_transform; i++) {
            int curr = Be_transform[i][0];//��ǰ�����к�
            int first = Be_transform[i][1];//��Ԫ��
            if (Xe_index[first] == -1) {
                int j;
                for (j = 0; j + 8 <= col; j += 8) {
                    __m256i v1 = _mm256_loadu_si256((__m256i*)(Be[curr] + j));
                    _mm256_storeu_si256((__m256i*)(Xe[xnum] + j), v1);
                }
                for (; j < col; j++) {
                    Xe[xnum][j] = Be[curr][j];
                }
                Xe_index[first] = xnum;
                xnum++;
            }
            else {//�����Ѿ��ж�Ӧ��Ԫ�ص���ת������Ԫ�ӣ���һ��������Ԫ
                Be_status[curr] = 1;
                finished_Elimination--;
            }
        }
    }
    clean();
}
//======================================================================

//=======================9:OpenMP=======================================
void elimination9()
{
    int xnum=rowx-rowb;
    initiate();
    int first;
    ini_xeindex();
    #pragma omp parallel num_threads(NUM_THREADS)
    while(finished_Elimination<rowb){
        cout<<finished_Elimination<<" ";
        need_Elimination=0;
        need_transform=0;
        #pragma omp single
        for(int i=0;i<rowb;i++){
            if(Be_status[i]==1){
                Be_Elimination[need_Elimination][0]=i;
                need_Elimination++;
            }
            else if(Be_status[i]==-1){
                Be_transform[need_transform][0]=i;
                Be_transform[need_transform][1]=find_first(Be[i],col);
                need_transform++;
                Be_status[i]=-2;
            }
        }
        #pragma omp for
        for(int i=0;i<need_Elimination;i++){
            int curr=Be_Elimination[i][0];//��ǰ��������Ԫ���к�
            Be_Elimination[i][1]=find_first(Be[curr],col);
            int first=Be_Elimination[i][1];//��Ԫ��
            int t=Xe_index[first];
            if(t!=-1){//��Ԫ�����и��У�ֱ����Ԫ
                for(int j=first;j<col;j++){
                    Be[curr][j]=Be[curr][j]^Xe[t][j];//������
                }
                first=find_first(Be[curr],col,first);
                //cout<<curr<<":"<<first<<" ";
                if(first==-1){//��Ԫ����б�ɿ�
                    Be_status[curr]=0;
                    //cout<<"2 ";
                    finished_Elimination++;
                }
            }
            else{//��Ԫ����û�и��У���������һ��ת������Ԫ��
                Be_status[curr]=-1;
                //cout<<"1 ";
                finished_Elimination++;
            }
        }
        #pragma omp single
        for(int i=0;i<need_transform;i++){
            int curr=Be_transform[i][0];//��ǰ�����к�
            int first=Be_transform[i][1];//��Ԫ��
            for(int j=0;j<col;j++){
                Xe[xnum][j]=Be[curr][j];
            }
            Xe_index[first]=xnum;
            xnum++;
        }
    }
    clean();
}
//======================================================================
int main(){
    rowx+=rowb;//��Ԫ��������
    ini_matrix();
    readfile();
    //=========================ƽ���㷨=====================================
    cout<<"common==============================================="<<endl;
    restart();
    long long head,tail,freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    common();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================1:SSE=====================================
    cout<<"1:SSE==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination1();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================2:AVX=====================================
    cout<<"2:AVX==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination2();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================3:����Ԫ���Ż�=====================================
    cout<<"3:����Ԫ���Ż�==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination3();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================4:AVXȫ�Ż�=====================================
    cout<<"4:AVXȫ�Ż�==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination4();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================5:SSEȫ�Ż�=====================================
    cout<<"5:SSEȫ�Ż�==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination5();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================6:pthread=====================================
    cout<<"6:pthread==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination6();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================7:pthread+AVX=====================================
    cout<<"7:pthread+AVX==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination7();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================8:openmp+AVX=====================================
    cout<<"8:openmp+AVX==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination8();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================9:openmp=====================================
    cout<<"9:openmp==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination9();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================
    return 0;
}
