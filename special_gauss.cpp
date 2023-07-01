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

int col=1011;//矩阵列数
int rowx=539;//消元子行数
int rowb=263;//被消元行数


int **X;//消元子矩阵
int **B;//被消元行矩阵
int **Xe;//实际操作的消元子矩阵
int **Be;//实际操作的被消元行矩阵
int *Xe_index;//消元子首元素位置-子矩阵对应索引

string pathx="D:\\data\\1011_539_263\\X.txt";//消元子文件路径
string pathb="D:\\data\\1011_539_263\\B.txt";//被消元行文件路径

void ini_matrix()//矩阵初始化为0
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

void string_to_num(string s,int r,int c,int **a)//将字符串转化为数字并将矩阵中对应位置置1
{
    string temp;
    int num;
    stringstream ss(s);//字符串流
    while(ss>>temp){
        stringstream t;
        t<<temp;
        t>>num;
        a[r][c-num-1]=1;
    }
}

void readfile(){//读取文件
    ifstream x;//消元子
    ifstream b;//被消元行
    x.open(pathx,ios::in);
    if(!x.is_open()){
        cout<<"消元子文件读取失败"<<endl;
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
        cout<<"被消元行文件读取失败"<<endl;
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

void restart()//将两个需要执行的矩阵恢复成初始状态
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

int find_first(int *a,int l,int start=0)//找到第一个1所在的位置
{
    for(int i=start;i<l;i++){
        if(a[i]==1)
            return i;
    }
    return -1;
}

int b_is_in_x(int *a,int l)//找到消元子中对应的行
{
    for(int i=0;i<rowb;i++){
        if(find_first(Xe[i],l)==find_first(a,l)){
            return i;
        }
    }
    return -1;
}
//===========================串行算法==============================
void common()//串行算法
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
            int t=b_is_in_x(Be[i],col);//找到消元子的目标行
            QueryPerformanceCounter((LARGE_INTEGER*)&tail);
            sumtime2+=(tail - head) * 1000.0 / freq;
            //cout<<t<<" ";
            if(t!=-1){
                count_xiao++;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(int j=0;j<col;j++){
                    Be[i][j]=Be[i][j]^Xe[t][j];//异或操作
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime+=(tail - head) * 1000.0 / freq;
            }
            else{
                count_zhuan++;
                for(int j=0;j<col;j++){
                    Xe[xnum][j]=Be[i][j];//转化成消元行
                }
                xnum++;
                break;
            }
        }
    }
    cout<<"消元次数:"<<count_xiao<<endl;
    cout<<"消元子转化次数："<<count_zhuan<<endl;
    cout<<"消元时间："<<sumtime<<"ms"<<endl;
    cout<<"找消元子时间："<<sumtime2<<"ms"<<endl;
}
//=============================1：SSE==================================
void elimination1()
{
    int xnum=rowx-rowb;
    int count_xiao=0;
    int count_zhuan=0;
    long long head,tail,freq;
    long double sumtime=0;
    for(int i=0;i<rowb;i++){
        while(find_first(Be[i],col)!=-1){
            int t=b_is_in_x(Be[i],col);//找到消元子的目标行
            if(t!=-1){
                int j=0;
                count_xiao++;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(j=0;j+4<=col;j+=4){//异或操作
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
    cout<<"消元次数："<<count_xiao<<endl;
    cout<<"消元子转化次数："<<count_zhuan<<endl;
    cout<<"消元时间："<<sumtime<<"ms"<<endl;
}
//=====================================================================

//=============================2：AVX==================================
void elimination2()
{
    int xnum=rowx-rowb;
    long long head,tail,freq;
    long double sumtime=0;
    for(int i=0;i<rowb;i++){
        while(find_first(Be[i],col)!=-1){
            int t=b_is_in_x(Be[i],col);//找到消元子的目标行
            if(t!=-1){
                int j=0;
                QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                QueryPerformanceCounter((LARGE_INTEGER*)&head);
                for(j=0;j+8<=col;j+=8){//异或操作
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
                    Xe[xnum][j]=Be[i][j];//转化成消元行
                }
                xnum++;
                break;
            }
        }
    }
     cout<<"消元时间："<<sumtime<<"ms"<<endl;
}
//=====================================================================

//===========================3:找消元子优化==============================
void ini_xeindex()
{
    for(int i=0;i<rowx;i++){
        Xe_index[i]=-1;
    }
    for(int i=0;i<rowx;i++){
        for(int j=0;j<col;j++){
            if(Xe[i][j]==1){
                Xe_index[j]=i;//首元素为j的消元子在消元子第i行
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

int find_first2(int *a,int l,int &first,int start=0)//找到第一个1所在的位置
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
    long double sumtime=0;//消元时间
    long double sumtime2=0;//找消元子时间
    long double sumtime3=0;//转换成消元子时间
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
                    Be[i][j]=Be[i][j]^Xe[t][j];//异或操作
                }
                QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                sumtime+=(tail - head) * 1000.0 / freq;
            }
            else{
                count_zhuan++;
                for(int j=first;j<col;j++){
                    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
                    QueryPerformanceCounter((LARGE_INTEGER*)&head);
                    Xe[xnum][j]=Be[i][j];//转化成消元行
                    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
                    sumtime3+=(tail - head) * 1000.0 / freq;
                }
                Xe_index[first]=xnum;
                xnum++;
                break;
            }
        }
    }
    cout<<"消元次数："<<count_xiao<<endl;
    cout<<"消元子转化次数："<<count_zhuan<<endl;
    cout<<"消元时间："<<sumtime<<"ms"<<endl;
    cout<<"找消元子时间："<<sumtime2<<"ms"<<endl;
    cout<<"被消元行转化成消元子时间："<<sumtime3<<"ms"<<endl;
}
//===============================================================

//=======================4:AVX全优化==============================
void elimination4()
{
    int xnum=rowx-rowb;
    int count_xiao=0;
    int count_zhuan=0;
    long long head,tail,freq;
    long double sumtime=0;//消元时间
    long double sumtime2=0;//找消元子时间
    long double sumtime3=0;//转换成消元子时间
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
                for(j=0;j+8<=col;j+=8){//异或操作
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
    cout<<"消元次数："<<count_xiao<<endl;
    cout<<"消元子转化次数："<<count_zhuan<<endl;
    cout<<"消元时间："<<sumtime<<"ms"<<endl;
    cout<<"找消元子时间："<<sumtime2<<"ms"<<endl;
    cout<<"被消元行转化成消元子时间："<<sumtime3<<"ms"<<endl;
}
//=================================================================

//=======================5:SSE全优化==============================
void elimination5()
{
    int xnum=rowx-rowb;
    int count_xiao=0;
    int count_zhuan=0;
    long long head,tail,freq;
    long double sumtime=0;//消元时间
    long double sumtime2=0;//找消元子时间
    long double sumtime3=0;//转换成消元子时间
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
                for(j=0;j+4<=col;j+=4){//异或操作
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
    cout<<"消元次数："<<count_xiao<<endl;
    cout<<"消元子转化次数："<<count_zhuan<<endl;
    cout<<"消元时间："<<sumtime<<"ms"<<endl;
    cout<<"找消元子时间："<<sumtime2<<"ms"<<endl;
    cout<<"被消元行转化成消元子时间："<<sumtime3<<"ms"<<endl;
}
//=================================================================

//=============================6:pthread=======================================
const int NUM_THREADS = 4;
pthread_barrier_t barrier_getBe;//获取本轮被消元行集合
pthread_barrier_t barrier_Elimination;//多线程消元
pthread_barrier_t barrier_Transform;//转化成消元子

int *Be_status;//每个被消元行当前状态（0：为空 1：需要消元 -1：需要转化成消元子 -2:转化完成）
int **Be_Elimination;//当前需要进行消元的被消元行及其首元素集合
int **Be_transform;//当前需要转化成消元子的被消元行及其首元素集合

int finished_Elimination;//当前消元完成的被消元行数
int need_Elimination;//本轮需要消元的被消元行数
int need_transform;//本轮需要转化的被消元行数
int xnum=rowx-rowb;//消元子数量

typedef struct {
	int t_id;
}threadParm_t;

void initiate()//初始化
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

void clean() //delete指针
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
        //0号进程获取本轮消元的被消元行行号
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
        //同步
        pthread_barrier_wait(&barrier_getBe);
        //每个线程分配若干行进行消元并更新状态数组
        for(int i=t_id;i<need_Elimination;i+=NUM_THREADS){
            int curr=Be_Elimination[i][0];//当前操作被消元行行号
            Be_Elimination[i][1]=find_first(Be[curr],col);
            int first=Be_Elimination[i][1];//首元素
            int t=Xe_index[first];
            if(t!=-1){//消元子中有该行，直接消元
                for(int j=first;j<col;j++){
                    Be[curr][j]=Be[curr][j]^Xe[t][j];//异或操作
                }
                first=find_first(Be[curr],col,first);
                //cout<<curr<<":"<<first<<" ";
                if(first==-1){//消元后该行变成空
                    Be_status[curr]=0;
                    //cout<<"2 ";
                    finished_Elimination++;
                }
            }
            else{//消元子中没有该行，该行在下一轮转化成消元子
                Be_status[curr]=-1;
                //cout<<"1 ";
                finished_Elimination++;
            }
        }
        //同步
        pthread_barrier_wait(&barrier_Elimination);
        //0号线程将消元子转化
        if(t_id==0){
            for(int i=0;i<need_transform;i++){
                int curr=Be_transform[i][0];//当前操作行号
                int first=Be_transform[i][1];//首元素
                if(Xe_index[first]==-1){
                    for(int j=0;j<col;j++){
                        Xe[xnum][j]=Be[curr][j];
                    }
                    Xe_index[first]=xnum;
                    xnum++;
                }
                else{//该轮已经有对应首元素的行转化成消元子，下一轮重新消元
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
    //初始化barrier
    pthread_barrier_init(&barrier_getBe,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Transform,NULL,NUM_THREADS);
    //创建线程
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
        //0号进程获取本轮消元的被消元行行号
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
        //同步
        pthread_barrier_wait(&barrier_getBe);
        //每个线程分配若干行进行消元并更新状态数组
        for(int i=t_id;i<need_Elimination;i+=NUM_THREADS){
            int curr=Be_Elimination[i][0];//当前操作被消元行行号
            Be_Elimination[i][1]=find_first(Be[curr],col);
            int first=Be_Elimination[i][1];//首元素
            int t=Xe_index[first];
            if(t!=-1){//消元子中有该行，直接消元
                int j;
                for(j=first;j+8<=col;j+=8){//异或操作
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
                if(first==-1){//消元后该行变成空
                    Be_status[curr]=0;
                    //cout<<"2 ";
                    finished_Elimination++;
                }
            }
            else{//消元子中没有该行，该行在下一轮转化成消元子
                Be_status[curr]=-1;
                //cout<<"1 ";
                finished_Elimination++;
            }
        }
        //同步
        pthread_barrier_wait(&barrier_Elimination);
        //0号线程将消元子转化
        if(t_id==0){
            for (int i = 0; i < need_transform; i++) {
                int curr = Be_transform[i][0];//当前操作行号
                int first = Be_transform[i][1];//首元素
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
                else {//该轮已经有对应首元素的行转化成消元子，下一轮重新消元
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
    //初始化barrier
    pthread_barrier_init(&barrier_getBe,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Transform,NULL,NUM_THREADS);
    //创建线程
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
            int curr=Be_Elimination[i][0];//当前操作被消元行行号
            Be_Elimination[i][1]=find_first(Be[curr],col);
            int first=Be_Elimination[i][1];//首元素
            int t=Xe_index[first];
            if(t!=-1){//消元子中有该行，直接消元
                int j;
                for(j=first;j+8<=col;j+=8){//异或操作
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
                if(first==-1){//消元后该行变成空
                    Be_status[curr]=0;
                    //cout<<"2 ";
                    finished_Elimination++;
                }
            }
            else{//消元子中没有该行，该行在下一轮转化成消元子
                Be_status[curr]=-1;
                //cout<<"1 ";
                finished_Elimination++;
            }
        }
        #pragma omp single
        for (int i = 0; i < need_transform; i++) {
            int curr = Be_transform[i][0];//当前操作行号
            int first = Be_transform[i][1];//首元素
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
            else {//该轮已经有对应首元素的行转化成消元子，下一轮重新消元
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
            int curr=Be_Elimination[i][0];//当前操作被消元行行号
            Be_Elimination[i][1]=find_first(Be[curr],col);
            int first=Be_Elimination[i][1];//首元素
            int t=Xe_index[first];
            if(t!=-1){//消元子中有该行，直接消元
                for(int j=first;j<col;j++){
                    Be[curr][j]=Be[curr][j]^Xe[t][j];//异或操作
                }
                first=find_first(Be[curr],col,first);
                //cout<<curr<<":"<<first<<" ";
                if(first==-1){//消元后该行变成空
                    Be_status[curr]=0;
                    //cout<<"2 ";
                    finished_Elimination++;
                }
            }
            else{//消元子中没有该行，该行在下一轮转化成消元子
                Be_status[curr]=-1;
                //cout<<"1 ";
                finished_Elimination++;
            }
        }
        #pragma omp single
        for(int i=0;i<need_transform;i++){
            int curr=Be_transform[i][0];//当前操作行号
            int first=Be_transform[i][1];//首元素
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
    rowx+=rowb;//消元子总行数
    ini_matrix();
    readfile();
    //=========================平凡算法=====================================
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

    //=========================3:找消元子优化=====================================
    cout<<"3:找消元子优化==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination3();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================4:AVX全优化=====================================
    cout<<"4:AVX全优化==============================================="<<endl;
    restart();
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&head);
    elimination4();
    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
    cout <<"sum:"<< (tail - head) * 1000.0 / freq <<"ms"<< endl;
    //==============================================================

    //=========================5:SSE全优化=====================================
    cout<<"5:SSE全优化==============================================="<<endl;
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
