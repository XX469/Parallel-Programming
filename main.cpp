#include<iostream>
#include<math.h>
#include<cstdlib>
#include<windows.h>
#include<stdlib.h>
#include<xmmintrin.h>
#include<fstream>
#include<immintrin.h>
#include<pthread.h>
#include<semaphore.h>

using namespace std;

const long long int N = 10000;
float m[N][N];
int n;
//创建矩阵
void Create_Matrix()
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
//==================平凡算法=========================
void Gaussian_common(){
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
//===================================================
//==================0:动态线程========================
//定义线程数据结构
typedef struct{
    int k;//消去的轮次
    int t_id;//线程ID
}threadParam_t1;

//定义线程函数
void *threadFunc0(void *param){
    threadParam_t1 *p=(threadParam_t1*)param;
    int k=p->k;//消去的轮次
    int t_id=p->t_id;//线程编号
    int i=k+t_id+1;//获取自己的计算任务
    for(int j=k+1;j<n;++j){
        m[i][j]=m[i][j]-m[i][k]*m[k][j];
    }
    m[i][k]=0;
    pthread_exit(NULL);
}

void Gaussian_elimination0(){
    for(int k=0;k<n;k++){//主线程除法操作
        for(int j=k+1;j<n;j++){
            m[k][j]=m[k][j]/m[k][k];//division
        }
        m[k][k]=1.0;
        //创建工作线程，进行消去操作
        int worker_count=n-1-k;//工作线程数量
        pthread_t* handles=new pthread_t[worker_count];//创建对应的handle
        threadParam_t1* param=new threadParam_t1[worker_count];//创建对应的线程数据结构
        for(int t_id=0;t_id<worker_count;t_id++){//分配任务
            param[t_id].k=k;
            param[t_id].t_id=t_id;
        }
        for(int t_id=0;t_id<worker_count;t_id++){//创建线程
            pthread_create(&handles[t_id],NULL,threadFunc0,(void*)&param[t_id]);
        }
        for(int t_id=0;t_id<worker_count;t_id++){//主线程挂起等待所有的工作线程完成此轮消去工
            pthread_join(handles[t_id],0);
        }
    }
}
//==============================================================================
//线程定义
const int NUM_THREADS = 4;
typedef struct {
	int t_id;
}threadParm_t;

//=====================1:静态线程 + 信号量同步========================
//信号量定义
sem_t sem_main;
sem_t sem_workerstart[NUM_THREADS]; // 每个线程有自己专属的信号量
sem_t sem_workerend[NUM_THREADS];

void *threadFunc1(void *param){
    threadParm_t *p = (threadParm_t*)param;
    int t_id=p->t_id;
    for(int k=0;k<n;++k){
        sem_wait(&sem_workerstart[t_id]); // 阻塞，等待主线完成除法操作
        //循环划分任务
        for(int i=k+1+t_id;i<n;i+= NUM_THREADS){
            for(int j=k+1;j<n;++j){
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[k][k]=1.0;
        }
        sem_post(&sem_main); // 唤醒主线程
        sem_wait(&sem_workerend[t_id]); //阻塞，等待主线程唤醒进入下一轮
    }
    pthread_exit(NULL);
}
void Gaussian_elimination1(){
    //初始化
    sem_init(&sem_main, 0, 0);
    for(int i=0;i<NUM_THREADS;++i){
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }
    //创建线程
    pthread_t handles[NUM_THREADS];
    threadParm_t param[NUM_THREADS];
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        param[t_id].t_id=t_id;
        pthread_create(&handles[t_id],NULL,threadFunc1,(void*)&param[t_id]);
    }
    for(int k=0;k<n;++k){
        //主线程做除法操作
        for(int j=k+1;j<n;j++){
            m[k][j]=m[k][j]/m[k][k];
        }
        m[k][k]=1.0;
        //开始唤醒工作线程
        for(int t_id=0;t_id<NUM_THREADS;++t_id){
            sem_post(&sem_workerstart[t_id]);
        }
        //主线程睡眠（等待所有的工作线程完成此轮消去任务）
        for(int t_id=0;t_id<NUM_THREADS;++t_id){
            sem_wait(&sem_main);
        }
        // 主线程再次唤醒工作线程进入下一轮次的消去任务
        for(int t_id=0;t_id<NUM_THREADS;++t_id){
            sem_post(&sem_workerend[t_id]);
        }
    }
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        pthread_join(handles[t_id],0);
    }
    sem_destroy(&sem_main);
}
//====================================================================

//=================2:静态线程 +barrier 同步========================
//barrier
pthread_barrier_t barrier_Division;
pthread_barrier_t barrier_Elimination;

//线程函数定义
void *threadFunc2(void *param){
    threadParm_t *p = (threadParm_t*)param;
    int t_id=p->t_id;
    for(int k=0;k<n;++k){
        if(t_id==0){
            for(int j=k+1;j<n;j++){
                m[k][j]=m[k][j]/m[k][k];
            }
            m[k][k]=1.0;
        }
        //第一个同步点
        pthread_barrier_wait(&barrier_Division);
        //循环划分任务
        for(int i=k+1+t_id;i<n;i+=NUM_THREADS){
            for(int j=k+1;j<n;++j){
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0.0;
        }
        // 第二个同步点
        pthread_barrier_wait(&barrier_Elimination);
    }
    pthread_exit(NULL);
}
void Gaussian_elimination2(){
    //初始化barrier
    pthread_barrier_init(&barrier_Division,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
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
    pthread_barrier_destroy(&barrier_Division);
    pthread_barrier_destroy(&barrier_Elimination);
}
//========================================================================

//======3:静态线程 + 信号量同步版本 + 三重循环全部纳入线程函数===========
//信号量定义
sem_t sem_leader;
sem_t sem_Division[NUM_THREADS-1];
sem_t sem_Elimination[NUM_THREADS-1];
//线程函数定义
void *threadFunc3(void *param){
    threadParm_t *p = (threadParm_t*)param;
    int t_id=p->t_id;
    for(int k=0;k<n;++k){
        if(t_id==0){
            for(int j=k+1;j<n;j++){
                m[k][j]=m[k][j]/m[k][k];
            }
            m[k][k]=1.0;
        }
        else{
            sem_wait(&sem_Division[t_id-1]);
        }
        if(t_id==0){
            for(int i=0;i<NUM_THREADS-1;++i){
                sem_post(&sem_Division[i]);
            }
        }
        for(int i=k+1+t_id;i<n;i+=NUM_THREADS){
            for(int j=k+1;j<n;++j){
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0.0;
        }
        if(t_id==0){
            for(int i=0;i<NUM_THREADS-1;++i){
                sem_wait(&sem_leader);
            }
            for(int i=0;i<NUM_THREADS-1;++i){
                sem_post(&sem_Elimination[i]);
            }
        }
        else{
            sem_post(&sem_leader);
            sem_wait(&sem_Elimination[t_id-1]);
        }
    }
    pthread_exit(NULL);
}
void Gaussian_elimination3(){
    //初始化barrier
    sem_init(&sem_leader,0,0);
    for(int i=0;i<NUM_THREADS-1;++i){
        sem_init(&sem_Division[i],0,0);
        sem_init(&sem_Elimination[i],0,0);
    }
    //创建线程
    pthread_t handles[NUM_THREADS];
    threadParm_t param[NUM_THREADS];
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        param[t_id].t_id=t_id;
        pthread_create(&handles[t_id],NULL,threadFunc3,(void*)&param[t_id]);
    }
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        pthread_join(handles[t_id],0);
    }
    sem_destroy(&sem_leader);
}
//============================================================================


//=================4:静态线程 +barrier 同步+按列划分========================
//线程函数定义
void *threadFunc4(void *param){
    threadParm_t *p = (threadParm_t*)param;
    int t_id=p->t_id;
    for(int k=0;k<n;++k){
        if(t_id==0){
            for(int j=k+1;j<n;j++){
                m[k][j]=m[k][j]/m[k][k];
            }
            m[k][k]=1.0;
        }
        //第一个同步点
        pthread_barrier_wait(&barrier_Division);
        //按列划分，每个线程分配若干列
        for(int i=k+1;i<n;i++){
            for(int j=k+1+t_id;j<n;j+=NUM_THREADS){
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0.0;
        }
        // 第二个同步点
        pthread_barrier_wait(&barrier_Elimination);
    }
    pthread_exit(NULL);
}
void Gaussian_elimination4(){
    //初始化barrier
    pthread_barrier_init(&barrier_Division,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
    //创建线程
    pthread_t handles[NUM_THREADS];
    threadParm_t param[NUM_THREADS];
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        param[t_id].t_id=t_id;
        pthread_create(&handles[t_id],NULL,threadFunc4,(void*)&param[t_id]);
    }
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        pthread_join(handles[t_id],0);
    }
    pthread_barrier_destroy(&barrier_Division);
    pthread_barrier_destroy(&barrier_Elimination);
}
//========================================================================

//=================5:静态线程 +barrier 同步+AVX========================
//线程函数定义
void *threadFunc5(void *param){
    threadParm_t *p = (threadParm_t*)param;
    __m256 vt,va;
    __m256 vaik,vakj,vaij,vx;
    int t_id=p->t_id;
    for(int k=0;k<n;++k){
        if(t_id==0){
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
        }
        //第一个同步点
        pthread_barrier_wait(&barrier_Division);
        //循环划分任务
        for(int i=k+1+t_id;i<n;i+=NUM_THREADS){
            int j;
            vaik=_mm256_set1_ps(m[i][k]);//dupToVector4(m[i,k])
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
            m[i][k]=0.0;
        }
        // 第二个同步点
        pthread_barrier_wait(&barrier_Elimination);
    }
    pthread_exit(NULL);
}
void Gaussian_elimination5(){
    //初始化barrier
    pthread_barrier_init(&barrier_Division,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
    //创建线程
    pthread_t handles[NUM_THREADS];
    threadParm_t param[NUM_THREADS];
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        param[t_id].t_id=t_id;
        pthread_create(&handles[t_id],NULL,threadFunc5,(void*)&param[t_id]);
    }
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        pthread_join(handles[t_id],0);
    }
    pthread_barrier_destroy(&barrier_Division);
    pthread_barrier_destroy(&barrier_Elimination);
}
//========================================================================

//==============================6:OpenMP================================
void Gaussian_elimination6(){
    double tmp;
    int i,j,k;
    #pragma omp parallel num_threads(NUM_THREADS) private(i,j,k,tmp)
    for(k=0;k<n;k++){
        #pragma omp single
        tmp=m[k][k];
        for(j=k+1;j<n;j++){
            m[k][j]=m[k][j]/tmp;//division
        }
        m[k][k]=1.0;
        #pragma omp for
        for(i=k+1;i<n;i++){
            tmp=m[i][k];
            for(j=k+1;j<n;j++){
                m[i][j]=m[i][j]-tmp*m[k][j];//subtraction
            }
            m[i][k]=0;
        }
    }
}
//======================================================================

//==============================7:OpenMP+AVX================================
void Gaussian_elimination7(){
    int i,j,k;
    __m256 vt,va;
    __m256 vaik,vakj,vaij,vx;
    #pragma omp parallel num_threads(NUM_THREADS) private(i,j,k,vaik,vakj,vaij,vx)
    for(k=0;k<n;k++){
        #pragma omp single
        vt = _mm256_set1_ps(m[k][k]);//dupTo4Float(m[k,k])
        for(j=k+1;j+8<=n;j+=8){
            va = _mm256_loadu_ps(&m[k][j]);//load4FloatFrom(&m[k,j])
            va = _mm256_div_ps(va,vt);//va/vt
            _mm256_storeu_ps(&m[k][j],va);//store4FloatTo(&m[k,j],va)
        }
        for(j=j ;j<n;j++){//The remaining part
            m[k][j]=m[k][j]/m[k][k];
        }
        m[k][k]=1.0;
        #pragma omp for
        for(i=k+1;i<n;i++){
            vaik=_mm256_set1_ps(m[i][k]);//dupToVector4(m[i,k])
            for(j=k+1;j+8<=n;j+=8){
                vakj= _mm256_loadu_ps(&m[k][j]);//load4FloatFrom(&m[k,j])
                vaij= _mm256_loadu_ps(&m[i][j]);//load4FloatFrom(&m[i,j])
                vx=_mm256_mul_ps(vakj,vaik);//mul vakj*vaik
                vaij=_mm256_sub_ps(vaij,vx);//sub vaij-vx
                _mm256_storeu_ps(&m[i][j],vaij);//store4FloatTo(&m[i,j],vaij)
            }
            for(j=j;j<n;j++){//The remaining part
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0;
        }
    }
}
//======================================================================

//==============================8:OpenMP+列划分================================
void Gaussian_elimination8(){
    double tmp;
    int i,j,k;
    #pragma omp parallel num_threads(NUM_THREADS) private(i,j,k,tmp)
    for(k=0;k<n;k++){
        #pragma omp single
        tmp=m[k][k];
        for(j=k+1;j<n;j++){
            m[k][j]=m[k][j]/tmp;//division
        }
        m[k][k]=1.0;
        #pragma omp for
        for(j=k+1;j<n;j++){
            for(i=n-1;i>=k+1;i--){
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0;
        }
    }
}
//======================================================================
int main()
{
    int datasize[9]={10,50,100,200,500,1000,1500,2000,4000};
    int threadsize[6]={2,4,6,8,10,12};
    long long head,tail,freq;
    cout<<"===============test for n:==================="<<endl;
    for (int i = 7; i < 8; i++) {
        n = datasize[i];
        cout << "n:" << n << endl;
        //=========================平凡算法===============================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_common();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //==================0:动态线程========================
        //Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        //Gaussian_elimination1();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        //cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //==================1:静态线程 + 信号量同步====================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination1();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //=================2:静态线程 +barrier 同步========================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination2();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //======3:静态线程 + 信号量同步版本 + 三重循环全部纳入线程函数===========
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination3();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //=================4:静态线程 +barrier 同步+按列划分========================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination4();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //=================5:静态线程 +barrier 同步+AVX========================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination5();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //=================6:OpenMP========================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination6();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //=================7:OpenMP+AVX========================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination7();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //=================8:OpenMP+按列划分========================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination8();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================
        cout<<"============================================================"<<endl;
    }
    return 0;
}
