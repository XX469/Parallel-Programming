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
//��������
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
//==================ƽ���㷨=========================
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
//==================0:��̬�߳�========================
//�����߳����ݽṹ
typedef struct{
    int k;//��ȥ���ִ�
    int t_id;//�߳�ID
}threadParam_t1;

//�����̺߳���
void *threadFunc0(void *param){
    threadParam_t1 *p=(threadParam_t1*)param;
    int k=p->k;//��ȥ���ִ�
    int t_id=p->t_id;//�̱߳��
    int i=k+t_id+1;//��ȡ�Լ��ļ�������
    for(int j=k+1;j<n;++j){
        m[i][j]=m[i][j]-m[i][k]*m[k][j];
    }
    m[i][k]=0;
    pthread_exit(NULL);
}

void Gaussian_elimination0(){
    for(int k=0;k<n;k++){//���̳߳�������
        for(int j=k+1;j<n;j++){
            m[k][j]=m[k][j]/m[k][k];//division
        }
        m[k][k]=1.0;
        //���������̣߳�������ȥ����
        int worker_count=n-1-k;//�����߳�����
        pthread_t* handles=new pthread_t[worker_count];//������Ӧ��handle
        threadParam_t1* param=new threadParam_t1[worker_count];//������Ӧ���߳����ݽṹ
        for(int t_id=0;t_id<worker_count;t_id++){//��������
            param[t_id].k=k;
            param[t_id].t_id=t_id;
        }
        for(int t_id=0;t_id<worker_count;t_id++){//�����߳�
            pthread_create(&handles[t_id],NULL,threadFunc0,(void*)&param[t_id]);
        }
        for(int t_id=0;t_id<worker_count;t_id++){//���̹߳���ȴ����еĹ����߳���ɴ�����ȥ��
            pthread_join(handles[t_id],0);
        }
    }
}
//==============================================================================
//�̶߳���
const int NUM_THREADS = 4;
typedef struct {
	int t_id;
}threadParm_t;

//=====================1:��̬�߳� + �ź���ͬ��========================
//�ź�������
sem_t sem_main;
sem_t sem_workerstart[NUM_THREADS]; // ÿ���߳����Լ�ר�����ź���
sem_t sem_workerend[NUM_THREADS];

void *threadFunc1(void *param){
    threadParm_t *p = (threadParm_t*)param;
    int t_id=p->t_id;
    for(int k=0;k<n;++k){
        sem_wait(&sem_workerstart[t_id]); // �������ȴ�������ɳ�������
        //ѭ����������
        for(int i=k+1+t_id;i<n;i+= NUM_THREADS){
            for(int j=k+1;j<n;++j){
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[k][k]=1.0;
        }
        sem_post(&sem_main); // �������߳�
        sem_wait(&sem_workerend[t_id]); //�������ȴ����̻߳��ѽ�����һ��
    }
    pthread_exit(NULL);
}
void Gaussian_elimination1(){
    //��ʼ��
    sem_init(&sem_main, 0, 0);
    for(int i=0;i<NUM_THREADS;++i){
        sem_init(&sem_workerstart[i], 0, 0);
        sem_init(&sem_workerend[i], 0, 0);
    }
    //�����߳�
    pthread_t handles[NUM_THREADS];
    threadParm_t param[NUM_THREADS];
    for(int t_id=0;t_id<NUM_THREADS;t_id++){
        param[t_id].t_id=t_id;
        pthread_create(&handles[t_id],NULL,threadFunc1,(void*)&param[t_id]);
    }
    for(int k=0;k<n;++k){
        //���߳�����������
        for(int j=k+1;j<n;j++){
            m[k][j]=m[k][j]/m[k][k];
        }
        m[k][k]=1.0;
        //��ʼ���ѹ����߳�
        for(int t_id=0;t_id<NUM_THREADS;++t_id){
            sem_post(&sem_workerstart[t_id]);
        }
        //���߳�˯�ߣ��ȴ����еĹ����߳���ɴ�����ȥ����
        for(int t_id=0;t_id<NUM_THREADS;++t_id){
            sem_wait(&sem_main);
        }
        // ���߳��ٴλ��ѹ����߳̽�����һ�ִε���ȥ����
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

//=================2:��̬�߳� +barrier ͬ��========================
//barrier
pthread_barrier_t barrier_Division;
pthread_barrier_t barrier_Elimination;

//�̺߳�������
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
        //��һ��ͬ����
        pthread_barrier_wait(&barrier_Division);
        //ѭ����������
        for(int i=k+1+t_id;i<n;i+=NUM_THREADS){
            for(int j=k+1;j<n;++j){
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0.0;
        }
        // �ڶ���ͬ����
        pthread_barrier_wait(&barrier_Elimination);
    }
    pthread_exit(NULL);
}
void Gaussian_elimination2(){
    //��ʼ��barrier
    pthread_barrier_init(&barrier_Division,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
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
    pthread_barrier_destroy(&barrier_Division);
    pthread_barrier_destroy(&barrier_Elimination);
}
//========================================================================

//======3:��̬�߳� + �ź���ͬ���汾 + ����ѭ��ȫ�������̺߳���===========
//�ź�������
sem_t sem_leader;
sem_t sem_Division[NUM_THREADS-1];
sem_t sem_Elimination[NUM_THREADS-1];
//�̺߳�������
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
    //��ʼ��barrier
    sem_init(&sem_leader,0,0);
    for(int i=0;i<NUM_THREADS-1;++i){
        sem_init(&sem_Division[i],0,0);
        sem_init(&sem_Elimination[i],0,0);
    }
    //�����߳�
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


//=================4:��̬�߳� +barrier ͬ��+���л���========================
//�̺߳�������
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
        //��һ��ͬ����
        pthread_barrier_wait(&barrier_Division);
        //���л��֣�ÿ���̷߳���������
        for(int i=k+1;i<n;i++){
            for(int j=k+1+t_id;j<n;j+=NUM_THREADS){
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0.0;
        }
        // �ڶ���ͬ����
        pthread_barrier_wait(&barrier_Elimination);
    }
    pthread_exit(NULL);
}
void Gaussian_elimination4(){
    //��ʼ��barrier
    pthread_barrier_init(&barrier_Division,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
    //�����߳�
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

//=================5:��̬�߳� +barrier ͬ��+AVX========================
//�̺߳�������
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
        //��һ��ͬ����
        pthread_barrier_wait(&barrier_Division);
        //ѭ����������
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
        // �ڶ���ͬ����
        pthread_barrier_wait(&barrier_Elimination);
    }
    pthread_exit(NULL);
}
void Gaussian_elimination5(){
    //��ʼ��barrier
    pthread_barrier_init(&barrier_Division,NULL,NUM_THREADS);
    pthread_barrier_init(&barrier_Elimination,NULL,NUM_THREADS);
    //�����߳�
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

//==============================8:OpenMP+�л���================================
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
        //=========================ƽ���㷨===============================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_common();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //==================0:��̬�߳�========================
        //Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        //Gaussian_elimination1();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        //cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //==================1:��̬�߳� + �ź���ͬ��====================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination1();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //=================2:��̬�߳� +barrier ͬ��========================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination2();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //======3:��̬�߳� + �ź���ͬ���汾 + ����ѭ��ȫ�������̺߳���===========
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination3();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //=================4:��̬�߳� +barrier ͬ��+���л���========================
        Create_Matrix();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        Gaussian_elimination4();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        cout << (tail - head) * 1000.0 / freq << endl;
        //=================================================================

        //=================5:��̬�߳� +barrier ͬ��+AVX========================
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

        //=================8:OpenMP+���л���========================
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
