#include <iostream>
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <sys/time.h>
#include <arm_neon.h>

using namespace std;

const long long int N = 10000;
float m[N][N];
int n;
int NUM_THREADS = 8;

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
//==================1:块划分（行）========================
double MPI_Function1(){
    double start_time,end_time;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //报告识别调用进程的rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //报告进程数
    if(rank==0){
        Create_Matrix();
    }
    start_time=MPI_Wtime();
    int rownum=n/size; //每个进程负责子矩阵行数
    if(rank==0){//0号进程分发任务
        for(int i=1;i<size;i++){//计算起始行和结束行
            int r1=i*rownum;
            int r2;
            if(i==size-1)
                r2=n;
            else
                r2=(i+1)*rownum;
            MPI_Send(&m[r1][0], (r2-r1)*n, MPI_FLOAT,i, 0, MPI_COMM_WORLD);
        }
    }
    else{//非0号进程接收任务
        if(rank==size-1){
            MPI_Recv(&m[rank*rownum][0],(n-rank*rownum)*n,MPI_FLOAT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        else{
            MPI_Recv(&m[rank*rownum][0],rownum*n,MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    int r1=rank*rownum;
    int r2=min((rank+1)*rownum,n);
    for(int k=0;k<n;k++){
        //除法
        if(k>=r1&&k<r2){
            for(int j=k+1;j<n;j++){
                m[k][j]/=m[k][k];
            }
            m[k][k]=1;
            for(int j=0;j<size;j++){
                MPI_Send(&m[k][0], n, MPI_FLOAT, j, 0, MPI_COMM_WORLD);
            }
        }
        else{
            MPI_Recv(&m[k][0], n, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //消元
        for(int i=max(k+1,r1);i<r2;i++){
            for(int j=k+1;j<n;j++){
                m[i][j]=m[i][j]-m[i][k]*m[k][j];
            }
            m[i][k]=0;
        }
    }
    end_time=MPI_Wtime();
    double time_sum = (end_time - start_time) * 1000;
    return time_sum;
}
//======================================================================
//=======================2:循环划分========================================
double MPI_Function2() {
    double start_time, end_time;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //报告识别调用进程的rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //报告进程数
    if (rank == 0) {
        Create_Matrix();
    }
    start_time = MPI_Wtime();
    int numtasks = 0;//每个进程的任务量
    if (rank < n % size)
        numtasks = n / size + 1;
    else
        numtasks = n / size;
    //存储每个进程对应任务的数组并由0号进程发送
    float* task = new float[numtasks * n];
    if (rank == 0) {
        for (int t = 1; t < size; t++) {
            for (int i = t; i < n; i += size) {
                for (int j = 0; j < n; j++) {
                    task[i / size * n + j] = m[i][j];
                }
            }
            int num = t < n% size ? n / size + 1 : n / size;//对应进程的任务量
            MPI_Send(task, num * n, MPI_FLOAT, t, 0, MPI_COMM_WORLD);
        }
    }
    //非0号进程接收任务
    else {
        MPI_Recv(&m[rank][0], numtasks*n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //将数据从接收的行复制到其他行
        for (int i = 0; i < numtasks; i++) {
            for (int j = 0; j < n; j++) {
                m[rank + i * size][j] = m[rank + i][j];
            }
        }
    }

    for (int k = 0; k < n; k++) {
        //除法
        if (rank == k%size) {
            for (int j = k + 1; j < n; j++) {
                m[k][j] /= m[k][k];
            }
            m[k][k] = 1;
            for (int t = 0; t < size; t++) {
                if (t != rank) {
                    MPI_Send(&m[k][0], n, MPI_FLOAT, t, 1, MPI_COMM_WORLD);
                }
            }
        }
        else {//其他进程接收除法结果
            MPI_Recv(&m[k][0], n, MPI_FLOAT, k % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //消元
        for (int i = k + 1; i < n; i ++) {
            if ((i-rank)%size==0) {
                for (int j = k + 1; j < n; j++) {
                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
                }
                m[i][k] = 0;
            }
        }
    }
    end_time = MPI_Wtime();
    double time_sum= (end_time - start_time) * 1000;
    return time_sum;
}
//=================================================================

//=======================3:MPI+SIMD========================================
double MPI_Function3() {
    double start_time, end_time;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //报告识别调用进程的rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //报告进程数
    if (rank == 0) {
        Create_Matrix();
    }
    start_time = MPI_Wtime();
    int numtasks = 0;//每个进程的任务量
    if (rank < n % size)
        numtasks = n / size + 1;
    else
        numtasks = n / size;
    //存储每个进程对应任务的数组并由0号进程发送
    float* task = new float[numtasks * n];
    if (rank == 0) {
        for (int t = 1; t < size; t++) {
            for (int i = t; i < n; i += size) {
                for (int j = 0; j < n; j++) {
                    task[i / size * n + j] = m[i][j];
                }
            }
            int num = t < n% size ? n / size + 1 : n / size;//对应进程的任务量
            MPI_Send(task, num * n, MPI_FLOAT, t, 0, MPI_COMM_WORLD);
        }
    }
    //非0号进程接收任务
    else {
        MPI_Recv(&m[rank][0], numtasks * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //将数据从接收的行复制到其他行
        for (int i = 0; i < numtasks; i++) {
            for (int j = 0; j < n; j++) {
                m[rank + i * size][j] = m[rank + i][j];
            }
        }
    }
    float32x4_t vt, va, vaik, vakj, vaij, vx;
    for (int k = 0; k < n; k++) {
        //除法
        if (rank == k % size) {
            int j;
            vt = vdupq_n_f32(m[k][k]);
            for (j = k + 1; j + 4 <= n; j += 4) {
                va = vld1q_f32(&m[k][j]);
                va = vdivq_f32(va, vt);
                vst1q_f32(&m[k][j], va);
            }
            for (int t = j; t < n; t++) {
                m[k][t] = m[k][t] / m[k][k];
            }
            m[k][k] = 1.0;
            for (int t = 0; t < size; t++) {
                if (t != rank) {
                    MPI_Send(&m[k][0], n, MPI_FLOAT, t, 1, MPI_COMM_WORLD);
                }
            }
        }
        else {//其他进程接收除法结果
            MPI_Recv(&m[k][0], n, MPI_FLOAT, k % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //消元
        for (int i = k + 1; i < n; i++) {
            if ((i - rank) % size == 0) {
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
    end_time = MPI_Wtime();
    double time_sum = (end_time - start_time) * 1000;
    return time_sum;
}
//==================================================================
//=======================4:MPI+OpenMP========================================
double MPI_Function4() {
    double start_time, end_time;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //报告识别调用进程的rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //报告进程数
    if (rank == 0) {
        Create_Matrix();
    }
    start_time = MPI_Wtime();
    int numtasks = 0;//每个进程的任务量
    if (rank < n % size)
        numtasks = n / size + 1;
    else
        numtasks = n / size;
    //存储每个进程对应任务的数组并由0号进程发送
    float* task = new float[numtasks * n];
    if (rank == 0) {
        for (int t = 1; t < size; t++) {
            for (int i = t; i < n; i += size) {
                for (int j = 0; j < n; j++) {
                    task[i / size * n + j] = m[i][j];
                }
            }
            int num = t < n% size ? n / size + 1 : n / size;//对应进程的任务量
            MPI_Send(task, num * n, MPI_FLOAT, t, 0, MPI_COMM_WORLD);
        }
    }
    //非0号进程接收任务
    else {
        MPI_Recv(&m[rank][0], numtasks * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //将数据从接收的行复制到其他行
        for (int i = 0; i < numtasks; i++) {
            for (int j = 0; j < n; j++) {
                m[rank + i * size][j] = m[rank + i][j];
            }
        }
    }
    int i, j, k;
    double tmp;
#pragma omp parallel num_threads(NUM_THREADS) private(i,j,k,tmp)
    for (k = 0; k < n; k++) {
        //除法
#pragma omp single
        tmp = m[k][k];
        if (rank == k % size) {
            for (j = k + 1; j < n; j++) {
                m[k][j] = m[k][j] / tmp;//division
            }
            m[k][k] = 1;
            for (int t = 0; t < size; t++) {
                if (t != rank) {
                    MPI_Send(&m[k][0], n, MPI_FLOAT, t, 1, MPI_COMM_WORLD);
                }
            }
        }
        else {//其他进程接收除法结果
            MPI_Recv(&m[k][0], n, MPI_FLOAT, k % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //消元
#pragma omp for
        for (int i = k + 1; i < n; i++) {
            if ((i - rank) % size == 0) {
                tmp = m[i][k];
                for (int j = k + 1; j < n; j++) {
                    m[i][j] = m[i][j] - tmp * m[k][j];//subtraction
                }
                m[i][k] = 0;
            }
        }
    }
    end_time = MPI_Wtime();
    double time_sum = (end_time - start_time) * 1000;
    return time_sum;
}
//=================================================================
//=======================5:MPI+OpenMP+SIMD========================================
double MPI_Function5() {
    double start_time, end_time;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //报告识别调用进程的rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //报告进程数
    if (rank == 0) {
        Create_Matrix();
    }
    start_time = MPI_Wtime();
    int numtasks = 0;//每个进程的任务量
    if (rank < n % size)
        numtasks = n / size + 1;
    else
        numtasks = n / size;
    //存储每个进程对应任务的数组并由0号进程发送
    float* task = new float[numtasks * n];
    if (rank == 0) {
        for (int t = 1; t < size; t++) {
            for (int i = t; i < n; i += size) {
                for (int j = 0; j < n; j++) {
                    task[i / size * n + j] = m[i][j];
                }
            }
            int num = t < n% size ? n / size + 1 : n / size;//对应进程的任务量
            MPI_Send(task, num * n, MPI_FLOAT, t, 0, MPI_COMM_WORLD);
        }
    }
    //非0号进程接收任务
    else {
        MPI_Recv(&m[rank][0], numtasks * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //将数据从接收的行复制到其他行
        for (int i = 0; i < numtasks; i++) {
            for (int j = 0; j < n; j++) {
                m[rank + i * size][j] = m[rank + i][j];
            }
        }
    }
    int i, j, k;
    double tmp;
    float32x4_t vt, va, vaik, vakj, vaij, vx;
#pragma omp parallel num_threads(NUM_THREADS) private(i,j,k,tmp)
    for (k = 0; k < n; k++) {
        //除法
#pragma omp single
        vt = vdupq_n_f32(m[k][k]);
        if (rank == k % size) {
            for (j = k + 1; j + 4 <= n; j += 4) {
                va = vld1q_f32(&m[k][j]);
                va = vdivq_f32(va, vt);
                vst1q_f32(&m[k][j], va);
            }
            for (int t = j; t < n; t++) {
                m[k][t] = m[k][t] / m[k][k];
            }
            m[k][k] = 1;
            for (int t = 0; t < size; t++) {
                if (t != rank) {
                    MPI_Send(&m[k][0], n, MPI_FLOAT, t, 1, MPI_COMM_WORLD);
                }
            }
        }
        else {//其他进程接收除法结果
            MPI_Recv(&m[k][0], n, MPI_FLOAT, k % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //消元
#pragma omp for
        for (int i = k + 1; i < n; i++) {
            if ((i - rank) % size == 0) {
                vaik = vdupq_n_f32(m[i][k]);
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
    end_time = MPI_Wtime();
    double time_sum = (end_time - start_time) * 1000;
    return time_sum;
}
//=================================================================
//==========================6：流水线划分=================================
double MPI_Function6() {
    double start_time, end_time;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //报告识别调用进程的rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //报告进程数
    if (rank == 0) {
        Create_Matrix();
    }
    start_time = MPI_Wtime();
    int numtasks = 0;//每个进程的任务量
    if (rank < n % size)
        numtasks = n / size + 1;
    else
        numtasks = n / size;
    //存储每个进程对应任务的数组并由0号进程发送
    float* task = new float[numtasks * n];
    if (rank == 0) {
        for (int t = 1; t < size; t++) {
            for (int i = t; i < n; i += size) {
                for (int j = 0; j < n; j++) {
                    task[i / size * n + j] = m[i][j];
                }
            }
            int num = t < n% size ? n / size + 1 : n / size;//对应进程的任务量
            MPI_Send(task, num * n, MPI_FLOAT, t, 0, MPI_COMM_WORLD);
        }
    }
    //非0号进程接收任务
    else {
        MPI_Recv(&m[rank][0], numtasks * n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //将数据从接收的行复制到其他行
        for (int i = 0; i < numtasks; i++) {
            for (int j = 0; j < n; j++) {
                m[rank + i * size][j] = m[rank + i][j];
            }
        }
    }
    int Next = (rank + 1) % size;//下一个进程
    int Previous = (rank + size - 1) % size;//上一个进程
    for (int k = 0; k < n; k++) {
        //除法
        if (rank == k % size) {
            for (int j = k + 1; j < n; j++) {
                m[k][j] /= m[k][k];
            }
            m[k][k] = 1;
            //将除法结果广播给下一个进程
            MPI_Send(&m[k][0], n, MPI_FLOAT, Next, 1, MPI_COMM_WORLD);
        }
        else {//其他进程接收上一个进程传送的结果
            MPI_Recv(&m[k][0], n, MPI_FLOAT, Previous, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (Next != k % size) {//继续传给下一个
                MPI_Send(&m[k][0], n, MPI_FLOAT, Next, 1, MPI_COMM_WORLD);
            }
        }
        //消元
        for (int i = k + 1; i < n; i++) {
            if ((i - rank) % size == 0) {
                for (int j = k + 1; j < n; j++) {
                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
                }
                m[i][k] = 0;
            }
        }
    }
    end_time = MPI_Wtime();
    double time_sum = (end_time - start_time) * 1000;
    return time_sum;
}
//======================================================================================

int main()
{
    MPI_Init(nullptr,nullptr);
    int datasize[8]={50,100,200,500,1000,1500,2000,3000};
    for (int i = 0; i < 1; i++) {
        n = datasize[i];
        int rank,size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank==0){
            cout << "n:" << n << endl;
        }
        struct timeval start {};
        struct timeval end {};
        //=========================平凡算法===============================
        Create_Matrix();
        gettimeofday(&start, nullptr);
        Gaussian_common();
        gettimeofday(&end, nullptr);
        if(rank==0){
            cout <<"common:"<<((end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec)) * 1.0 / 1000 << endl;
        }
        //==========================1:块划分（行）===============================
        MPI_Barrier(MPI_COMM_WORLD);//进程同步
        double time1 = MPI_Function1();
        if(rank==0){
            cout<<"1:" << time1 << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);//进程同步
        //=================================================================
        //==========================2:循环划分===============================
        
        double time2 = MPI_Function2();
        if (rank == 0) {
            cout <<"2:" << time2 << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);//进程同步
        //=================================================================
        //==========================3:MPI+SIMD===============================
        
        double time3 = MPI_Function3();
        if (rank == 0) {
            cout <<"3:" << time3 << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);//进程同步
        //=================================================================
        //==========================4:MPI+OpenMP===============================

        double time4 = MPI_Function4();
        if (rank == 0) {
            cout <<"4:" << time4 << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);//进程同步
        //=================================================================
        //==========================5:MPI+OpenMP+SIMD===============================

        double time5 = MPI_Function5();
        if (rank == 0) {
            cout << "5:" << time5 << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);//进程同步
        //=================================================================
        //==========================6:流水线划分===============================

        double time6 = MPI_Function6();
        if (rank == 0) {
            cout << "6:" << time6 << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);//进程同步
        //=================================================================
    }
    MPI_Finalize();//终止
    return 0;
}

