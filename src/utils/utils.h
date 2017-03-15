#include <windows.h>
#include <vector>
#include <iostream>
#pragma once

class Utils
{

public:

    template <typename T>
    static T* allocate_vector(int N)
    {
        T* ret=new T[N];
        return ret;
    }

    template <typename T>
    static void fill_vector(int N, T* vect, T fill_value)
    {
        int i;
        for(i=0;i<N;i++)
            vect[i]=fill_value;
    }

    template <typename T>
    static void deallocate_vector(T* vect)
    {
        delete vect;
    }

    template <typename T>
    static T** allocate_matrix(int N, int M)
    {
        T** ret=allocate_vector<T*>(N);
        int i;
        for(i=0;i<N;i++)
            ret[i]=allocate_vector<T>(M);
        return ret;
    }

    template <typename T>
    static void fill_matrix(int N, int M, T** mat, T fill_value)
    {
        int i,j;
        for(i=0;i<N;i++)
            for(j=0;j<M;j++)
                mat[i][j]=fill_value;
    }

    template <typename T>
    static void deallocate_matrix(T** mat, int N)
    {
        int i;
        for(i=0;i<N;i++)
            deallocate_vector<T>(mat[i]);
        deallocate_vector<T*>(mat);
    }

    template <typename T>
    static T* liniarize_matrix(T** mat, int N, int M)
    {
        T* vect=allocate_vector<T>(N*M);
        int i,j;
        for(i=0;i<N;i++)
            for(j=0;j<M;j++)
                vect[i*M+j]=mat[i][j];
        return vect;
    }

    template <typename T>
    static T** deliniarize_matrix(T* vect, int N, int M)
    {
        T** ret=allocate_matrix<T>(N,M);
        int i,j;
        for(i=0;i<N;i++)
            for(j=0;j<M;j++)
                ret[i][j]=vect[i*M+j];
        return ret;
    }

    template <typename T>
    static void border_matrix(T** &mat,int &N, int &M, T border_value)
    {
        N=N+2;M=M+2;
        T** ret=allocate_matrix<T>(N,M);
        int i,j;
        for(i=0;i<N;i++)
            for(j=0;j<M;j++)
                if(i==0 || i== N-1 || j==0 || j == M-1)
                    ret[i][j]=border_value;
                else
                    ret[i][j]=mat[i-1][j-1];

        //deallocate_matrix<T>(mat,N-2);
        mat=ret;
    }

    template <typename T>
    static void afis_matrice(T** mat, int N, int M)
    {
        int i,j;
        for(i=0;i<N;i++)
            {for(j=0;j<M;j++)
                std::cout<<mat[i][j]<<" ";
            std::cout<<std::endl;
            std::cout<<"--------------------------";
            std::cout<<std::endl;
            }
        std::cout<<std::endl;
        std::cout<<"++++++++++++++++++++++++++";
        std::cout<<std::endl;
    }
};


