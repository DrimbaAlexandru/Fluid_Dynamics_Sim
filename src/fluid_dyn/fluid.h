#include "../utils/utils.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>

#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}

class fluid_space
{
private:
    float factor;
    int phase;
    int N;
    long size;
    float *u, *v, *u_prev, *v_prev;
    float *dens, *dens_prev;
    float _dt,_visc,_diff;
    float *a, *b;

    int IX(int i, int j)
    {
        return (i)+(N+2)*(j);
    }

    void initialize()
    {
        u=Utils::allocate_vector<float>(size);
        Utils::fill_vector<float>(size,u,0);

        v=Utils::allocate_vector<float>(size);
        Utils::fill_vector<float>(size,v,0);

        u_prev=Utils::allocate_vector<float>(size);
        Utils::fill_vector<float>(size,u_prev,0);

        v_prev=Utils::allocate_vector<float>(size);
        Utils::fill_vector<float>(size,v_prev,0);

        dens=Utils::allocate_vector<float>(size);
        Utils::fill_vector<float>(size,dens,0);

        dens_prev=Utils::allocate_vector<float>(size);
        Utils::fill_vector<float>(size,dens_prev,0);

        a=Utils::allocate_vector<float>(size);
        Utils::fill_vector<float>(size,a,0);

        b=Utils::allocate_vector<float>(size);
        Utils::fill_vector<float>(size,b,0);


         srand (time(NULL));
    }

    void diffuse ( int N, int b, float * x, float * x0, float diff)
    {
        int i, j, k;
        float a=_dt*diff*N*N;
        for ( k=0 ; k<20 ; k++ )
        {
            for ( i=1 ; i<=N ; i++ )
            {
                for ( j=1 ; j<=N ; j++ )
                {
                    x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+
                    x[IX(i,j-1)]+x[IX(i,j+1)]))/(1+4*a);
                }
            }
            set_bnd ( N, b, x );
        }
    }

    void advect ( int N, int b, float * d, float * d0, float * u, float * v)
    {
        int i, j, i0, j0, i1, j1;
        float x, y, s0, t0, s1, t1, dt0;
        dt0 = _dt*N;
        for ( i=1 ; i<=N ; i++ )
        {
            for ( j=1 ; j<=N ; j++ )
            {
                x = i-dt0*u[IX(i,j)];
                y = j-dt0*v[IX(i,j)];
                if (x<0.5) x=0.5;
                if (x>N+0.5) x=N+ 0.5;
                i0=(int)x;
                i1=i0+1;
                if (y<0.5) y=0.5;
                if (y>N+0.5) y=N+ 0.5;
                j0=(int)y;
                j1=j0+1;
                s1 = x-i0;
                s0 = 1-s1;
                t1 = y-j0;
                t0 = 1-t1;
                d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
            }
        }
        set_bnd ( N, b, d );
    }

    void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff)
    {
        add_source ( N, x, x0);
        SWAP ( x0, x );
        diffuse ( N, 0, x, x0, diff );
        SWAP ( x0, x );
        advect ( N, 0, x, x0, u, v );
    }

    void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc )
    {
        add_source ( N, u, u0);
        add_source ( N, v, v0);
        SWAP ( u0, u );
        diffuse ( N, 1, u, u0, visc);
        SWAP ( v0, v );
        diffuse ( N, 2, v, v0, visc);
        project ( N, u, v, u0, v0 );
        SWAP ( u0, u );
        SWAP ( v0, v );
        advect ( N, 1, u, u0, u0, v0);
        advect ( N, 2, v, v0, u0, v0);
        project ( N, u, v, u0, v0 );
    }

    void project ( int N, float * u, float * v, float * p, float * div )
    {
        int i, j, k;
        float h;
        h = 1.0/N;
        for ( i=1 ; i<=N ; i++ )
        {
            for ( j=1 ; j<=N ; j++ )
            {
                div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+
                v[IX(i,j+1)]-v[IX(i,j-1)]);
                p[IX(i,j)] = 0;
            }
        }
        set_bnd ( N, 0, div );
        set_bnd ( N, 0, p );
        for ( k=0 ; k<20 ; k++ )
        {
            for ( i=1 ; i<=N ; i++ )
            {
                for ( j=1 ; j<=N ; j++ )
                {
                    p[IX(i,j)] = (div[IX(i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+
                                  p[IX(i,j-1)]+p[IX(i,j+1)])/4;
                }
            }
            set_bnd ( N, 0, p );
        }
        for ( i=1 ; i<=N ; i++ )
        {
            for ( j=1 ; j<=N ; j++ )
            {
                u[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
                v[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
            }
        }
        set_bnd ( N, 1, u );
        set_bnd ( N, 2, v );
    }

    void set_bnd ( int N, int b, float * x )
    {
        int i;
        for ( i=1 ; i<=N ; i++ )
        {
            x[IX(0,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
            x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
            x[IX(i,0 )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
            x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
        }
        x[IX(0,0 )] = 0.5*(x[IX(1,0 )]+x[IX(0,1)]);
        x[IX(0,N+1)] = 0.5*(x[IX(1,N+1)]+x[IX(0,N )]);
        x[IX(N+1,0 )] = 0.5*(x[IX(N,0 )]+x[IX(N+1,1)]);
        x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N )]);
    }

    void add_velocity(float *_u, float *_v)
    {
        int i;
        for(i=0;i<size;i++)
        {
            u[i]+=_u[i];
            v[i]+=_v[i];
            u_prev[i]+=_u[i];
            v_prev[i]+=_v[i];
        }
    }
    void add_source ( int N, float *x, float *s)
    {
        int i;
        for ( i=0 ; i<size ; i++ )
            x[i] += _dt*s[i];
    }


public:

    fluid_space( int N, float time_step, float viscosity, float diffusion_rate)
    {
        this->N=N;
        size=(N+2)*(N+2);
        initialize();
        _dt=time_step;
        _visc=viscosity;
        _diff=diffusion_rate;
        factor=1;
        phase=0;
    }

    void add_density(float** src)
    {
        int dim1=N,dim2=N;

        Utils::border_matrix<float>(src,dim1,dim2,0);
        float* vect=Utils::liniarize_matrix<float>(src,dim1,dim2);
        int i;

        for(i=0;i<dim1*dim2;i++)
        {
            dens[i]+=vect[i];
            dens_prev[i]+=vect[i];
        }

        Utils::deallocate_vector<float>(vect);
        Utils::deallocate_matrix<float>(src,dim1);
    }

    float** get_density()
    {
        float** ret=Utils::allocate_matrix<float>(N,N);
        float** aux=Utils::deliniarize_matrix<float>(dens,N+2,N+2);
        int i,j;
        for(i=0;i<N;i++)
            for(j=0;j<N;j++)
                ret[i][j]=aux[i+1][j+1];
        Utils::deallocate_matrix<float>(aux,N+2);
        return ret;
    }

    void step()
    {
        vel_step ( N, u, v, u_prev, v_prev, _visc );
        dens_step ( N, dens, dens_prev, u, v, _diff);

        int i;
        for(i=N;i<size-N;i++)
        {
            a[i]=0;
            b[i]=-2;
            //a[i]+=(float)(rand()%32000)/16000-1;
            //b[i]+=(float)(rand()%32000)/16000-1;
        }

        add_velocity(a,b);

        /*phase++;
        if(phase<20)
            _diff*=1.1;
        else _diff =0.0001;*/
    }

    ~fluid_space()
    {
        delete u;
        delete v;
        delete u_prev;
        delete v_prev;
        delete dens;
        delete dens_prev;
        delete a;
        delete b;
    }
};
