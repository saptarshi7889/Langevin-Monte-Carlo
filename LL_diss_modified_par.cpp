//
//  main.cpp
//  Langevin dynmaics for XXZ spin chain with dissipation
//
//  Created by smajumdar on 26/03/2021.
//

#include <iostream>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<cmath>
#include<random>
#include<iomanip>
#include"Eigen/SVD"
#include<fstream>
#include"mpi.h"

#define Pi 3.1415926535897
#define  MASTER 0

using namespace std;
using namespace Eigen;


void Laplace_X(const Ref<const MatrixXd> a, Ref<MatrixXd> b, int R, int C)
{
    int i,j;
    for(i=0;i<R;i++)
    {
        for(j=0; j<C; j++)
        {
           if(j==0)
               b(i,j)=a(i,j+1)+a(i,C-1)-(2.0*a(i,j));
           else if(j==C-1)
                b(i,j)=a(i,0)+a(i,j-1)-(2.0*a(i,j));
           else
               b(i,j)=a(i,j+1)+a(i,j-1)-(2.0*a(i,j));
        }
    }
} //function to calculate the X-laplacian of a matrix a[]

void Laplace_T(const Ref<const MatrixXd> a, Ref<MatrixXd> b, int R, int C)

{
    int i,j;
    for(i=0;i<R;i++)
    {
        for(j=0;j<C;j++)
        {
           if(i==0)
               b(i,j)=a(i+1,j)+a(R-1,j)-(2.0*a(i,j));
           else if(i==R-1)
                b(i,j)=a(0,j)+a(i-1,j)-(2.0*a(i,j));
           else
               b(i,j)=a(i+1,j)+a(i-1,j)-(2.0*a(i,j));
        }
    }
} //function to calculate the Tau-laplacian of a matrix a[]

void RoughnessX(const Ref<const MatrixXd> a, double b[], int R, int C)
{
    int i,j,k,l;
    for(i=0;i<R;i++)
    {
        for(j=0;j<C;j++)
        {
            for(k=1;k<C;k++)
            {
                l=(j+k+C)%C;
                b[k-1]=b[k-1]+(pow((a(i,l)-a(i,j)),2))/(R*C);
            }
        }
    }
} // function to calculate the spatial roughness correlation function from matrix a[] with translational symmetry and periodic boundary condition

void RoughnessT(const Ref<const MatrixXd> a, double b[], int R, int C)
{
    int i,j,k,l;
    for(j=0;j<C;j++)
    {
        for(i=0;i<R;i++)
        {
            for(k=1;k<R;k++)
            {
                l=(i+k+R)%R;
                //b[k-1]=b[k-1]+(pow((a(i,j+k)-a(i,j)),2))/(C-k);
                b[k-1]=b[k-1]+(pow((a(l,j)-a(i,j)),2))/(C*R);
            }
        }
    }
} // function to calculate the imaginary time roughness correlation function from matrix a[] with translational symmetry and periodic boundary condition

void Kernel(double b[], int R, double a)
{
  int i,j;
  double sum,s1;
  for(i=1;i<=R-1;i++)
    {
     sum=0.0;
     for(j=-R/2;j<R/2;j++)
      {
       s1=pow((i+j*R),2)-(1.0/4.0);
       sum=sum+(a/(s1*Pi*Pi));
      } 
    b[i-1]=sum;
    }
}


void Diss(const Ref<const MatrixXd> a, Ref<MatrixXd> b,double c[], int R, int C)
{
    int i,j,k;
    //double ker;
    for(j=0;j<C;j++)
    {
        for(i=0;i<R;i++)
        {
            for(k=0;k<R;k++)
            {
                if(k!=i)
                {
                    b(i,j)=b(i,j)+c[abs(i-k)-1]*sin(2.0*(a(k,j)-a(i,j))); // for langevin
                    //b(i,j)=b(i,j)+c[abs(i-k)-1]*2.0*(a(k,j)-a(i,j)); // for harmonic approximation in langevin
                }
            }
        }
    }
}
// Function to calculate the dissipative matrix (- delta S_diss/delta Phi) from a particular configuration Phi(t). The argument c[] is the array of kernel.

int main(int argc, char *argv[]) {
    int   numtasks, taskid, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Get_processor_name(hostname, &len);
    printf ("Hello from task %d on %s!\n", taskid, hostname);
    if (taskid == MASTER)
       printf("MASTER: Number of MPI tasks is: %d\n",numtasks);
    
    
    double u = 1.0,K = 0.45,alpha=1.0,s=1.0; //Luttinger liquid parameters defined 5.196152422
    double m1 = u/(Pi*K);
    double m2 = 1/(Pi*u*K); //multiplicative factors with Pi etc. taken care of.
   // double m3 = alpha/(2.0*Pi*sin(Pi*s/2)); // only to be multiplied when using the function Kernel() and not Kernel1()!
    int L=20,beta=20,i,j;//time-space size
    MatrixXd Phi(beta,L),Phi2_X(beta,L),Phi2_T(beta,L),Diss_act(beta,L);
    //Defining the configuration matrix Phi, space and time laplacian of Phi and the the dissipative matrix (-delta S_diss/delta Phi). Here x is the column of matrix (j) and tau is the row of matrix (i). j goes from 0 to L-1 and i goes from 0 to beta-1.
    for(i=0;i<beta;i++)
    {
        for(j=0;j<L;j++)
        {
         Phi(i,j)=0.0;
        }
    } // initial configuration Phi=constant=0.0;
    
    // introduce noise
    random_device rd;
    mt19937 e2(rd());
    normal_distribution<> dist(0,sqrt(2.0)); //generates number from gaussian distribution of mean 0 and variance 2
    
    //creating a file for taking data output
    double t_init=0.0, t_fin=10000.0, delta_t=0.01, t, t_eq=100.0;
    int skip=(int)(t_eq/delta_t);
    char filename[100],filename1[100];
    sprintf(filename,"diss_roughness_alpha%gs%gsize%dtime%ddeltat%g_kernel3_harmonic.dat",alpha,s,L,(int)t_fin,delta_t);
    sprintf(filename1,"diss_roughness_alpha%gs%gsize%dtime%ddeltat%g_kernel3_hconfig.dat",alpha,s,L,(int)t_fin,delta_t);
    ofstream myfile,myfile1;
    myfile.open(filename, fstream::out);
    myfile1.open(filename1, fstream::out);
    
    // Langevin dynamics
    
    double B_X[L-1],B_T[beta-1]; //defining arrays for spatial and imaginary time roughness functions, respectively.
    for(i=0;i<L-1;i++)
    B_X[i]=B_T[i]=0.0;
    int c=0;//,c1=1;;//counting
    
    double Diss_Ker[beta-1];
    Kernel(Diss_Ker,beta,alpha);
    
    /*for(i=0;i<beta-1;i++)
    cout<<"kernel for "<<i+1<<" distance is "<<Diss_Ker[i]<<endl;*/
    
    //dynamics time loop
    for(t=t_init;t<=t_fin;t=t+delta_t)
    {
        for(i=0;i<beta;i++)
        {
            for(j=0;j<L;j++)
            {
                Diss_act(i,j)=0.0;
            }
            
        } // in the beginning of each loop, dissipation matrix is initialised with 0
        
        Laplace_X(Phi,Phi2_X,beta,L);
        Laplace_T(Phi,Phi2_T,beta,L);
        Diss(Phi,Diss_act,Diss_Ker,beta,L);
        //calculating the required functions of Phi
        
        for(i=0;i<beta;i++)
        {
            for(j=0;j<L;j++)
            {
                Phi(i,j)=Phi(i,j)+delta_t*(m1*Phi2_X(i,j)+m2*Phi2_T(i,j)+Diss_act(i,j))+sqrt(delta_t)*dist(e2);
            }//obtaining the new configuration
        }
        
        // calculating the roughness starting from only after a certain time and then once in 50 time iterations
        if(((int)t>0)&&(c%skip==0))
        {
            myfile1<<"The configuration at time "<<t<<"is: \n"<<Phi<<endl;
            RoughnessX(Phi,B_X,beta,L);
            RoughnessT(Phi,B_T,beta,L);
        }
        c++;
        
    }
    
    int no_of_obs=int(t_fin/t_eq);
    cout<<"number of observations "<<no_of_obs<<endl;
    for(i=0;i<L-1;i++)
    myfile<<i+1<<"\t"<<B_X[i]/no_of_obs<<"\t"<<B_T[i]/no_of_obs<<"\n"; // taking the output in the file with roughness values averaged over no. of observations
    myfile.close();
    myfile1.close();
    MPI_Finalize();
    return 0;
}

