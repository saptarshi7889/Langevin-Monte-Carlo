//
//  main.cpp
//  Langevin dynmaics for XXZ spin chain with dissipation
//
//  Created by smajumdar on 26/03/2021.
//

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<cmath>
#include<random>
#include<iomanip>
#include"Eigen/SVD"
#include<fstream>
#include<string>
#include<cstring>
#include<unordered_map>




#define Pi 3.1415926535897

using namespace std;
using namespace Eigen;

unordered_map<string,string> process_arguments(int argc, char** argv)
{
    std::string delimiter = "=";
    unordered_map<string,string> dict;
    for(int i = 1; i < argc; i++)
    {
        string str = string(argv[i]);
        size_t pos = str.find(delimiter);
        string name = str.substr(0,pos);
        string value = str.substr(pos+1);
        dict.insert(make_pair(name,value));
    }
    return dict;
}


void Laplace_X(const Ref<const MatrixXd> a, Ref<MatrixXd> b, int R, int C)
{
    int i,j;
    for(i=0;i<R;++i)
    {
        for(j=0; j<C; ++j)
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
    for(i=0;i<R;++i)
    {
        for(j=0;j<C;++j)
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
    for(i=0;i<R;++i)
    {
        for(j=0;j<C;++j)
        {
            for(k=1;k<C;++k)
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
    for(j=0;j<C;++j)
    {
        for(i=0;i<R;++i)
        {
            for(k=1;k<R;++k)
            {
                l=(i+k+R)%R;
                b[k-1]=b[k-1]+(pow((a(l,j)-a(i,j)),2))/(C*R);
            }
        }
    }
} // function to calculate the imaginary time roughness correlation function from matrix a[] with translational symmetry and periodic boundary condition

void OrderParameter(const Ref<const MatrixXd> a, double &b, int R, int C)
{
    int i,j;
    for(j=0;j<C;++j)
    {
        for(i=0;i<R;++i)
        {
            //b=b+(cos(a(i,j))/(C*R));
              b=b+cos(a(i,j));
        }
    }
}

double CenterOfMass(const Ref<const MatrixXd> a, int R, int C)
{
    int i,j;
    double sum=0.0;
    for(j=0;j<C;++j)
    {
        for(i=0;i<R;++i)
        {
            sum=sum+a(i,j);
        }
    }
    sum = sum/(R*C);
    return sum;
}

void Kernel(double b[],int R, double a,double c)
{
    int i,j;
    a = a/(Pi*Pi);
    double sum,s1;
    for(i=1;i<=R-1;++i)
    {
       sum=0.0;
       for(j=-R/2;j<R/2;++j)
        {
           s1=beta((double)(i+j*R)-(c/2),(c+1));
           sum=sum+s1;
        }
      b[i-1]=a*sum;
    }
} //exact kernel in real space

void Diss(const Ref<const MatrixXd> a, Ref<MatrixXd> b,double c[], int R, int C)
{
    int i,j,k;
    b=ArrayXXd::Zero(R,C);
    for(j=0;j<C;++j)
    {
        for(i=0;i<R;++i)
        {
            for(k=0;k<R;++k)
            {
                if(k!=i)
                {
                  b(i,j)=b(i,j)+c[abs(i-k)-1]*sin(2.0*(a(k,j)-a(i,j)));
                  //b(i,j)=b(i,j)+c[abs(i-k)-1]*2.0*(a(k,j)-a(i,j));
                }
            }
        }
    }
}
// Function to calculate the dissipative matrix (- delta S_diss/delta Phi) from a particular configuration Phi(t). The argument c[] is the array of kernel.

int main(int argc, char *argv[]) {
    
    int simulation_id = 0;
    if (argc >1)
    {
        cout<<"Trying to load the file \n";
        unordered_map<string,string> args = process_arguments(argc,argv);
            for (std::pair<string, string> element : args)
            {
             if (element.first == "id")
               {
                   simulation_id=stoi(element.second);
               }
            }
    }
    double u = 1.0,K = 1,s=1.5,alpha=1.0; //Luttinger liquid parameters defined
    double m1 = u/(Pi*K);
    double m2 = 1/(Pi*u*K); //multiplicative factors with Pi etc. taken care of.
    int L=16,Beta=16,i,j;//time-space size
    MatrixXd Phi(Beta,L),Phi2_X(Beta,L),Phi2_T(Beta,L),Phi_star(Beta,L),Diss_act(Beta,L),F1(Beta,L),F2(Beta,L),Noise(Beta,L),Phi_K(Beta,L);
    //Defining the configuration matrix Phi, space and time laplacian of Phi and the the dissipative matrix (-delta S_diss/delta Phi). Here x is the column of matrix (j) and tau is the row of matrix (i). j goes from 0 to L-1 and i goes from 0 to beta-1.
    Phi=ArrayXXd::Zero(Beta,L);// initial configuration Phi=constant=0.0;
    
    // introduce noise
    random_device rd;
    mt19937 e2(rd());
    normal_distribution<> dist(0.0,1.0); //generates number from gaussian distribution of mean 0 and variance 1
    
    double t_init=0.0, t_fin=16.0, delta_t=0.05, t, t_eq2=16.0;
    double t_eq1 = t_eq2*4.0;
    int skip=(int)(t_eq2/delta_t);
    char filename[200],filename1[200],filename2[200],filename3[200],filename4[200];
    
    sprintf(filename1,"time_alpha%gK%gs%gsize%dtime%ddeltat%g_%d_langevin_image.dat",alpha,K,s,L,(int)t_fin,delta_t,(simulation_id+1));
    sprintf(filename2,"order_parameter_alpha%gK%gs%gsize%dtime%ddeltat%g_%d_langevin_image.dat",alpha,K,s,L,(int)t_fin,delta_t,(simulation_id+1));
    sprintf(filename3,"diss_roughness_alpha%gK%gs%gsize%dtime%ddeltat%g_%d_langevin_image.dat",alpha,K,s,L,(int)t_fin,delta_t,(simulation_id+1));
    sprintf(filename4,"num_of_config_alpha%gK%gs%gsize%dtime%ddeltat%g_%d_langevin_image.dat",alpha,K,s,L,(int)t_fin,delta_t,(simulation_id+1));
    ofstream myfile,myfile1,myfile2,myfile3,myfile4;
    myfile1.open(filename1, fstream::out);
    
    double OP=0.0,Phi_COM, B_X[L-1],B_T[Beta-1]; //defining arrays for spatial and imaginary time roughness functions, respectively.
   for(i=0;i<L-1;++i)
   B_X[i]=B_T[i]=0.0;
    int c=0,c1=0;//counting
    
    //kernel initiated//
    double Diss_Ker[Beta-1];
    Kernel(Diss_Ker,Beta,alpha,s);
    char filename_kernel[200];
    sprintf(filename_kernel,"kernel_alpha%gK%gs%gsize%d_image_simul.dat",alpha,K,s,L); //time_count
    ofstream myfile_kernel;
    myfile_kernel.open(filename_kernel, fstream::out);
     for(i=0;i<Beta-1;i++)
     myfile_kernel<<Diss_Ker[i]<<endl;
    myfile_kernel.close();
    
    cout<<"Size: "<<L<<" alpha: "<<alpha<<" K: "<<K<<" delta t :"<<delta_t<<endl;
        for(t=t_init;t<=t_fin;t=t+delta_t)
        {
            for(i=0;i<Beta;++i)
            {
                for(j=0;j<L;++j)
                {
                    Noise(i,j)=dist(e2);
                }
            }
        //dynamics time loop
            Phi_K=ArrayXXd::Zero(Beta,L);
            
            Laplace_X(Phi,Phi2_X,Beta,L);
            Laplace_T(Phi,Phi2_T,Beta,L);
            Diss(Phi,Diss_act,Diss_Ker,Beta,L);
            //calculating the required functions of Phi
            F1=m1*Phi2_X+m2*Phi2_T+Diss_act;//K1 defined
            
        Phi_K=Phi+(delta_t*F1)+(sqrt(2.0*delta_t)*Noise);
           Laplace_X(Phi_K,Phi2_X,Beta,L);
           Laplace_T(Phi_K,Phi2_T,Beta,L);
           Diss(Phi_K,Diss_act,Diss_Ker,Beta,L);
           F2=m1*Phi2_X+m2*Phi2_T+Diss_act;//K2 defined
                     
            
            
            Phi=Phi+((delta_t/2.0)*(F1+F2))+(sqrt(2.0*delta_t)*Noise);//obtaining the new configuration
            
            
            // calculating the roughness starting from only after a certain time and then once in 50 time iterations
            if((t>=t_eq1)&&(c1%skip==0))
            {
                myfile1<<"the time is "<<t<<endl;
                Phi_COM = CenterOfMass(Phi, Beta, L);
                for(i=0;i<Beta;++i)
                  {
                   for(j=0;j<L;++j)
                      {
                          Phi_star(i,j) = Phi(i,j)-Phi_COM;
                      }
                  }
                RoughnessX(Phi_star,B_X,Beta,L);
                RoughnessT(Phi_star,B_T,Beta,L);
                OrderParameter(Phi_star,OP,Beta,L);
                ++c;
                sprintf(filename,"configuration_alpha%gK%gs%gsize%dtime%ddeltat%gcore%d_%d_langevin_image.dat",alpha,K,s,L,(int)t_fin,delta_t,simulation_id+1,c);
                myfile.open(filename,fstream::out);
                myfile<<Phi_star<<endl;
                myfile.close(); //configuration saved
                myfile2.open(filename2,fstream::out);
                myfile2<<alpha<<"\t"<<OP/(c*Beta*L)<<endl; // OP saved
                myfile2.close();
                myfile3.open(filename3,fstream::out);
                for(i=0;i<L-1;++i)
                 {
                myfile3<<i+1<<"\t"<<B_X[i]/c<<"\t"<<B_T[i]/c<<"\n"; // taking the output in the file with roughness values averaged over no. of observations
                 }
                myfile3.close(); // B(r) saved
                myfile4.open(filename4,fstream::out);
                myfile4<<c<<endl;
                myfile4.close();
            }
            ++c1;
        }
        
        //int no_of_obs=int(t_fin/t_eq);
        cout<<"number of observations "<<c<<endl;
        myfile1.close();
    return 0;
}


