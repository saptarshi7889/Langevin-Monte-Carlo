//
//  main.cpp
//  Langevin dynmaics for XXZ spin chain with dissipation
//
//  Created by smajumdar on 26/03/2021.
//

#include "kernel.hpp"

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




constexpr double Pi = langevin::PI;

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
              b=b+cos(2.0*a(i,j));
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

void Kernel(double b[], int R, double alpha, double s)
{
    langevin::kernel(b, R, alpha, s);
}

void Diss_comm(const Ref<const MatrixXd> a, Ref<MatrixXd> b,const double c[], int R, int C, double g1)
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
                  b(i,j)=b(i,j)+c[abs(i-k)-1]*(sin(2.0*(a(k,j)-a(i,j)))-sin(2.0*(a(k,j)+a(i,j))));
                  //b(i,j)=b(i,j)+c[abs(i-k)-1]*2.0*(a(k,j)-a(i,j));
                }
            }
        b(i,j)=b(i,j)-g1*(sin(4.0*a(i,j)));      
        }
    }
}
// Function to calculate the dissipative matrix (- delta S_diss/delta Phi) from a particular configuration Phi(t). The argument c[] is the array of kernel.

#include "simulation.hpp"

int main(int argc, char* argv[])
{
    const auto force = [](const MatrixXd& phi, MatrixXd& result,
                          const vector<double>& coefficients, const langevin::Parameters& p) {
        Diss_comm(phi, result, coefficients.data(), p.Beta, p.L, langevin::LOCAL_FORCE_FACTOR * p.g / (Pi * Pi));
    };
    return langevin::run(argc, argv, true, false, force);
}
