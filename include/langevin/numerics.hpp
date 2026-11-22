#pragma once

//
//  main.cpp
//  Langevin dynmaics for XXZ spin chain with dissipation
//
//  Created by smajumdar on 26/03/2021.
//

#include <Eigen/Core>
#include <cmath>

namespace langevin::detail {
using Eigen::MatrixXd;
using Eigen::Ref;
inline constexpr double ORDER_HARMONIC = 2.0;
inline constexpr double BASE_HARMONIC = 1.0;

inline void Laplace_X(const Ref<const MatrixXd> a, Ref<MatrixXd> b, int R, int C)
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

inline void Laplace_T(const Ref<const MatrixXd> a, Ref<MatrixXd> b, int R, int C)

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

inline void RoughnessX(const Ref<const MatrixXd> a, double b[], int R, int C)
{
    int i,j,k,l;
    for(i=0;i<R;++i)
    {
        for(j=0;j<C;++j)
        {
            for(k=1;k<C;++k)
            {
                l=(j+k+C)%C;
                b[k-1]=b[k-1]+(std::pow((a(i,l)-a(i,j)),2))/(R*C);
            }
        }
    }
} // function to calculate the spatial roughness correlation function from matrix a[] with translational symmetry and periodic boundary condition

inline void RoughnessT(const Ref<const MatrixXd> a, double b[], int R, int C)
{
    int i,j,k,l;
    for(j=0;j<C;++j)
    {
        for(i=0;i<R;++i)
        {
            for(k=1;k<R;++k)
            {
                l=(i+k+R)%R;
                b[k-1]=b[k-1]+(std::pow((a(l,j)-a(i,j)),2))/(C*R);
            }
        }
    }
} // function to calculate the imaginary time roughness correlation function from matrix a[] with translational symmetry and periodic boundary condition

inline void OrderParameter(const Ref<const MatrixXd> a, double &b, int R, int C, bool commensurate)
{
    int i,j;
    for(j=0;j<C;++j)
    {
        for(i=0;i<R;++i)
        {
            //b=b+(cos(a(i,j))/(C*R));
              b=b+std::cos((commensurate ? ORDER_HARMONIC : BASE_HARMONIC)*a(i,j));
        }
    }
}

inline void Diss(const Ref<const MatrixXd> a, Ref<MatrixXd> b,const double c[], int R, int C)
{
    int i,j,k;
    b=Eigen::ArrayXXd::Zero(R,C);
    for(j=0;j<C;++j)
    {
        for(i=0;i<R;++i)
        {
            for(k=0;k<R;++k)
            {
                if(k!=i)
                {
                  b(i,j)=b(i,j)+c[std::abs(i-k)-1]*std::sin(2.0*(a(k,j)-a(i,j)));
                  //b(i,j)=b(i,j)+c[abs(i-k)-1]*2.0*(a(k,j)-a(i,j));
                }
            }
        }
    }
}
// Function to calculate the dissipative matrix (- delta S_diss/delta Phi) from a particular configuration Phi(t). The argument c[] is the array of kernel.

inline void Diss_comm(const Ref<const MatrixXd> a, Ref<MatrixXd> b,const double c[], int R, int C, double g1)
{
    int i,j,k;
    b=Eigen::ArrayXXd::Zero(R,C);
    for(j=0;j<C;++j)
    {
        for(i=0;i<R;++i)
        {
            for(k=0;k<R;++k)
            {
                if(k!=i)
                {
                  b(i,j)=b(i,j)+c[std::abs(i-k)-1]*(std::sin(2.0*(a(k,j)-a(i,j)))-std::sin(2.0*(a(k,j)+a(i,j))));
                  //b(i,j)=b(i,j)+c[abs(i-k)-1]*2.0*(a(k,j)-a(i,j));
                }
            }
        b(i,j)=b(i,j)-g1*(std::sin(4.0*a(i,j)));      
        }
    }
}
// Function to calculate the dissipative matrix (- delta S_diss/delta Phi) from a particular configuration Phi(t). The argument c[] is the array of kernel.

} // namespace langevin::detail
