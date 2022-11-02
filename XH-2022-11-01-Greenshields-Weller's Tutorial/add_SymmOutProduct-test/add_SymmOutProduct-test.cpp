#include <unistd.h>
#include <iostream>
#include "Vec_Tools.hpp"
#include "Vector_3.hpp"
#include "Matrix_3x3.hpp"
#include "SymmMatrix_3x3.hpp"
#include "HDF5_Reader.hpp"
#include "PostVectSolution.hpp"
#include "Tensor4_3D.hpp"

int main(int argc, char *argv[])
{
  std::cout<<"\nTo test add_SymmOutProduct(val, A, B)\n"
  <<"D_ijkl = val*(A_ij*B_kl + B_ij*A_kl):\n\n";
  sleep(2);

// 普通随机测试
  std::cout<<"Ordinary random test: 5 times\n\n";
  sleep(1);
  for(int temp{1}; temp < 6; temp++)
  {
    std::cout<<"Loop "<<temp<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 A;
    A.gen_rand(); 
    std::cout<<"matrix A =\n";
    A.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 B;
    B.gen_rand();
    std::cout<<"matrix B =\n";
    B.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    double val = ( rand() % 1000 ) * 1.0e-3 - 0.5;
    std::cout<<"val = \n"<<std::setprecision(9)<<val<<"\n";
    std::cout<<"\n";
    sleep(1);

    Tensor4_3D D;
    D.add_SymmOutProduct(val, A, B);
    std::cout<<"tensor D = \n";
    D.print_in_mat();
    std::cout<<"\n";
    sleep(1);
  }

  std::cout<<"Extreme cases test:\n\n";
  sleep(1);
// 设一个矩阵的分量值较大
  std::cout<<"(1): A's components are multiplied by a larger number: 5 times\n";
  sleep(2);
  for(int temp{1}; temp < 6; temp++)
  {
    std::cout<<"Loop "<<temp<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 A;
    A.gen_rand();
    A *= 1.0e3;
    std::cout<<"matrix A =\n";
    A.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 B;
    B.gen_rand();
    std::cout<<"matrix B =\n";
    B.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    double val = ( rand() % 1000 ) * 1.0e-3 - 0.5;
    std::cout<<"val = \n"<<std::setprecision(9)<<val<<"\n";
    std::cout<<"\n";
    sleep(1);

    Tensor4_3D D;
    D.add_SymmOutProduct(val, A, B);
    std::cout<<"tensor D = \n";
    D.print_in_mat();
    std::cout<<"\n";
    sleep(1);
  }

// 设一个矩阵的分量值较小
  std::cout<<"(2): A's components are multiplied with a smaller number: 5 times\n";
  sleep(2);
  for(int temp{1}; temp < 6; temp++)
  {
    std::cout<<"Loop "<<temp<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 A;
    A.gen_rand();
    A *= 1.0e-3;
    std::cout<<"matrix A =\n";
    A.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 B;
    B.gen_rand();
    std::cout<<"matrix B =\n";
    B.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    double val = ( rand() % 1000 ) * 1.0e-3 - 0.5;
    std::cout<<"val = \n"<<std::setprecision(9)<<val<<"\n";
    std::cout<<"\n";
    sleep(1);

    Tensor4_3D D;
    D.add_SymmOutProduct(val, A, B);
    std::cout<<"tensor D = \n";
    D.print_in_mat();
    std::cout<<"\n";
    sleep(1);
  }

// 设val较大
  std::cout<<"(3): Larger val: 5 times\n";
  sleep(2);
  for(int temp{1}; temp < 6; temp++)
  {
    std::cout<<"Loop "<<temp<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 A;
    A.gen_rand();
    std::cout<<"matrix A =\n";
    A.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 B;
    B.gen_rand();
    std::cout<<"matrix B =\n";
    B.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    double val = (( rand() % 1000 ) * 1.0e-3 - 0.5) * 1.0e5;
    std::cout<<"val = \n"<<std::setprecision(9)<<val<<"\n";
    std::cout<<"\n";
    sleep(1);

    Tensor4_3D D;
    D.add_SymmOutProduct(val, A, B);
    std::cout<<"tensor D = \n";
    D.print_in_mat();
    std::cout<<"\n";
    sleep(1);
  }

// 设val较小
  std::cout<<"(4): Smaller val: 5 times\n";
  sleep(2);
  for(int temp{1}; temp < 6; temp++)
  {
    std::cout<<"Loop "<<temp<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 A;
    A.gen_rand();
    std::cout<<"matrix A =\n";
    A.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    Matrix_3x3 B;
    B.gen_rand();
    std::cout<<"matrix B =\n";
    B.print();
    std::cout<<"\n";
    sleep(1);

    srand(time(NULL));
    double val = (( rand() % 1000 ) * 1.0e-3 - 0.5) * 1.0e-3;
    std::cout<<"val = \n"<<std::setprecision(9)<<val<<"\n";
    std::cout<<"\n";
    sleep(1);

    Tensor4_3D D;
    D.add_SymmOutProduct(val, A, B);
    std::cout<<"tensor D = \n";
    D.print_in_mat();
    std::cout<<"\n";
    sleep(1);
  }

  std::cout<<"Test finished.\n";
  return EXIT_SUCCESS;
}

// EOF