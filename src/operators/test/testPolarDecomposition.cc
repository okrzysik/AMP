#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "operators/mechanics/UpdatedLagrangianUtils.h" 

#include <iostream>


void myTest(AMP::UnitTest *ut) {

  double A[3][3] = { {-2, 2, -2},
    {-1, 1, 3},
    {2, 0, -1} };

  double R[3][3];
  double U[3][3];
  double F[3][3];
  double R2[3][3], U2[3][3], F2[3][3];

  AMP::Operator::polarDecomposeRU(A, R, U);
  AMP::Operator::matMatMultiply(R, U, F);
  AMP::Operator::polarDecompositionFeqRU_Simo(A, R2, U2);
  AMP::Operator::matMatMultiply(R2, U2, F2);

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      std::cout<<"U["<<i<<"]["<<j<<"]="<<U[i][j]<<" U2["<<i<<"]["<<j<<"]="<<U2[i][j]<<std::endl;
      std::cout<<"R["<<i<<"]["<<j<<"]="<<R[i][j]<<" R2["<<i<<"]["<<j<<"]="<<R2[i][j]<<std::endl;
    }
  }

  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      std::cout<<"F["<<i<<"]["<<j<<"]="<<F[i][j]<<" F2["<<i<<"]["<<j<<"]="<<F2[i][j]<<std::endl;

  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(!AMP::Operator::softEquals(U[i][j] , U[j][i])) {
        ut->failure("U is not symmetric.");
      }
    }//end for j
  }//end for i

  double eVal[3];
  AMP::Operator::eigenValues(U, eVal);
  for(int i = 0; i < 3; i++) {
    if(eVal[i] <= 0) {
      ut->failure("U is not positive definite.");
    }
  }

  double eValSum = eVal[0] + eVal[1] + eVal[2];
  double Utrace = U[0][0] + U[1][1] + U[2][2];
  if(!AMP::Operator::softEquals(eValSum, Utrace)) {
    ut->failure("Trace(U) != sum of EigenValues.");
  }

  double eValProd = eVal[0] * eVal[1] * eVal[2]; 
  double detU = AMP::Operator::matDeterminant(U);
  if(!AMP::Operator::softEquals(eValProd, detU)) {
    ut->failure("det(U) != product of EigenValues.");
  }

  double Rtran[3][3];
  AMP::Operator::matTranspose(R, Rtran);

  double RtranR[3][3];
  AMP::Operator::matMatMultiply(Rtran, R, RtranR);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(i == j) {
        if(!AMP::Operator::softEquals(RtranR[i][j] , 1.0)) {
          ut->failure("R is not orthogonal.");
        }
      } else {
        if(!AMP::Operator::softEquals(RtranR[i][j] , 0.0)) {
          ut->failure("R is not orthogonal.");
        }
      }
    }//end for j
  }//end for i

  double RU[3][3];
  AMP::Operator::matMatMultiply(R, U, RU);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      if(!AMP::Operator::softEquals(RU[i][j] , A[i][j])) {
        ut->failure("RU != A.");
      }
    }//end for j
  }//end for i

  ut->passes("testPolarDecomposition");
}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
        myTest(&ut);
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
        ut.failure("ERROR: While testing");
    } catch( ... ) {
        std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
        ut.failure("ERROR: While testing");
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   



