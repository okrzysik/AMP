
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>

#include "ml_include.h"

struct GlobalData {
  int N;
  double** mat;
} myData;

void createMatrix() {
  typedef double* doublePtr;
  myData.mat = new doublePtr[myData.N];
  for(int i = 0; i < myData.N; ++i) {
    myData.mat[i] = new double[myData.N];
  }//end for i
}

void computeMatrix() {
}

void freeMatrix() {
  for(int i = 0; i < myData.N; ++i) {
    delete [] (myData.mat[i]);
    myData.mat[i] = NULL;
  }//end for i
  delete [] (myData.mat);
  myData.mat = NULL;
}

int myMatVec(ML_Operator *data, int in_length, double in[], int out_length, double out[]) {
  for(int i = 0; i < out_length; ++i) {
    out[i] = 0.0;
    for(int j = 0; j < in_length; ++j) {
      out[i] += ((myData.mat[i][j])*in[j]);
    }//end for j
  }//end for i
  return 0;
}

int myGetRow(ML_Operator *data, int N_requested_rows, int requested_rows[],
    int allocated_space, int columns[], double values[], int row_lengths[]) {
  int spaceRequired = 0;
  int cnt = 0;
  for(int i = 0; i < N_requested_rows; ++i) {
    int row = requested_rows[i];
    std::vector<unsigned int> cols;
    std::vector<double> vals;

    for(int j = 0; j < myData.N; ++j) {
      if(fabs(myData.mat[row][j]) > 1.0e-15) {
        cols.push_back(j);
        vals.push_back(myData.mat[row][j]);
      }
    }//end for j

    spaceRequired += cols.size();
    if(allocated_space >= spaceRequired) {
      for(size_t j = 0; j < cols.size(); ++j) {
        columns[cnt] = cols[j];
        values[cnt] = vals[j];
        ++cnt;
      }//end for j
      row_lengths[i] = cols.size();
    } else {
      return 0;
    }
  }//end for i
  return 1;
}

int main(int argc, char *argv[])
{
  myData.N = atoi(argv[1]);
  createMatrix();
  computeMatrix();

  const int numGrids = 10;
  const int numPDEs = 1;
  const int maxIterations = 30;

  ML_set_random_seed(123456);
  ML* ml_object;
  ML_Create(&ml_object, numGrids);

  ML_Init_Amatrix(ml_object, 0, myData.N, myData.N, &myData);
  ML_Set_Amatrix_Getrow(ml_object, 0, &myGetRow, NULL, myData.N);
  ML_Set_Amatrix_Matvec(ml_object, 0, &myMatVec);
  ML_Set_MaxIterations(ml_object, maxIterations);
  ML_Set_ResidualOutputFrequency(ml_object, 1);
  ML_Set_PrintLevel(10);
  ML_Set_OutputLevel(ml_object, 10);

  ML_Aggregate* agg_object;
  ML_Aggregate_Create(&agg_object);
  agg_object->num_PDE_eqns = numPDEs;
  agg_object->nullspace_dim = 0; 
  //agg_object->nullspace_dim = numPDEs;
  ML_Aggregate_Set_MaxCoarseSize(agg_object, 128);
  ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(agg_object);

  const int nlevels = ML_Gen_MGHierarchy_UsingAggregation(ml_object, 0, ML_INCREASING, agg_object);
  std::cout<<"Number of actual levels: "<<nlevels<<std::endl;

  for(int lev = 0; lev < (nlevels - 1); ++lev) {
    ML_Gen_Smoother_SymGaussSeidel(ml_object, lev, ML_BOTH, 2, 1.0);
  }
  ML_Gen_Smoother_Amesos(ml_object, (nlevels - 1), ML_AMESOS_KLU, -1, 0.0);

  ML_Gen_Solver(ml_object, ML_MGV, 0, (nlevels-1));

  double* solArr = new double[myData.N];
  double* rhsArr = new double[myData.N];

  for(int i = 0; i < myData.N; i++) {
    solArr[i] = (static_cast<double>(rand()))/(static_cast<double>(RAND_MAX));
  }//end for i

  myMatVec(NULL, myData.N, solArr, myData.N, rhsArr);

  for(int i = 0; i < myData.N; i++) {
    solArr[i] = 0.0;
  }//end for i

  ML_Iterate(ml_object, solArr, rhsArr);

  ML_Aggregate_Destroy(&agg_object);
  ML_Destroy(&ml_object);

  delete [] solArr;
  delete [] rhsArr;
  freeMatrix();
}  



