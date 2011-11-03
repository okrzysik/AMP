
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

#define __PI__ 3.1415926535897932

int main(int argc, char** argv) {

  if(argc < 10) {
    std::cout<<"Usage: exe le me ne pe a b c r output_file"<<std::endl;
    exit(0);
  }

  int le = atoi(argv[1]);
  int me = atoi(argv[2]);
  int ne = atoi(argv[3]);
  int pe = atoi(argv[4]);

  double a = atof(argv[5]); //X-dimension
  double b = atof(argv[6]); //Z-dimension
  double c = atof(argv[7]); //Y-dimension
  double r = atof(argv[8]);

  std::vector<double> lYarr(le + 1);
  std::vector<double> mXarr(me + 1);
  std::vector<double> nXarr(ne + 1);
  std::vector<double> nZarr(ne + 1);
  std::vector<double> pZarr(pe + 1);
  std::vector<double> rMxArr(me + 1);
  std::vector<double> rMzArr(me + 1);
  std::vector<double> rPxArr(pe + 1);
  std::vector<double> rPzArr(pe + 1);
  std::vector<double> dNxArr(ne + 1);
  std::vector<double> dNzArr(ne + 1);

  for(int li = 0; li <= le; li++) {
    lYarr[li] = static_cast<double>(li)*c/static_cast<double>(le);
  }//end for li

  for(int mi = 0; mi <= me; mi++) {
    mXarr[mi] = static_cast<double>(mi)*a/static_cast<double>(me);
    double th = (__PI__/4.0) + (static_cast<double>(mi)*(__PI__)/(4.0*static_cast<double>(me)));
    rMxArr[mi] = r*cos(th);
    rMzArr[mi] = r*sin(th);
  }//end for mi

  double dX0 = r*cos(__PI__/4.0);
  double dZ0 = r*sin(__PI__/4.0);
  double dX1 = a;
  double dZ1 = b;
  for(int ni = 0; ni <= ne; ni++) {
    nXarr[ni] = r  + (static_cast<double>(ni)*(a - r)/static_cast<double>(ne));
    nZarr[ni] = r  + (static_cast<double>(ni)*(b - r)/static_cast<double>(ne));
    dNxArr[ni] = dX0 + ((dX1 - dX0)*static_cast<double>(ni)/static_cast<double>(ne));
    dNzArr[ni] = dZ0 + ((dZ1 - dZ0)*static_cast<double>(ni)/static_cast<double>(ne));
  }//end for ni

  for(int pi = 0; pi <= pe; pi++) {
    pZarr[pi] = static_cast<double>(pi)*b/static_cast<double>(pe);
    double th = static_cast<double>(pi)*(__PI__)/(4.0*static_cast<double>(pe));
    rPxArr[pi] = r*cos(th);
    rPzArr[pi] = r*sin(th);
  }//end for pi

  FILE* fp = fopen(argv[9], "w");

  fprintf(fp, "Mesh { \n");

  int numPts = 4*(le + 1)*(ne + 1)*(me + pe);
  fprintf(fp, "NumberOfNodes = %d \n", numPts);

  int nodeCnt = 0;

  fprintf(fp, "\n");

  int numElem = 4*le*ne*(me + pe);
  fprintf(fp, "NumberOfElements = %d \n", numElem);

  int elemCnt = 0;

  fprintf(fp,"\n");

  fprintf(fp, "NumberOfBoundaryNodeIds = 2 \n\n");

  fprintf(fp, "NumberOfBoundaryNodes1 = %d \n", ( ((2*me) + 1)*(le + 1) ) );
  fprintf(fp, "BoundaryNodeId1 = ");

  fprintf(fp,"\n\n");

  fprintf(fp, "NumberOfBoundaryNodes2 = %d \n", ( ((2*me) + 1)*(le + 1) ) );
  fprintf(fp, "BoundaryNodeId2 = ");

  fprintf(fp,"\n\n");

  fprintf(fp, "NumberOfBoundarySideIds = 0 \n\n");

  fprintf(fp, "} \n");

  fclose(fp);

}



