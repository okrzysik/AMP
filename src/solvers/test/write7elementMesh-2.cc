
#include <cstdlib>
#include <cstdio>
#include <iostream>

int main(int argc, char**argv) {

  if(argc == 1) {
    std::cout<<"Usage: argv[0] output_file "<<std::endl;
    exit(0);
  }

  FILE* fp = fopen(argv[1], "w");

  fprintf(fp, "Mesh { \n");

  fprintf(fp, "NumberOfNodes = 16 \n");

  fprintf(fp, "Point0 = 0.0, 0.0, 0.0 \n");
  fprintf(fp, "Point1 = 1.0, 0.0, 0.0 \n");
  fprintf(fp, "Point2 = 1.0, 1.0, 0.0 \n");
  fprintf(fp, "Point3 = 0.0, 1.0, 0.0 \n");
  fprintf(fp, "Point4 = 0.0, 0.0, 1.0 \n");
  fprintf(fp, "Point5 = 1.0, 0.0, 1.0 \n");
  fprintf(fp, "Point6 = 1.0, 1.0, 1.0 \n");
  fprintf(fp, "Point7 = 0.0, 1.0, 1.0 \n");
  fprintf(fp, "Point8 = 0.249, 0.342, 0.192 \n");
  fprintf(fp, "Point9 = 0.826, 0.288, 0.288 \n");
  fprintf(fp, "Point10 = 0.850, 0.649, 0.263 \n");
  fprintf(fp, "Point11 = 0.273, 0.750, 0.230 \n");
  fprintf(fp, "Point12 = 0.320, 0.186, 0.643 \n");
  fprintf(fp, "Point13 = 0.677, 0.305, 0.683 \n");
  fprintf(fp, "Point14 = 0.788, 0.693, 0.644 \n");
  fprintf(fp, "Point15 = 0.165, 0.745, 0.702 \n");

  fprintf(fp, "\n");

  fprintf(fp, "NumberOfElements = 7 \n");

  fprintf(fp, "Elem0 = 0, 1, 2, 3, 8, 9, 10, 11 \n");
  fprintf(fp, "Elem1 = 12, 13, 14, 15, 4, 5, 6, 7 \n");
  fprintf(fp, "Elem2 = 0, 8, 11, 3, 4, 12, 15, 7 \n");
  fprintf(fp, "Elem3 = 9, 1, 2, 10, 13, 5, 6, 14 \n");
  fprintf(fp, "Elem4 = 0, 1, 9, 8, 4, 5, 13, 12 \n");
  fprintf(fp, "Elem5 = 11, 10, 2, 3, 15, 14, 6, 7 \n");
  fprintf(fp, "Elem6 = 8, 9, 10, 11, 12, 13, 14, 15 \n");

  fprintf(fp,"\n");

  fprintf(fp, "NumberOfBoundaryNodeIds = 8 \n\n");

  for(int i = 0; i < 8; i++) {
    fprintf(fp, "NumberOfBoundaryNodes%d = 1 \n", (i+1));
    fprintf(fp, "BoundaryNodeId%d = %d",(i + 1), i);
    fprintf(fp,"\n\n");
  }//end for i

  fprintf(fp, "NumberOfBoundarySideIds = 0 \n\n");

  fprintf(fp, "} \n");

  fclose(fp);

}



