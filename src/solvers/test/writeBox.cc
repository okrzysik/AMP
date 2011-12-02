
#include <cstdlib>
#include <cstdio>
#include <iostream>

#define __NODE__(xi, yi, zi, nx, ny) ((xi) + ((yi)*(nx)) + ((zi)*(nx)*(ny)))

int main(int argc, char**argv) {

  if(argc < 8) {
    std::cout<<"Usage: argv[0] nx ny nz Lx Ly Lz output_file"<<std::endl;
    exit(0);
  }

  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int nz = atoi(argv[3]);

  double Lx = atof(argv[4]);
  double Ly = atof(argv[5]);
  double Lz = atof(argv[6]);

  double hx = Lx/(static_cast<double>(nx - 1));
  double hy = Ly/(static_cast<double>(ny - 1));
  double hz = Lz/(static_cast<double>(nz - 1));

  FILE* fp = fopen(argv[7], "w");

  fprintf(fp, "Mesh { \n");

  fprintf(fp, "NumberOfNodes = %d \n", (nx*ny*nz));

  for(int zi = 0, pi = 0; zi < nz; zi++) {
    for(int yi = 0; yi < ny; yi++) {
      for(int xi = 0; xi < nx; xi++, pi++) {
        double p[3];
        p[0] = (static_cast<double>(xi))*hx;
        p[1] = (static_cast<double>(yi))*hy;
        p[2] = (static_cast<double>(zi))*hz;
        fprintf(fp, "Point%d = %lf, %lf, %lf \n", pi, p[0], p[1], p[2]);
      }//end for xi
    }//end for yi
  }//end for zi

  fprintf(fp, "\n");

  int numElem = (nx - 1)*(ny - 1)*(nz - 1);
  fprintf(fp, "NumberOfElements = %d \n", numElem);

  for(int zi = 0, pi = 0; zi < (nz - 1); zi++) {
    for(int yi = 0; yi < (ny - 1); yi++) {
      for(int xi = 0; xi < (nx - 1); xi++, pi++) {
        int p[8];
        p[0] = __NODE__(xi, yi, zi, nx, ny);
        p[1] = __NODE__((xi + 1), yi, zi, nx, ny);
        p[2] = __NODE__((xi + 1), (yi + 1), zi, nx, ny);
        p[3] = __NODE__(xi, (yi + 1), zi, nx, ny);
        p[4] = __NODE__(xi, yi, (zi + 1), nx, ny);
        p[5] = __NODE__((xi + 1), yi, (zi + 1), nx, ny);
        p[6] = __NODE__((xi + 1), (yi + 1), (zi + 1), nx, ny);
        p[7] = __NODE__(xi, (yi + 1), (zi + 1), nx, ny);
        fprintf(fp, "Elem%d = %d, %d, %d, %d, %d, %d, %d, %d \n", pi, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
      }//end for xi
    }//end for yi
  }//end for zi

  fprintf(fp,"\n");

  fprintf(fp, "NumberOfBoundaryNodeIds = 2 \n\n");

  //z = Lz surface
  fprintf(fp, "NumberOfBoundaryNodes1 = %d \n", (nx*ny));
  fprintf(fp, "BoundaryNodeId1 = ");

  for(int yi = 0, cnt = 0; yi < ny; yi++) {
    for(int xi = 0; xi < nx; xi++, cnt++) {
      int num = __NODE__(xi, yi, (nz - 1), nx, ny);
      fprintf(fp, "%d", num);
      if(cnt < ((nx*ny) - 1)) {
        fprintf(fp, ", ");
      }
    }//end for yi
  }//end for zi

  fprintf(fp,"\n\n");

  //z = 0 surface
  fprintf(fp, "NumberOfBoundaryNodes2 = %d \n", (nx*ny));
  fprintf(fp, "BoundaryNodeId2 = ");

  for(int yi = 0, cnt = 0; yi < ny; yi++) {
    for(int xi = 0; xi < nx; xi++, cnt++) {
      int num = __NODE__(xi, yi, 0, nx, ny);
      fprintf(fp, "%d", num);
      if(cnt < ((nx*ny) - 1)) {
        fprintf(fp, ", ");
      }
    }//end for yi
  }//end for zi

  fprintf(fp,"\n\n");

  fprintf(fp, "NumberOfBoundarySideIds = 0 \n\n");

  fprintf(fp, "} \n");

  fprintf(fp, " Lx = %lf\n Ly = %lf\n Lz = %lf\n", Lx, Ly, Lz);
  fprintf(fp, " nx = %d\n ny = %d\n nz = %d\n", nx, ny, nz); 

  fprintf(fp,"\n\n");

  fclose(fp);

}



