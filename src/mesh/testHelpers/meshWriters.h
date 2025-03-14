#ifndef included_AMP_MeshWriters
#define included_AMP_MeshWriters

#include <string>


namespace AMP::Mesh::MeshWriters {


void writeBox( int nx, int ny, int nz, double Lx, double Ly, double Lz, const std::string &file );
void writeDispValsForPatchTest( const std::string &file );
void writeDistortedElement( const std::string &file );
void writePlateWithHole( int le,
                         int me,
                         int ne,
                         int pe,
                         double a,
                         double b,
                         double c,
                         double r,
                         const std::string &file );
void write2elementMesh(
    double a, int ny, int nz, double Lx, double Ly, double Lz, const std::string &file );
void write7elementMesh( int NumberOfBoundaryNodeIds, const std::string &file );
void writeConstrainedMesh(
    int nx, int ny, int nz, double Lx, double Ly, double Lz, const std::string &file );

} // namespace AMP::Mesh::MeshWriters


#endif
