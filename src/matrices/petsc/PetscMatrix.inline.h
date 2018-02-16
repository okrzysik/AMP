
namespace AMP {
namespace LinearAlgebra {

inline PetscMatrix::PetscMatrix() {}

inline PetscMatrix::PetscMatrix( MatrixParameters::shared_ptr params )
    : Matrix( params ), d_MatCreatedInternally( false ), d_Mat( nullptr )
{
}

inline PetscMatrix::~PetscMatrix() {}

inline Mat &PetscMatrix::getMat() { return d_Mat; }

inline Mat PetscMatrix::getMat() const { return d_Mat; }


} // namespace LinearAlgebra
} // namespace AMP
