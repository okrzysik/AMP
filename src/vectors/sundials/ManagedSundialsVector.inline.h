
namespace AMP {
namespace LinearAlgebra {


inline ManagedVector *ManagedSundialsVector::getNewRawPtr() const
{
    return new ManagedSundialsVector(
        std::dynamic_pointer_cast<VectorParameters>( d_pParameters ) );
}


inline std::string ManagedSundialsVector::type() const
{
    std::string retVal = "Managed SUNDIALS Vector";
    retVal += ManagedVector::type();
    return retVal;
}


inline void ManagedSundialsVector::assemble() {}
} // namespace LinearAlgebra
} // namespace AMP
