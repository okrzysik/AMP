namespace AMP {
namespace LinearAlgebra {

  inline
  VectorEngine::~VectorEngine () 
  {
  }

  inline
  VectorEngineParameters::shared_ptr    VectorEngine::getEngineParameters() const 
  { 
    return d_Params; 
  }

}
}

