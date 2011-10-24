namespace AMP {
namespace LinearAlgebra {

  inline
  StridedVariable::StridedVariable ( const std::string &name , 
                                     size_t offset , 
                                     size_t stride )
         : SubsetVariable ( name )
         , d_Indexer ( new StridedIndexer ( offset , stride ) ) 
  {
  }

  inline
  VectorIndexer::shared_ptr  StridedVariable::getIndexer () 
  { 
    return d_Indexer; 
  }

  inline
  size_t   StridedVariable::DOFsPerObject () const 
  { 
    return 1; 
  }

}
}
