
namespace AMP {
namespace Mesh {

  template <typename MESH , 
            template <typename M> class CONTACT_ALGO , 
            template <typename M1 , typename I1> class NORMAL_FCN ,
            template <typename M2 , typename I2 , typename I4> class COARSE_SEARCH , 
            template <typename M3 , typename MAP , typename I3 , typename I5> class FINE_SEARCH>
  void ContactManagerTmpl<MESH,CONTACT_ALGO,NORMAL_FCN,COARSE_SEARCH,FINE_SEARCH>::
          computeContactSets (  typename MESH::shared_ptr  mesh1 ,
                                typename MESH::shared_ptr  mesh2 ,
                                Side1Iterator              side1Begin ,
                                Side1Iterator              side1End ,
                                Side1FaceIterator          side1FaceBegin ,
                                Side1FaceIterator          side1FaceEnd ,
                                Side2Iterator              side2Begin ,
                                Side2Iterator              side2End ,
                                Side1Set                   &side1ActiveSet ,
                                Side2Set                   &side2ActiveSet ,
                                Side1ToSide2MultiMap       &side1ToSide2Map ,
                                Side2ToSide1MultiMap       &side2ToSide1Map )
  {


    typename CoarseSearch::Parameters::shared_ptr coarseParams ( 
                         new typename CoarseSearch::Parameters ( side1Begin ,
                                                                 side1End ,
                                                                 side2Begin ,
                                                                 side2End ) );
    coarseParams->d_mesh1 = mesh1;
    coarseParams->d_mesh2 = mesh2;
    CoarseSearch cs ( coarseParams );

    Side1ToSide2MultiMap  coarseAnswer1;
    Side2ToSide1MultiMap  coarseAnswer2;
    cs.findNeighbors ( coarseAnswer1 , coarseAnswer2 );

    typename FineSearch::Parameters::shared_ptr  fineParams ( new typename FineSearch::Parameters );

    NormalMap &normals = fineParams->d_Normals;
    NormalFunction::compute ( side1FaceBegin , side1FaceEnd , normals , mesh1 , NormalFunction::Node );
    fineParams->d_mesh1 = mesh1;
    fineParams->d_mesh2 = mesh2;
    FineSearch fs ( fineParams );
    fs.refineNeighbors ( coarseAnswer1 , coarseAnswer2 , side1ToSide2Map , side2ToSide1Map );

    typename Side1ToSide2MultiMap::iterator  curPair = side1ToSide2Map.begin();
    while ( curPair != side1ToSide2Map.end() )
    {
      Side1Type  obj1 = *(curPair->first);
      Side2Type  obj2 = *(curPair->second);
      side1ActiveSet.insert ( obj1 );
      side2ActiveSet.insert ( obj2 );
      curPair++;
    } 
  }

}
}


