
#define __PI__ 3.14159265

#define __SQUARE__( x ) ( ( x ) * ( x ) )

#define __SOL_FN__( f, x, y, z )                                           \
    ( ( __SQUARE__( sin( static_cast<double>( f ) * __PI__ * ( x ) ) ) ) + \
      ( __SQUARE__( sin( static_cast<double>( f ) * __PI__ * ( y ) ) ) ) + \
      ( __SQUARE__( sin( static_cast<double>( f ) * __PI__ * ( z ) ) ) ) )


AMP::LinearAlgebra::Vector::shared_ptr
getInitialGuess( AMP::LinearAlgebra::Vector::shared_ptr solution,
                 AMP::MeshManager::Adapter::shared_ptr meshAdapter,
                 AMP::Database::shared_ptr input_db )
{

    AMP::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap( 1 );

    AMP::MeshManager::Adapter::ElementIterator el     = meshAdapter->beginElement();
    AMP::MeshManager::Adapter::ElementIterator end_el = meshAdapter->endElement();

    std::vector<unsigned int> dof_indices;

    bool useConstantInitialGuess = input_db->getBool( "useConstantInitialGuess" );

    if ( useConstantInitialGuess ) {
        solution->setToScalar( 1 );
    } else {
        solution->setToScalar( 0 );

        const int numFreq = input_db->getInteger( "numFrequencies" );

        for ( ; el != end_el; ++el ) {

            dof_map->getDOFs( *el, dof_indices );

            for ( unsigned int k = 0; k < el->numNodes(); k++ ) {
                AMP::MeshManager::Adapter::Point pt = el->getPoint( k );

                double val = 0.0;
                for ( int freq = 1; freq <= numFreq; freq++ ) {
                    val += ( __SOL_FN__( freq, pt.x(), pt.y(), pt.z() ) );
                } // end for freq

                solution->setValueByGlobalID( dof_indices[k], val );
            } // end for k
        }     // end for el

        solution->makeConsistent();
    }

    double initialGuessScaleFactor = input_db->getDouble( "initialGuessScaleFactor" );

    solution->scale( initialGuessScaleFactor );

    /* Modify the vector for DIRICHLET boundary conditions */
    /* All boundary nodes are assumed to be Dirichlet. */
    el = meshAdapter->beginElement();
    for ( ; el != end_el; ++el ) {
        dof_map->getDOFs( *el, dof_indices );

        for ( unsigned int s = 0; s < el->getElem().n_sides(); s++ ) {
            if ( el->getElem().neighbor( s ) == NULL ) {
                AutoPtr<Elem> side( el->getElem().build_side( s ) );
                for ( unsigned int ns = 0; ns < side->n_nodes(); ns++ ) {
                    for ( unsigned int n = 0; n < el->getElem().n_nodes(); n++ ) {
                        if ( el->getElem().node( n ) == side->node( ns ) ) {
                            solution->setValueByGlobalID( dof_indices[n], 0.0 );
                        }
                    } // end for n
                }     // end for ns
            }
        } // end for s

    } // end for el

    solution->makeConsistent();
    return solution;
}


AMP::LinearAlgebra::Vector::shared_ptr
getInitialGuess( AMP::MeshManager::Adapter::shared_ptr meshAdapter,
                 AMP::Database::shared_ptr input_db )
{
    AMP::LinearAlgebra::Vector::shared_ptr solution = meshAdapter->getVector( 1 );
    return getInitialGuess( solution, meshAdapter, input_db );
}
