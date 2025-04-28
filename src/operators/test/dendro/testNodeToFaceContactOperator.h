#ifndef TESTNODETOFACECONTACTOPERATOR_H
#define TESTNODETOFACECONTACTOPERATOR_H

#include "AMP/mesh/euclidean_geometry_tools.h"
#include "AMP/mesh/latex_visualization_tools.h"
#include "newton_solver.h"

class my_problem_parameters_t
{
public:
    my_problem_parameters_t()
    {
        a       = 1.05;
        b       = 2150.0;
        c       = 200.0;
        d       = 1.0;
        R_fo    = 5.17e-3;
        T_fo    = 400.0;
        q_prime = 23.9e3;
        r       = 0.0;
    }
    double a, b, c, d, R_fo, T_fo, q_prime, r;
} my_p;

static void my_f( double const &T, double &f, void *params )
{
    my_problem_parameters_t *p = static_cast<my_problem_parameters_t *>( params );

    double R_fo    = p->R_fo;
    double T_fo    = p->T_fo;
    double q_prime = p->q_prime;
    double a       = p->a;
    double b       = p->b;
    double c       = p->c;
    double d       = p->d;
    double r       = p->r;

    f = a * ( T_fo - T ) + b / d * log( ( c + T_fo * d ) / ( c + T * d ) ) +
        q_prime / ( 4.0 * M_PI ) * ( 1 - ( r * r ) / ( R_fo * R_fo ) );
}

static void my_ijmf( double const &T, double const &f, double &ijmf, void *params )
{
    my_problem_parameters_t *p = static_cast<my_problem_parameters_t *>( params );

    double a = p->a;
    double b = p->b;
    double c = p->c;
    double d = p->d;

    double df = -a - b / ( c + T * d );
    ijmf      = -f / df;
}

static void computeFuelTemperature( std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                                    AMP::LinearAlgebra::Vector::shared_ptr temperatureField,
                                    double fuelOuterRadius,
                                    double fuelOuterRadiusTemperature,
                                    double linearHeatGenerationRate,
                                    double fuelThermalConductivity )
{
    auto meshIterator       = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex );
    auto meshIterator_begin = meshIterator.begin();
    auto meshIterator_end   = meshIterator.end();
    double epsilon          = 1.0e-14;
    double temperature;
    my_p.R_fo    = fuelOuterRadius;
    my_p.T_fo    = fuelOuterRadiusTemperature;
    my_p.q_prime = linearHeatGenerationRate;
    newton_solver_t<double> solver;
    solver.set( &my_f, &my_ijmf );
    std::vector<double> vertexCoord;
    std::vector<size_t> DOFsIndices;
    for ( meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator ) {
        auto vertexCoord = meshIterator->coord();
        auto radius =
            std::sqrt( vertexCoord[0] * vertexCoord[0] + vertexCoord[1] * vertexCoord[1] );
        AMP_ASSERT( radius <= fuelOuterRadius + epsilon );
        my_p.r = radius;
        solver.apply( temperature, &my_p );
        temperatureField->getDOFManager()->getDOFs( meshIterator->globalID(), DOFsIndices );
        AMP_ASSERT( DOFsIndices.size() == 1 );
        temperatureField->setLocalValuesByGlobalID( 1, &( DOFsIndices[0] ), &temperature );
    } // end for
}

static void computeCladTemperature( std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                                    AMP::LinearAlgebra::Vector::shared_ptr temperatureField,
                                    double cladInnerRadius,
                                    double cladOuterRadius,
                                    double linearHeatGenerationRate,
                                    double cladOuterRadiusTemperature,
                                    double cladThermalConductivity )
{
    auto meshIterator       = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex );
    auto meshIterator_begin = meshIterator.begin();
    auto meshIterator_end   = meshIterator.end();
    std::vector<size_t> DOFsIndices;
    double epsilon = 1.0e-14;
    for ( meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator ) {
        auto vertexCoord = meshIterator->coord();
        auto radius =
            std::sqrt( vertexCoord[0] * vertexCoord[0] + vertexCoord[1] * vertexCoord[1] );
        AMP_ASSERT( ( radius >= cladInnerRadius - epsilon ) &&
                    ( radius <= cladOuterRadius + epsilon ) );
        auto temperature = cladOuterRadiusTemperature + linearHeatGenerationRate / ( 2.0 * M_PI ) /
                                                            cladThermalConductivity *
                                                            std::log( cladOuterRadius / radius );
        temperatureField->getDOFManager()->getDOFs( meshIterator->globalID(), DOFsIndices );
        AMP_ASSERT( DOFsIndices.size() == 1 );
        temperatureField->setLocalValuesByGlobalID( 1, &( DOFsIndices[0] ), &temperature );
    } // end for
}

static void
makeConstraintsOnFuel( AMP::Mesh::MeshIterator it,
                       double fuelOuterRadius,
                       std::max<AMP::Mesh::MeshElementID, std::map<size_t, double>> &constraints,
                       bool option,
                       double dishRadius = 0.0 )
{
    AMP::Mesh::MeshIterator it_begin = it.begin();
    AMP::Mesh::MeshIterator it_end   = it.end();
    double epsilon                   = 1.0e-10, radius;
    std::vector<double> coord;
    std::max<size_t, double>;
    tmp;
    for ( it = it_begin; it != it_end; ++it ) {
        coord  = it->coord();
        radius = std::sqrt( std::pow( coord[0], 2 ) + std::pow( coord[1], 2 ) );
        AMP_ASSERT( radius <= fuelOuterRadius + epsilon );
        tmp.clear();
        if ( std::abs( fuelOuterRadius - radius ) > epsilon ) {
            if ( std::abs( coord[0] - 0.0 ) < epsilon ) {
                tmp.insert( std::pair<size_t, double>( 0, 0.0 ) );
                if ( std::abs( coord[1] - 0.0 ) < epsilon ) {
                    tmp.insert( std::pair<size_t, double>( 1, 0.0 ) );
                    if ( option ) {
                        tmp.insert( std::pair<size_t, double>( 2, 0.0 ) );
                    }
                } // end if
            } else if ( std::abs( coord[1] - 0.0 ) < epsilon ) {
                tmp.insert( std::pair<size_t, double>( 1, 0.0 ) );
            } // end if
            if ( option && ( dishRadius > 0.0 ) ) {
                if ( std::abs( dishRadius - radius ) < epsilon ) {
                    tmp.insert( std::pair<size_t, double>( 2, 0.0 ) );
                } // end if
            }     // end if
            if ( !tmp.empty() ) {
                constraints.insert( std::pair<AMP::Mesh::MeshElementID, std::map<size_t, double>>(
                    it->globalID(), tmp ) );
            } // end if
        }     // end if
    }         // end for
}

static void
makeConstraintsOnClad( AMP::Mesh::MeshIterator it,
                       double cladInnerRadius,
                       double cladOuterRadius,
                       double cladHeight,
                       std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> &constraints )
{
    auto it_begin  = it.begin();
    auto it_end    = it.end();
    double epsilon = 1.0e-10, radius;
    std::vector<double> coord;
    std::map<size_t, double> tmp;
    for ( it = it_begin; it != it_end; ++it ) {
        coord  = it->coord();
        radius = std::sqrt( std::pow( coord[0], 2 ) + std::pow( coord[1], 2 ) );
        AMP_ASSERT( ( radius >= cladInnerRadius - epsilon ) &&
                    ( radius <= cladOuterRadius + epsilon ) );
        tmp.clear();
        if ( ( std::abs( coord[2] - 0.5 * cladHeight ) < epsilon ) ||
             ( std::abs( coord[2] + 0.5 * cladHeight ) < epsilon ) ) {
            if ( std::abs( radius - cladInnerRadius ) > epsilon ) {
                if ( std::abs( coord[0] - 0.0 ) < epsilon ) {
                    tmp.insert( std::pair<size_t, double>( 0, 0.0 ) );
                } else if ( std::abs( coord[1] - 0.0 ) < epsilon ) {
                    tmp.insert( std::pair<size_t, double>( 1, 0.0 ) );
                } // end if
            }     // end if
        } else if ( std::abs( coord[2] - 0.0 ) < epsilon ) {
            if ( std::abs( radius - cladOuterRadius ) < epsilon ) {
                tmp.insert( std::pair<size_t, double>( 2, 0.0 ) );
            } // end if
        }     // end if
        if ( !tmp.empty() ) {
            constraints.insert( std::pair<AMP::Mesh::MeshElementID, std::map<size_t, double>>(
                it->globalID(), tmp ) );
        } // end if
    }     // end for
}


static void applyCustomDirichletCondition(
    AMP::LinearAlgebra::Vector::shared_ptr rhs,
    AMP::LinearAlgebra::Vector::shared_ptr &cor,
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> const &constraints,
    std::shared_ptr<AMP::LinearAlgebra::Matrix> mat )
{
    auto dofManager = rhs->getDOFManager();
    std::vector<size_t> dofIndices;
    if ( mat.get() != NULL ) {
        char xyz[]   = "xyz";
        size_t count = 0;
        std::vector<double> coord;
        auto dir = rhs->clone();
        dir->zero();
        for ( auto it = constraints.begin(); it != constraints.end(); ++it ) {
            AMP_ASSERT( ( meshAdapter->getElement( it->first ) ).isOnSurface() );
            dofManager->getDOFs( it->first, dofIndices );
            coord = ( meshAdapter->getElement( it->first ) ).coord();
            for ( auto jt = it->second.begin(); jt != it->second.end(); ++jt ) {
                std::cout << ++count << " (" << coord[0] << ", " << coord[1] << ", " << coord[2]
                          << ") " << jt->second << "  " << xyz[jt->first] << "\n";
                dir->setLocalValueByGlobalID( dofIndices[jt->first], jt->second );
            } // end for
        }     // end for
        dir->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
        mat->mult( dir, cor );
        for ( auto it = constraints.begin(); it != constraints.end(); ++it ) {
            dofManager->getDOFs( it->first, dofIndices );
            for ( auto jt = it->second.begin(); jt != it->second.end(); ++jt ) {
                for ( size_t k = 0; k < dofIndices.size(); ++k ) {
                    if ( jt->first == k ) {
                        mat->setValueByGlobalID( dofIndices[k], dofIndices[k], 1.0 );
                    } else {
                        mat->setValueByGlobalID( dofIndices[jt->first], dofIndices[k], 0.0 );
                        mat->setValueByGlobalID( dofIndices[k], dofIndices[jt->first], 0.0 );
                    } // end if
                }     // end for k
            }         // end for jt
            auto neighbors = ( meshAdapter->getElement( it->first ) ).getNeighbors();
            std::vector<size_t> neighborsDofIndices;
            for ( size_t n = 0; n < neighbors.size(); ++n ) {
                dofManager->getDOFs( neighbors[n]->globalID(), neighborsDofIndices );
                for ( size_t k = 0; k < neighborsDofIndices.size(); ++k ) {
                    for ( auto jt = it->second.begin(); jt != it->second.end(); ++jt ) {
                        mat->setValueByGlobalID(
                            dofIndices[jt->first], neighborsDofIndices[k], 0.0 );
                        mat->setValueByGlobalID(
                            neighborsDofIndices[k], dofIndices[jt->first], 0.0 );
                    } // end for jt
                }     // end for k
            }         // end for n
        }             // end for it
        mat->makeConsistent();
    } // end if
    rhs->subtract( rhs, cor );
    for ( auto it = constraints.begin(); it != constraints.end(); ++it ) {
        dofManager->getDOFs( it->first, dofIndices );
        for ( auto jt = it->second.begin(); jt != it->second.end(); ++jt ) {
            rhs->setLocalValueByGlobalID( dofIndices[jt->first], jt->second );
        } // end for jt
    }     // end for it
}

static void shrinkMesh( std::shared_ptr<AMP::Mesh::Mesh> mesh, double const shrinkFactor )
{
    auto dofManager = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 0, 3, false );
    auto dispVar = std::make_shared<AMP::LinearAlgebra::Variable>( "disp" );
    auto dispVec = AMP::LinearAlgebra::createVector( dofManager, dispVar, false );
    dispVec->zero();
    auto meshIterator       = mesh->getIterator( AMP::Mesh::GeomType::Vertex );
    auto meshIterator_begin = meshIterator.begin();
    auto meshIterator_end   = meshIterator.end();
    std::vector<double> vertexCoord( 3 ), vertexDisp( 3 );
    std::vector<size_t> dofIndices;
    for ( meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator ) {
        vertexCoord = meshIterator->coord();
        dofManager->getDOFs( meshIterator->globalID(), dofIndices );
        AMP_ASSERT( dofIndices.size() == 3 );
        std::transform( vertexCoord.begin(),
                        vertexCoord.end(),
                        vertexDisp.begin(),
                        std::bind1st( std::multiplies<double>(), -shrinkFactor ) );
        dispVec->setLocalValuesByGlobalID( 3, &( dofIndices[0] ), &( vertexDisp[0] ) );
    } // end for
    mesh->displaceMesh( dispVec );
}

static void rotateMesh( std::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    auto dofManager = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 0, 3, false );
    auto dispVar = std::make_shared<AMP::LinearAlgebra::Variable>( "disp" );
    auto dispVec = AMP::LinearAlgebra::createVector( dofManager, dispVar, false );
    dispVec->zero();
    auto meshIterator       = mesh->getIterator( AMP::Mesh::GeomType::Vertex );
    auto meshIterator_begin = meshIterator.begin();
    auto meshIterator_end   = meshIterator.end();
    for ( meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator ) {
        auto oldGeomType::VertexCoord = meshIterator->coord();
        auto dofIndices;
        dofManager->getDOFs( meshIterator->globalID(), dofIndices );
        AMP_ASSERT( dofIndices.size() == 3 );
        auto newGeomType::VertexCoord( oldGeomType::VertexCoord );
        // 15 degrees around the z axis
        //    rotate_points(2, M_PI / 12.0, 1, &(newGeomType::VertexCoord[0]));
        rotate_points( 2, M_PI / 2.0, 1, &( newGeomType::VertexCoord[0] ) );
        std::vector<double> vertexDisp( 3, 0.0 );
        make_vector_from_two_points( &( oldGeomType::VertexCoord[0] ),
                                     &( newGeomType::VertexCoord[0] ),
                                     &( vertexDisp[0] ) );
        dispVec->setLocalValuesByGlobalID( 3, &( dofIndices[0] ), &( vertexDisp[0] ) );
    } // end for
    mesh->displaceMesh( dispVec );
}

struct dummyClass {
    bool operator()( std::pair<AMP::Mesh::MeshElementID, size_t> const &lhs,
                     std::pair<AMP::Mesh::MeshElementID, size_t> const &rhs )
    {
        return lhs.first < rhs.first;
    }
};


static void computeStressTensor( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                 AMP::LinearAlgebra::Vector::shared_ptr displacementField,
                                 AMP::LinearAlgebra::Vector::shared_ptr sigma_xx,
                                 AMP::LinearAlgebra::Vector::shared_ptr sigma_yy,
                                 AMP::LinearAlgebra::Vector::shared_ptr sigma_zz,
                                 AMP::LinearAlgebra::Vector::shared_ptr sigma_yz,
                                 AMP::LinearAlgebra::Vector::shared_ptr sigma_xz,
                                 AMP::LinearAlgebra::Vector::shared_ptr sigma_xy,
                                 AMP::LinearAlgebra::Vector::shared_ptr sigma_eff,
                                 double YoungsModulus,
                                 double PoissonsRatio,
                                 double referenceTemperature        = 273.0,
                                 double thermalExpansionCoefficient = 2.0e-6,
                                 AMP::LinearAlgebra::Vector::shared_ptr temperatureField =
                                     AMP::LinearAlgebra::Vector::shared_ptr() )
{

    sigma_eff->setToScalar( -1.0 );

    double constitutiveMatrix[36];
    compute_constitutive_matrix( YoungsModulus, PoissonsRatio, constitutiveMatrix );

    double verticesLocalCoordinates[24] = { -1.0, -1.0, -1.0, +1.0, -1.0, -1.0, +1.0, +1.0,
                                            -1.0, -1.0, +1.0, -1.0, -1.0, -1.0, +1.0, +1.0,
                                            -1.0, +1.0, +1.0, +1.0, +1.0, -1.0, +1.0, +1.0 };

    std::set<std::pair<AMP::Mesh::MeshElementID, size_t>, dummyClass> verticesGlobalIDsAndCount;

    size_t numLocalVertices = mesh->numLocalElements( AMP::Mesh::GeomType::Vertex );
    size_t countVertices    = 0;
    std::vector<size_t> verticesInHowManyGeomType::CellElements( numLocalVertices, 1 );
    std::vector<size_t> verticesDOFIndex( numLocalVertices );
    auto meshIterator       = mesh->getIterator( AMP::Mesh::GeomType::Cell );
    auto meshIterator_begin = meshIterator.begin();
    auto meshIterator_end   = meshIterator.end();
    for ( meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator ) {
        auto volumeVertices = meshIterator->getElements( AMP::Mesh::GeomType::Vertex );
        AMP_ASSERT( volumeVertices.size() == 8 );
        std::vector<AMP::Mesh::MeshElementID> volumeVerticesGlobalIDs( 8 );
        double volumeVerticesCoordinates[24];
        for ( size_t v = 0; v < 8; ++v ) {
            std::vector<double> vertexCoord = volumeVertices[v].coord();
            std::copy(
                vertexCoord.begin(), vertexCoord.end(), &( volumeVerticesCoordinates[3 * v] ) );
            volumeVerticesGlobalIDs[v] = volumeVertices[v].globalID();
        } // end for v
        hex8_element_t volumeElement( volumeVerticesCoordinates );

        std::vector<size_t> displacementIndices;
        displacementField->getDOFManager()->getDOFs( volumeVerticesGlobalIDs, displacementIndices );
        AMP_ASSERT( displacementIndices.size() == 24 );
        double displacementValues[24];
        displacementField->getValuesByGlobalID(
            24, &( displacementIndices[0] ), &( displacementValues[0] ) );

        double temperatureValues[8];
        if ( temperatureField.get() != NULL ) {
            std::vector<size_t> temperatureIndices;
            temperatureField->getDOFManager()->getDOFs( volumeVerticesGlobalIDs,
                                                        temperatureIndices );
            AMP_ASSERT( temperatureIndices.size() == 8 );
            temperatureField->getValuesByGlobalID(
                8, &( temperatureIndices[0] ), &( temperatureValues[0] ) );
        } else {
            std::fill( temperatureValues, temperatureValues + 8, referenceTemperature );
        } // end if

        double strainTensor[6];
        double stressTensor[6];
        for ( size_t v = 0; v < 8; ++v ) {
            volumeElement.compute_strain_tensor(
                &( verticesLocalCoordinates[3 * v] ), displacementValues, strainTensor );
            if ( temperatureField.get() != NULL ) {
                for ( size_t i = 0; i < 3; ++i ) {
                    strainTensor[i] -= ( thermalExpansionCoefficient *
                                         ( temperatureValues[v] - referenceTemperature ) );
                } // end for i
            }     // end if
            compute_stress_tensor( constitutiveMatrix, strainTensor, stressTensor );
            double vonMisesStress = compute_von_mises_stress( stressTensor );

            std::vector<size_t> indices;
            sigma_xx->getDOFManager()->getDOFs( volumeVerticesGlobalIDs[v], indices );
            AMP_ASSERT( indices.size() == 1 );
            // good job buddy
            std::pair<std::set<std::pair<AMP::Mesh::MeshElementID, size_t>, dummyClass>::iterator,
                      bool>
                dummy =
                    verticesGlobalIDsAndCount.insert( std::pair<AMP::Mesh::MeshElementID, size_t>(
                        volumeVerticesGlobalIDs[v], countVertices ) );
            if ( dummy.second ) {
                verticesDOFIndex[countVertices] = indices[0];
                ++countVertices;
                sigma_xx->setLocalValueByGlobalID( indices[0], stressTensor[0] );
                sigma_yy->setLocalValueByGlobalID( indices[0], stressTensor[1] );
                sigma_zz->setLocalValueByGlobalID( indices[0], stressTensor[2] );
                sigma_yz->setLocalValueByGlobalID( indices[0], stressTensor[3] );
                sigma_xz->setLocalValueByGlobalID( indices[0], stressTensor[4] );
                sigma_xy->setLocalValueByGlobalID( indices[0], stressTensor[5] );
                sigma_eff->setLocalValueByGlobalID( indices[0], vonMisesStress );
            } else {
                // sigh...
                ++verticesInHowManyGeomType::CellElements[dummy.first->second];
                sigma_xx->addLocalValueByGlobalID( indices[0], stressTensor[0] );
                sigma_yy->addLocalValueByGlobalID( indices[0], stressTensor[1] );
                sigma_zz->addLocalValueByGlobalID( indices[0], stressTensor[2] );
                sigma_yz->addLocalValueByGlobalID( indices[0], stressTensor[3] );
                sigma_xz->addLocalValueByGlobalID( indices[0], stressTensor[4] );
                sigma_xy->addLocalValueByGlobalID( indices[0], stressTensor[5] );
                sigma_eff->addLocalValueByGlobalID( indices[0], vonMisesStress );

                /*        double deltaStressTensor[6];
                        deltaStressTensor[0] = sigma_xx->getLocalValueByGlobalID(indices[0]) -
                   stressTensor[0];
                        deltaStressTensor[1] = sigma_yy->getLocalValueByGlobalID(indices[0]) -
                   stressTensor[1];
                        deltaStressTensor[2] = sigma_zz->getLocalValueByGlobalID(indices[0]) -
                   stressTensor[2];
                        deltaStressTensor[3] = sigma_yz->getLocalValueByGlobalID(indices[0]) -
                   stressTensor[3];
                        deltaStressTensor[4] = sigma_xz->getLocalValueByGlobalID(indices[0]) -
                   stressTensor[4];
                        deltaStressTensor[5] = sigma_xy->getLocalValueByGlobalID(indices[0]) -
                   stressTensor[5];
                        double deltaVonMisesStress = sigma_eff->getLocalValueByGlobalID(indices[0])
                   - vonMisesStress;
                        std::cout<<"xx  "<<deltaStressTensor[0]<<"
                   "<<100.0*deltaStressTensor[0]/std::abs(stressTensor[0])<<" %  |  ";
                        std::cout<<"yy  "<<deltaStressTensor[1]<<"
                   "<<100.0*deltaStressTensor[1]/std::abs(stressTensor[1])<<" %  |  ";
                        std::cout<<"zz  "<<deltaStressTensor[2]<<"
                   "<<100.0*deltaStressTensor[2]/std::abs(stressTensor[2])<<" %  |  ";
                        std::cout<<"yz  "<<deltaStressTensor[3]<<"
                   "<<100.0*deltaStressTensor[3]/std::abs(stressTensor[3])<<" %  |  ";
                        std::cout<<"xz  "<<deltaStressTensor[4]<<"
                   "<<100.0*deltaStressTensor[4]/std::abs(stressTensor[4])<<" %  |  ";
                        std::cout<<"xy  "<<deltaStressTensor[5]<<"
                   "<<100.0*deltaStressTensor[5]/std::abs(stressTensor[5])<<" %  |  ";
                        std::cout<<"eff  "<<deltaVonMisesStress<<"
                   "<<100.0*deltaVonMisesStress/std::abs(vonMisesStress)<<" % \n";
                */
            } // end if

        } // end for v
    }     // end for
    AMP_ASSERT( verticesGlobalIDsAndCount.size() == numLocalVertices );
    AMP_ASSERT( countVertices == numLocalVertices );
    AMP_ASSERT( find( verticesInHowManyGeomType::CellElements.begin(),
                      verticesInHowManyGeomType::CellElements.end(),
                      0 ) == verticesInHowManyGeomType::CellElements.end() );

    for ( size_t v = 0; v < numLocalVertices; ++v ) {
        if ( verticesInHowManyGeomType::CellElements[v] > 1 ) {
            //      AMP_ASSERT(verticesInHowManyGeomType::CellElements[v] < 9);
            sigma_xx->setLocalValueByGlobalID(
                verticesDOFIndex[v],
                sigma_xx->getLocalValueByGlobalID( verticesDOFIndex[v] ) /
                    static_cast<double>( verticesInHowManyGeomType::CellElements[v] ) );
            sigma_yy->setLocalValueByGlobalID(
                verticesDOFIndex[v],
                sigma_yy->getLocalValueByGlobalID( verticesDOFIndex[v] ) /
                    static_cast<double>( verticesInHowManyGeomType::CellElements[v] ) );
            sigma_zz->setLocalValueByGlobalID(
                verticesDOFIndex[v],
                sigma_zz->getLocalValueByGlobalID( verticesDOFIndex[v] ) /
                    static_cast<double>( verticesInHowManyGeomType::CellElements[v] ) );
            sigma_yz->setLocalValueByGlobalID(
                verticesDOFIndex[v],
                sigma_yz->getLocalValueByGlobalID( verticesDOFIndex[v] ) /
                    static_cast<double>( verticesInHowManyGeomType::CellElements[v] ) );
            sigma_xz->setLocalValueByGlobalID(
                verticesDOFIndex[v],
                sigma_xz->getLocalValueByGlobalID( verticesDOFIndex[v] ) /
                    static_cast<double>( verticesInHowManyGeomType::CellElements[v] ) );
            sigma_xy->setLocalValueByGlobalID(
                verticesDOFIndex[v],
                sigma_xy->getLocalValueByGlobalID( verticesDOFIndex[v] ) /
                    static_cast<double>( verticesInHowManyGeomType::CellElements[v] ) );
            sigma_eff->setLocalValueByGlobalID(
                verticesDOFIndex[v],
                sigma_eff->getLocalValueByGlobalID( verticesDOFIndex[v] ) /
                    static_cast<double>( verticesInHowManyGeomType::CellElements[v] ) );
        } // end if
    }     // end for v
}


static void
computeStressTensor( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                     AMP::LinearAlgebra::Vector::shared_ptr displacementField,
                     AMP::LinearAlgebra::Vector::shared_ptr sigmaXX,
                     AMP::LinearAlgebra::Vector::shared_ptr sigmaYY,
                     AMP::LinearAlgebra::Vector::shared_ptr sigmaZZ,
                     AMP::LinearAlgebra::Vector::shared_ptr sigmaYZ,
                     AMP::LinearAlgebra::Vector::shared_ptr sigmaXZ,
                     AMP::LinearAlgebra::Vector::shared_ptr sigmaXY,
                     AMP::LinearAlgebra::Vector::shared_ptr sigmaEff,
                     std::shared_ptr<AMP::Operator::MechanicsMaterialModel> mechanicsMaterialModel,
                     double referenceTemperature        = 273.0,
                     double thermalExpansionCoefficient = 2.0e-6,
                     AMP::LinearAlgebra::Vector::shared_ptr temperatureField =
                         AMP::LinearAlgebra::Vector::shared_ptr() )
{

    AMP::LinearAlgebra::VS_Mesh vectorSelector( mesh );
    auto subsetDisplacementField = displacementField->select(
        vectorSelector, ( displacementField->getVariable() )->getName() );
    auto subsetSigmaXX          = sigmaXX->select( vectorSelector );
    auto subsetSigmaYY          = sigmaYY->select( vectorSelector );
    auto subsetSigmaZZ          = sigmaZZ->select( vectorSelector );
    auto subsetSigmaYZ          = sigmaYZ->select( vectorSelector );
    auto subsetSigmaXZ          = sigmaXZ->select( vectorSelector );
    auto subsetSigmaXY          = sigmaXY->select( vectorSelector );
    auto subsetSigmaEff         = sigmaEff->select( vectorSelector );
    auto subsetTemperatureField = temperatureField->select( vectorSelector );

    double youngsModulus =
        std::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>( mechanicsMaterialModel )
            ->getYoungsModulus();
    double poissonsRatio =
        std::dynamic_pointer_cast<AMP::Operator::IsotropicElasticModel>( mechanicsMaterialModel )
            ->getPoissonsRatio();

    computeStressTensor( mesh,
                         subsetDisplacementField,
                         subsetSigmaXX,
                         subsetSigmaYY,
                         subsetSigmaZZ,
                         subsetSigmaYZ,
                         subsetSigmaXZ,
                         subsetSigmaXY,
                         subsetSigmaEff,
                         youngsModulus,
                         poissonsRatio,
                         referenceTemperature,
                         thermalExpansionCoefficient,
                         subsetTemperatureField );
}


static void drawVerticesOnBoundaryID( std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                                      int boundaryID,
                                      std::ostream &os,
                                      double const *point_of_view,
                                      const std::string &option = "" )
{
    auto boundaryIterator =
        meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, boundaryID );
    auto boundaryIterator_begin = boundaryIterator.begin();
    auto boundaryIterator_end   = boundaryIterator.end();
    std::vector<double> vertexCoordinates;
    os << std::setprecision( 6 ) << std::fixed;
    for ( boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end;
          ++boundaryIterator ) {
        vertexCoordinates = boundaryIterator->coord();
        AMP_ASSERT( vertexCoordinates.size() == 3 );
        draw_point( &( vertexCoordinates[0] ), option, os );
    } // end for
}

static void drawGeomType::FacesOnBoundaryID( std::shared_ptr<AMP::Mesh::Mesh> meshAdapter,
                                             int boundaryID,
                                             std::ostream &os,
                                             double const *point_of_view,
                                             const std::string &option = "" )
{
    auto boundaryIterator =
        meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Face, boundaryID );
    auto boundaryIterator_begin = boundaryIterator.begin();
    auto boundaryIterator_end   = boundaryIterator.end();
    std::vector<AMP::Mesh::MeshElement> faceVertices;
    std::vector<double> faceGeomType::VertexCoordinates;
    double faceData[12]          = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double const *faceDataPtr[4] = { faceData, faceData + 3, faceData + 6, faceData + 9 };
    os << std::setprecision( 6 ) << std::fixed;
    for ( boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end;
          ++boundaryIterator ) {
        faceVertices = boundaryIterator->getElements( AMP::Mesh::GeomType::Vertex );
        AMP_ASSERT( faceVertices.size() == 4 );
        for ( size_t i = 0; i < 4; ++i ) {
            faceGeomType::VertexCoordinates = faceVertices[i].coord();
            AMP_ASSERT( faceGeomType::VertexCoordinates.size() == 3 );
            std::copy( faceGeomType::VertexCoordinates.begin(),
                       faceGeomType::VertexCoordinates.end(),
                       faceData + 3 * i );
        } // end for i
        triangle_t t( faceDataPtr[0], faceDataPtr[1], faceDataPtr[2] );

        //    if (compute_scalar_product(point_of_view, t.get_normal()) > 0.0) {
        if ( true ) {
            os << "\\draw[" << option << "]\n";
            write_face( faceDataPtr, os );
        } // end if
    }     // end for
}

static void myPCG( AMP::LinearAlgebra::Vector::shared_ptr rhs,
                   AMP::LinearAlgebra::Vector::shared_ptr sol,
                   AMP::Operator::Operator::shared_ptr op,
                   std::shared_ptr<AMP::Solver::SolverStrategy> pre,
                   size_t maxIters,
                   double relTol,
                   double absTol,
                   bool verbose     = false,
                   std::ostream &os = std::cout )
{
    auto res    = sol->clone();
    auto dir    = sol->clone();
    auto ext    = sol->clone();
    auto oldSol = sol->clone();
    auto oldRes = sol->clone();
    auto oldDir = sol->clone();
    auto matVec = sol->clone();
    auto nullVec;

    op->apply( nullVec, sol, matVec, 1.0, 0.0 );
    oldRes->subtract( rhs, matVec );
    pre->apply( oldRes, ext );
    oldDir->copyVector( ext );
    oldSol->copyVector( sol );
    double initialResNorm = oldRes->L2Norm();
    AMP::AMP_MPI comm     = sol->getComm();
    int rank              = comm.getRank();
    verbose               = verbose && !rank;
    double tol            = absTol + relTol * initialResNorm;
    if ( verbose ) {
        os << std::setprecision( 15 ) << "  iter=0  itialResNorm=" << initialResNorm << "\n";
    }
    for ( size_t iter = 0; iter < maxIters; ++iter ) {
        if ( verbose ) {
            os << "  iter=" << iter + 1 << "  ";
        }
        op->apply( nullVec, oldDir, matVec, 1.0, 0.0 );
        double extDOToldRes    = ext->dot( oldRes );
        double oldDirDOTmatVec = oldDir->dot( matVec );
        double alpha           = extDOToldRes / oldDirDOTmatVec;
        if ( verbose ) {
            os << "alpha=" << alpha << "  ";
        }
        if ( verbose ) {
            os << "oldDirDOTmatVec=" << oldDirDOTmatVec << "  ";
        }
        sol->axpy( alpha, oldDir, oldSol );
        res->axpy( -alpha, matVec, oldRes );
        double resNorm = res->L2Norm();
        if ( verbose ) {
            os << "resNorm=" << resNorm << "  ";
        }
        if ( resNorm < tol ) {
            os << "\n";
            break;
        }
        pre->apply( res, ext );
        double extDOTres = ext->dot( res );
        double beta      = extDOTres / extDOToldRes;
        if ( verbose ) {
            os << "beta=" << beta << "\n";
        }
        dir->axpy( beta, oldDir, ext );
        oldSol->copyVector( sol );
        oldRes->copyVector( res );
        oldDir->copyVector( dir );
    } // end for
}

#endif // TESTNODETOFACECONTACTOPERATOR_H
