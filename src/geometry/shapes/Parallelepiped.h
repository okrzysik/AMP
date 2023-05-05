#ifndef included_AMP_Geometry_Parallelepiped
#define included_AMP_Geometry_Parallelepiped

#include "AMP/geometry/LogicalGeometry.h"

#include <array>
#include <vector>


namespace AMP::Geometry {


/**
 * \class Parallelepiped
 * \brief A class used to describe a Parallelepiped
 * \details  This class provides routines for reading, accessing and writing Parallelepiped.
 *     See https://en.wikipedia.org/wiki/Parallelepiped
 */
class Parallelepiped : public LogicalGeometry
{
public:
    /**
     * \brief Construct a Parallelepiped geometry
     * \param db        Input database
     */
    explicit Parallelepiped( std::shared_ptr<const AMP::Database> db );

    //! Construct from restart
    Parallelepiped( int64_t );

public: // Functions inherited from Geometry
    std::string getName() const override final { return "Parallelepiped"; }
    bool isConvex() const override final { return true; }
    Point nearest( const Point &pos ) const override final;
    double distance( const Point &pos, const Point &dir ) const override final;
    bool inside( const Point &pos ) const override final;
    int NSurface() const override final { return 6; }
    int surface( const Point & ) const override final;
    Point surfaceNorm( const Point & ) const override final;
    Point logical( const Point &x ) const override final;
    Point physical( const Point &x ) const override final;
    Point centroid() const override final;
    std::pair<Point, Point> box() const override final;
    double volume() const override final;
    void displace( const double *x ) override final;
    std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    virtual std::vector<int>
    getLogicalGridSize( const std::vector<double> &res ) const override final;
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;
    bool operator==( const Geometry &rhs ) const override final;
    void writeRestart( int64_t ) const override;

protected:
    // Internal data
    std::array<double, 3> d_a;      // First edge vector
    std::array<double, 3> d_b;      // Second edge vector
    std::array<double, 3> d_c;      // Third edge vector
    std::array<double, 3> d_offset; // Offset
    std::array<double, 9> d_M_inv;  // Inverse matrix used to compute logical coordinates
    double d_V;                     // Cached volume
    AMP::Mesh::Point d_n_ab;        // Normal to the plane defined by a-b
    AMP::Mesh::Point d_n_ac;        // Normal to the plane defined by a-c
    AMP::Mesh::Point d_n_bc;        // Normal to the plane defined by b-c

private:
    // Private constructor
    Parallelepiped();
};


} // namespace AMP::Geometry

#endif
