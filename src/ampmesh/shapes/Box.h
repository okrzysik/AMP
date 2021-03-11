#ifndef included_AMP_Geometry_Box
#define included_AMP_Geometry_Box

#include "AMP/ampmesh/LogicalGeometry.h"

#include <array>
#include <vector>


namespace AMP {
namespace Geometry {


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
template<std::size_t NDIM>
class Box final : public LogicalGeometry
{
public:
    /**
     * \brief Construct a Box geometry
     * \param db        Input database
     */
    explicit Box( std::shared_ptr<AMP::Database> db );

    /**
     * \brief Construct a box
     * \param range     The range of the box [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     */
    explicit Box( const std::vector<double> &range );

public: // Default constructors
    explicit Box( Box<NDIM> &&range )      = default;
    explicit Box( const Box<NDIM> &range ) = default;
    Box<NDIM> &operator=( Box<NDIM> &&range ) = default;
    Box<NDIM> &operator=( const Box<NDIM> &range ) = default;

public: // Functions inherited from Geometry
    std::string getName() const override final;
    bool isConvex() const override final { return true; }
    Point nearest( const Point &pos ) const override final;
    double distance( const Point &pos, const Point &dir ) const override final;
    bool inside( const Point &pos ) const override final;
    int NSurface() const override final { return 2 * NDIM; }
    int surface( const Point &x ) const override final;
    Point surfaceNorm( const Point &x ) const override final;
    Point logical( const Point &x ) const override final;
    Point physical( const Point &x ) const override final;
    Point centroid() const override final;
    std::pair<Point, Point> box() const override final;
    double volume() const override final;
    void displace( const double *x ) override final;
    std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    virtual std::vector<int>
    getLogicalGridSize( const std::vector<double> &res ) const override final;
    std::vector<bool> getPeriodicDim() const override final;
    std::vector<int> getLogicalSurfaceIds() const override final;
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;
    bool operator==( const Geometry &rhs ) const override final;

protected:
    // Internal data
    std::array<double, 6> d_range;

private:
    // Private constuctor
    Box();
};


/**
 * \class Geometry
 * \brief A class used to abstract away geometry information from an application or mesh.
 * \details  This class provides routines for reading, accessing and writing geometries.
 */
template<std::size_t NDIM>
class Grid final : public LogicalGeometry
{
public:
    /**
     * \brief Construct a Grid geometry
     * \param db        Input database
     */
    explicit Grid( std::shared_ptr<AMP::Database> db );

    /**
     * \brief Construct a grid geometry
     * \param coord     Optional list of coordinates to adjust the logical - physical mapping
     */
    explicit Grid( const std::vector<std::vector<double>> &coord );

public: // Default constructors
    explicit Grid( Grid<NDIM> &&rhs )      = default;
    explicit Grid( const Grid<NDIM> &rhs ) = default;
    Grid<NDIM> &operator=( Grid<NDIM> &&rhs ) = default;
    Grid<NDIM> &operator=( const Grid<NDIM> &rhs ) = default;

public: // Functions inherited from Geometry
    std::string getName() const override final;
    bool isConvex() const override final { return true; }
    Point nearest( const Point &pos ) const override final;
    double distance( const Point &pos, const Point &dir ) const override final;
    bool inside( const Point &pos ) const override final;
    int NSurface() const override final { return 2 * NDIM; }
    int surface( const Point &x ) const override final;
    Point surfaceNorm( const Point &x ) const override final;
    Point logical( const Point &x ) const override final;
    Point physical( const Point &x ) const override final;
    Point centroid() const override final;
    std::pair<Point, Point> box() const override final;
    double volume() const override final;
    void displace( const double *x ) override final;
    std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    virtual std::vector<int>
    getLogicalGridSize( const std::vector<double> &res ) const override final;
    std::vector<bool> getPeriodicDim() const override final;
    std::vector<int> getLogicalSurfaceIds() const override final;
    std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;
    bool operator==( const Geometry &rhs ) const override final;

protected:
    // Internal data
    std::array<double, 6> d_range;
    std::vector<double> d_coord[NDIM];

private:
    // Private constuctor
    Grid();
};


} // namespace Geometry
} // namespace AMP

#endif
