#ifndef included_AMP_Geometry_Box
#define included_AMP_Geometry_Box

#include "AMP/ampmesh/LogicalGeometry.h"

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
     * \param coord     Optional list of coordinates to adjust the logical - physical mapping
     */
    explicit Box( const std::vector<double> &range );

public: // Default constructors
    explicit Box( Box<NDIM> &&range )      = default;
    explicit Box( const Box<NDIM> &range ) = default;
    Box<NDIM> &operator=( Box<NDIM> &&range ) = default;
    Box<NDIM> &operator=( const Box<NDIM> &range ) = default;

public: // Functions inherited from Geometry
    virtual std::string getName() const override final;
    virtual bool isConvex() const override final { return true; }
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int NSurface() const override final { return 2 * NDIM; }
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual Point centroid() const override final;
    virtual std::pair<Point, Point> box() const override final;
    virtual void displace( const double *x ) override final;
    virtual std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    virtual std::vector<bool> getPeriodicDim() const override final;
    virtual std::vector<int> getLogicalSurfaceIds() const override final;
    virtual std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;

protected:
    // Internal data
    double d_range[6];

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
     * \param range     The range of the box [xmin, xmax, ymin, ymax, zmin, zmax, ...]
     * \param coord     Optional list of coordinates to adjust the logical - physical mapping
     */
    explicit Grid( const std::vector<std::vector<double>> &coord );

public: // Default constructors
    explicit Grid( Grid<NDIM> &&rhs )      = default;
    explicit Grid( const Grid<NDIM> &rhs ) = default;
    Grid<NDIM> &operator=( Grid<NDIM> &&rhs ) = default;
    Grid<NDIM> &operator=( const Grid<NDIM> &rhs ) = default;

public: // Functions inherited from Geometry
    virtual std::string getName() const override final;
    virtual bool isConvex() const override final { return true; }
    virtual double distance( const Point &pos, const Point &dir ) const override final;
    virtual bool inside( const Point &pos ) const override final;
    virtual int NSurface() const override final { return 2 * NDIM; }
    virtual int surface( const Point &x ) const override final;
    virtual Point surfaceNorm( const Point &x ) const override final;
    virtual Point logical( const Point &x ) const override final;
    virtual Point physical( const Point &x ) const override final;
    virtual Point centroid() const override final;
    virtual std::pair<Point, Point> box() const override final;
    virtual void displace( const double *x ) override final;
    virtual std::vector<int> getLogicalGridSize( const std::vector<int> &x ) const override final;
    virtual std::vector<bool> getPeriodicDim() const override final;
    virtual std::vector<int> getLogicalSurfaceIds() const override final;
    virtual std::unique_ptr<AMP::Geometry::Geometry> clone() const override final;

protected:
    // Internal data
    double d_range[6];
    std::vector<double> d_coord[NDIM];

private:
    // Private constuctor
    Grid();
};


} // namespace Geometry
} // namespace AMP

#endif
