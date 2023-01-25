#ifndef included_AMP_RestartManager
#define included_AMP_RestartManager

#include "AMP/IO/HDF5.h"
#include "AMP/utils/TypeTraits.h"

#include <map>
#include <string>


namespace AMP::IO {


//! Class to manage reading/writing restart data
class RestartManager final
{
public:
    class DataStore
    {
    public:
        virtual uint64_t getHash() const                               = 0;
        virtual const std::string &getName() const                     = 0;
        virtual void write( hid_t fid, const std::string &name ) const = 0;
    };
    template<class TYPE>
    class DataStoreType : public DataStore
    {
    public:
        DataStoreType( const std::string &, std::shared_ptr<const TYPE>, RestartManager * );
        uint64_t getHash() const override;
        const std::string &getName() const override { return d_name; }
        void write( hid_t fid, const std::string &name ) const override;

    private:
        std::string d_name;
        std::shared_ptr<const TYPE> d_data;
    };

public:
    //! Create a writer for restart data
    RestartManager();

    //! Create a reader for restart data
    RestartManager( const std::string &filename );

    //! Copy constructor
    RestartManager( const RestartManager & ) = delete;

    //! Move operator
    RestartManager( RestartManager && ) = default;

    //! Move assignment
    RestartManager &operator=( RestartManager && ) = default;

    //! Destructor
    ~RestartManager();

    /**
     * \brief  Write the data
     * \details  Write all of the data currently registered with the manager to the disk
     * @param[in] filename  Filename to use
     */
    void write( const std::string &filename, Compression compress = Compression::None );

    /**
     * \brief  Close the restart file
     * \details  Close the restart file
     * @param[in] filename  Filename to use
     */
    void close( const std::string &filename );

    /**
     * \brief  Register data with the restart manager
     * \details This function registers an object with the restart manager
     * @param[in] name      Name to use for object
     * @param[in] data      Data to register
     */
    template<class TYPE>
    void registerData( const std::string &name, const TYPE &data );

    /**
     * \brief  Register data with the restart manager
     * \details This function registers an object with the restart manager
     * @param[in] data      Data to register
     */
    void registerData( std::shared_ptr<DataStore> data );

    /**
     * \brief  Register a communicator with the restart manager
     * \details This function registers a communicator (based on ranks) with the restart manager
     * @param[in] comm      Communicator to register
     */
    void registerComm( const AMP::AMP_MPI &comm );

    /**
     * \brief  Get data from the restart manager
     * \details This function will get a registered/loaded object from the restart manager
     * @param[in] name      Name to use for object
     */
    template<class TYPE>
    std::shared_ptr<TYPE> getData( const std::string &name );

    /**
     * \brief  Get data from the restart manager
     * \details This function will get a registered/loaded object from the restart manager
     * @param[in] hash      Object ID
     */
    template<class TYPE>
    std::shared_ptr<TYPE> getData( uint64_t hash );


    /**
     * \brief  Get the communicator from the restart manager
     * \details This function will get a registered/loaded object from the restart manager
     * @param[in] hash      Object ID
     */
    AMP_MPI getComm( uint64_t hash );


private:
    template<class TYPE>
    std::shared_ptr<RestartManager::DataStore> create( const std::string &,
                                                       std::shared_ptr<const TYPE> );

    void writeCommData( const std::string &file, Compression compress );
    void readCommData( const std::string &file );
    std::string hash2String( uint64_t );

private:
    bool d_writer;
    hid_t d_fid;                                           // fid for reading
    std::map<uint64_t, std::shared_ptr<DataStore>> d_data; // Object data
    std::map<std::string, uint64_t> d_names;               // Names to associate with the hashes
    std::map<uint64_t, AMP_MPI> d_comms;                   // Registered communicators
};


} // namespace AMP::IO


#include "AMP/IO/RestartManager.hpp"


#endif
