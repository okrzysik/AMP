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
        virtual ~DataStore() = default;
        inline uint64_t getHash() const { return d_hash; }
        inline const std::string &getName() const { return d_name; }
        virtual void write( hid_t fid, const std::string &name ) const = 0;

    protected:
        uint64_t d_hash = 0;
        std::string d_name;
    };
    template<class TYPE>
    class DataStoreType : public DataStore
    {
    public:
        DataStoreType( hid_t fid, uint64_t hash, RestartManager *manager )
        {
            d_hash = hash;
            d_data = read( fid, hash2String( hash ), manager );
        }
        DataStoreType( const std::string &, std::shared_ptr<const TYPE>, RestartManager * );
        virtual ~DataStoreType() = default;
        void write( hid_t fid, const std::string &name ) const override;
        std::shared_ptr<TYPE> read( hid_t fid, const std::string &name, RestartManager * ) const;
        auto getData() { return std::const_pointer_cast<TYPE>( d_data ); }

    protected:
        std::shared_ptr<const TYPE> d_data;
    };
    using DataStorePtr = std::shared_ptr<DataStore>;

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
     * \brief  Register data with the restart manager
     * \details This function registers an object with the restart manager
     * @param[in] name      Name to use for object
     * @param[in] data      Data to register
     */
    template<class TYPE>
    void registerData( const TYPE &data, const std::string &name = "" );

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
    RestartManager::DataStorePtr create( const std::string &, std::shared_ptr<const TYPE> );

    void writeCommData( const std::string &file, Compression compress );
    void readCommData( const std::string &file );
    static std::string hash2String( uint64_t );

private:
    bool d_writer;
    hid_t d_fid;                             // fid for reading
    std::map<uint64_t, DataStorePtr> d_data; // Object data
    std::map<std::string, uint64_t> d_names; // Names to associate with the hashes
    std::map<uint64_t, AMP_MPI> d_comms;     // Registered communicators
};


} // namespace AMP::IO


#include "AMP/IO/RestartManager.hpp"


#endif
