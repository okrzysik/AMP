#ifndef included_AMP_RestartManager
#define included_AMP_RestartManager

#include "AMP/IO/HDF5.h"

#include <map>
#include <string>


namespace AMP {


//! Class to manage reading/writing restart data
class RestartManager final
{
public:
    /**
     * \brief  Write the data
     * \details  Write all of the data currently registered with the manager to the disk
     * @param[in] filename  Filename to use
     */
    void write( const std::string &filename, Compression compress = Compression::None );

    /**
     * \brief  Open a restart file
     * \details  Write all of the data currently registered with the manager to the disk
     * @param[in] filename  Filename to use
     */
    hid_t open( const std::string &filename );

    /**
     * \brief  Register data with the restart manager
     * \details This function registers an object with the restart manager
     * @param[in] name      Name to use for object
     * @param[in] data      Data to register
     */
    template<class TYPE>
    void registerData( const std::string &name, std::shared_ptr<TYPE> data );

    /**
     * \brief  Get data from the restart manager
     * \details This function will get a registered/loaded object from the restart manager
     * @param[in] fid       Handle to data stored in restart manager
     * @param[in] name      Name to use for object
     * @param[in] comm      Communicator to use for the object (if parallel)
     */
    template<class TYPE>
    std::shared_ptr<TYPE>
    getData( hid_t fid, const std::string &name, const AMP::AMP_MPI &comm = AMP_COMM_WORLD );


private:
    class DataStore
    {
    public:
        virtual AMP_MPI getComm() const                          = 0;
        virtual void write( hid_t fid, const std::string &name ) = 0;
    };
    template<class TYPE>
    class DataStoreType : public DataStore
    {
    public:
        DataStoreType( std::shared_ptr<TYPE> data );
        AMP_MPI getComm() const override;
        void write( hid_t fid, const std::string &name ) override;

    private:
        std::shared_ptr<TYPE> d_data;
    };

private:
    std::map<std::string, std::unique_ptr<DataStore>> d_data;
};


} // namespace AMP


#include "AMP/IO/RestartManager.hpp"

#endif
