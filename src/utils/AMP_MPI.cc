// This file impliments a wrapper class for MPI functions

// Include mpi.h if used
#if defined(USE_MPI) || defined(USE_EXT_MPI)
    #include "mpi.h"
    #ifndef USE_MPI
        #define USE_MPI
    #endif
#endif

// Include AMP headers
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/ProfilerApp.h"

// Include all other headers
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Utilities.h"

// Include OS specific headers
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    // We are using windows
    #define USE_WINDOWS
    #include <windows.h>
    #include <process.h>
    #define sched_yield() Sleep(0)
#else
    // We are using linux
    #define USE_LINUX
    #include <sys/time.h>
    #include <sched.h>
    #include <pthread.h>
#endif


// Global variable to track create new unique comms (dup and split)
#ifndef USE_MPI
    MPI_Comm uniqueGlobalComm=11;
#endif



namespace AMP{


// Some special structs to work with MPI
#ifdef USE_MPI
    struct IntIntStruct { int j; int i; };
    struct LongIntStruct { long int j; int i; };
    struct FloatIntStruct { float f; int i; };
    struct DoubleIntStruct { double d; int i; };
#endif


// Initialized the static member variables
volatile unsigned int AMP_MPI::N_MPI_Comm_created=0;
volatile unsigned int AMP_MPI::N_MPI_Comm_destroyed=0;
int AMP_MPI::profile_level=127;


// Define a type for use with size_t
#ifdef USE_MPI
    static MPI_Datatype MPI_SIZE_T = 0x0;
    static MPI_Datatype getSizeTDataType( )
    {
        int size_int, size_long, size_longlong;
        MPI_Type_size( MPI_UNSIGNED, &size_int );
        MPI_Type_size( MPI_UNSIGNED_LONG, &size_long );
        MPI_Type_size( MPI_LONG_LONG, &size_longlong );
        if ( sizeof(size_t) == size_int ) {
            return MPI_UNSIGNED;
        } else if ( sizeof(size_t) == size_long ) {
            return MPI_UNSIGNED_LONG;
        } else if ( sizeof(size_t) == size_longlong ) {
            AMP_WARNING("Using signed long long datatype for size_t in MPI");
            return MPI_LONG_LONG;   // Note: this is not unsigned
        } else {
            AMP_ERROR("No suitable datatype found");
        }
        return 0;
    }
#endif


// Static data for asyncronous communication without MPI
// Note: these routines may not be thread-safe yet
#ifndef USE_MPI
    static const int mpi_max_tag = 0x003FFFFF;
    struct Isendrecv_struct {
        const char *data;   // Pointer to data
        int status;         // Status: 1-sending, 2-recieving
    };
    std::map<MPI_Request,Isendrecv_struct>  global_isendrecv_list;
    static MPI_Request getRequest( MPI_Comm comm, int tag )
    {
        AMP_ASSERT(tag>=0&&tag<=mpi_max_tag);
        MPI_Request request = 0;
        if ( sizeof(MPI_Request)==4 && sizeof(MPI_Comm)==4 ) {
            AMP_ASSERT(sizeof(unsigned int)==4);
            unsigned int hash = comm*0x9E3779B9;                        // 2^32*0.5*(sqrt(5)-1)
            unsigned int key  = (hash&0xFFFFFFFF)>>22;                  // Get a key 0-1024
            request = (MPI_Request) tag + (key<<22);
        } else if ( sizeof(MPI_Request)==8 && sizeof(MPI_Comm)==4 ) {
            AMP_ASSERT(sizeof(unsigned long int)==8);
            request = (MPI_Request) ((unsigned long)tag + (((unsigned long)comm)<<32));
        } else if ( sizeof(MPI_Request)==4 && sizeof(MPI_Comm)==8 ) {
            AMP_ASSERT(sizeof(unsigned long int)==8);
            unsigned long int hash = comm*0x9E3779B97F4A7C15;           // 2^64*0.5*(sqrt(5)-1)
            unsigned long int key  = (hash&0xFFFFFFFF)>>22;             // Get a key 0-1024
            request = (MPI_Request) ((unsigned int)tag + (key<<32));
        } else if ( sizeof(MPI_Request)==8 && sizeof(MPI_Comm)==8 ) {
            AMP_ASSERT(sizeof(unsigned long int)==8);
            unsigned long int hash = comm*0x9E3779B97F4A7C15;           // 2^64*0.5*(sqrt(5)-1)
            unsigned long int key  = (hash&0xFFFFFFFF);                 // Get a key 0-2^32-1
            request = (MPI_Request) ((unsigned int)tag + (key<<32));
        } else {
            char text[50];
            sprintf(text,"Not Programmed (%i,%i)",(int)sizeof(MPI_Request),(int)sizeof(MPI_Comm));
            AMP_ERROR(std::string(text));
        }
        return request;
    }
#endif


/******************************************************************
* Some helper function to safely increment/decrement an integer   *
* in a multi-threaded environment.                                *
* These functions return the post increment/decrement count.      *
******************************************************************/
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    static HANDLE global_mpi_lock_queue;    // Mutex lock for changing the queue
    static inline int increment_count( int* volatile count )
    {
        WaitForSingleObject(global_mpi_lock_queue,INFINITE);
        int tmp = ++(*count);
        ReleaseMutex(global_mpi_lock_queue);
        return tmp;
    }
    static inline int decrement_count( int* volatile count )
    {
        WaitForSingleObject(global_mpi_lock_queue,INFINITE);
        int tmp = --(*count);
        ReleaseMutex(global_mpi_lock_queue);
        return tmp;
    }
#else
    static pthread_mutex_t global_mpi_lock_queue;    // Mutex lock for changing the queue
    static inline int increment_count( int* volatile count )
    {
        pthread_mutex_lock(&global_mpi_lock_queue);
        int tmp = ++(*count);
        pthread_mutex_unlock(&global_mpi_lock_queue);
        return tmp;
    }
    static inline int decrement_count( int* volatile count )
    {
        pthread_mutex_lock(&global_mpi_lock_queue);
        int tmp = --(*count);
        pthread_mutex_unlock(&global_mpi_lock_queue);
        return tmp;
    }
#endif


/************************************************************************
*  Empty constructor                                                    *
************************************************************************/
AMP_MPI::AMP_MPI() 
{
    // Initialize the data members to a defaul communicator of self
    #ifdef USE_MPI
        communicator = MPI_COMM_NULL;
        d_maxTag = 0x7FFFFFFF;
    #else
        communicator = AMP_COMM_NULL;
        d_maxTag = mpi_max_tag;
    #endif
    d_count = NULL;
    comm_rank = 0;
    comm_size = 1;
    d_isNull = true;
    d_currentTag = NULL;
    call_abort_in_serial_instead_of_exit = true;
    tmp_allignment = -1;
}


/************************************************************************
*  Empty deconstructor                                                  *
************************************************************************/
AMP_MPI::~AMP_MPI() 
{
    // Decrement the count if used
    int count = -1;
    if ( d_count!=NULL )
        count = decrement_count( d_count );
    if ( count == 0 ) {
        // We are holding that last reference to the MPI_Comm object, we need to free it
        #ifdef USE_MPI
            delete d_count;
            d_count = NULL;
            int err = MPI_Comm_free(&communicator);
            if ( err != MPI_SUCCESS )
                AMP_ERROR("Problem free'ing MPI_Comm object");
            communicator = AMP_COMM_NULL;
            ++N_MPI_Comm_destroyed;
        #else
            AMP_ERROR("Internal Error (why do we have a count in serial)");
        #endif
    }
    if ( d_currentTag==NULL ) {
        // No tag index
    } else if ( d_currentTag[1] > 1 ) {
        --(d_currentTag[1]);
    } else {
        delete [] d_currentTag;
    }
    comm_rank = 0;
    comm_size = 1;
    d_maxTag = 0;
    d_isNull = true;
    call_abort_in_serial_instead_of_exit = true;
}


/************************************************************************
*  Copy constructor                                                     *
************************************************************************/
AMP_MPI::AMP_MPI( const AMP::AMP_MPI& comm ) 
{
    // Initialize the data members to the existing comm object
    communicator = comm.communicator;
    comm_rank = comm.comm_rank;
    comm_size = comm.comm_size;
    d_isNull = comm.d_isNull;
    d_maxTag = comm.d_maxTag;
    d_currentTag = comm.d_currentTag;
    if ( d_currentTag != NULL )
        ++d_currentTag[1];
    call_abort_in_serial_instead_of_exit = comm.call_abort_in_serial_instead_of_exit;
    // Set and increment the count
    d_count = comm.d_count;
    if ( d_count != NULL )
        increment_count(d_count);
    tmp_allignment = -1;
}


/************************************************************************
*  Assignment operator                                                  *
************************************************************************/
AMP_MPI& AMP_MPI::operator=(const AMP::AMP_MPI& comm) 
{
    if (this == &comm) // protect against invalid self-assignment
        return *this;
    // Destroy the previous object
    this->~AMP_MPI();
    // Initialize the data members to the existing object
    this->communicator = comm.communicator;
    this->comm_rank = comm.comm_rank;
    this->comm_size = comm.comm_size;
    this->d_isNull = comm.d_isNull;
    this->d_maxTag = comm.d_maxTag;
    this->call_abort_in_serial_instead_of_exit = comm.call_abort_in_serial_instead_of_exit;
    this->d_currentTag = comm.d_currentTag;
    if ( this->d_currentTag != NULL )
        ++(this->d_currentTag[1]);
    // Set and increment the count
    this->d_count = comm.d_count;
    if ( this->d_count != NULL )
        increment_count(this->d_count);
    this->tmp_allignment = -1;
    return *this;
}


/************************************************************************
*  Constructor from existing MPI communicator                           *
************************************************************************/
AMP_MPI::AMP_MPI( MPI_Comm comm ) 
{
    #ifdef USE_MPI
        // We are using MPI, use the MPI communicator to initialize the data
        if ( comm==AMP_COMM_WORLD ) {
            communicator = AMP::AMPManager::comm_world.communicator;
        } else if ( comm==AMP_COMM_SELF ) {
            communicator = MPI_COMM_SELF;
        } else if ( comm==AMP_COMM_NULL ) {
            communicator = MPI_COMM_NULL;
        } else {
            communicator = comm;
        }
        if ( communicator!=MPI_COMM_NULL) {
            // Set the MPI_SIZE_T datatype if it has not been set
            if ( MPI_SIZE_T==0x0 )
                MPI_SIZE_T = getSizeTDataType();
            // Attach the error handler
            MPI_Comm_set_errhandler( communicator, AMP::AMPManager::mpierr );
            // Get the communicator properties
            MPI_Comm_rank(communicator, &comm_rank);
            MPI_Comm_size(communicator, &comm_size);
            int flag, *val;
            int ierr = MPI_Comm_get_attr(communicator,MPI_TAG_UB,&val,&flag);
            AMP_ASSERT(ierr==MPI_SUCCESS);
            if ( flag==0 ) { 
                d_maxTag = 0x7FFFFFFF;     // The tag is not a valid attribute (set to 2^31-1)
            } else {
                d_maxTag = *val;
                if ( d_maxTag<0 ) { d_maxTag = 0x7FFFFFFF; }    // The maximum tag is > a signed int (set to 2^31-1)
                AMP_INSIST(d_maxTag>=0x7FFF,"maximum tag size is < MPI standard");
            }
        } else {
            comm_rank = 1;
            comm_size = 0;
            d_maxTag = 0x7FFFFFFF;
        }
        d_isNull = communicator==MPI_COMM_NULL;
    #else
        // We are not using MPI, intialize based on the communicator
        communicator = comm;
        comm_rank = 0;
        comm_size = 1;
        d_maxTag = mpi_max_tag;
        d_isNull = communicator==AMP_COMM_NULL;
        if ( d_isNull )
            comm_size = 0;
    #endif
    if ( d_isNull ) {
        d_currentTag = NULL;
    } else {
        d_currentTag = new int[2];
        d_currentTag[0] = (d_maxTag<=0x10000) ? 1:0x1FFF;
        d_currentTag[1] = 1;
    }
    call_abort_in_serial_instead_of_exit = true;
    // We are creating a comm object from an MPI_Comm, the user is responsible for freeing the MPI_Comm object
    d_count = NULL;
    tmp_allignment = -1;
}


/************************************************************************
*  Intersect two communicators                                          *
************************************************************************/
#ifdef USE_MPI
AMP_MPI AMP_MPI::intersect( const AMP_MPI &comm1, const AMP_MPI &comm2 ) 
{
    MPI_Group group1=MPI_GROUP_EMPTY, group2=MPI_GROUP_EMPTY, group12=MPI_GROUP_EMPTY;
    if ( !comm1.isNull() )
        MPI_Comm_group ( comm1.communicator, &group1 );
    if ( !comm2.isNull() )
        MPI_Comm_group ( comm2.communicator, &group2 );
    MPI_Group_intersection( group1, group2, &group12 );
    int compare1, compare2;
    MPI_Group_compare ( group1, group12, &compare1 );
    MPI_Group_compare ( group2, group12, &compare2 );
    AMP_MPI new_comm(AMP_COMM_NULL);
    int size;
    MPI_Group_size( group12, &size );
    if ( compare1!=MPI_UNEQUAL && size!=0 ) {
        // The intersection matches comm1
        new_comm = comm1;
    } else if ( compare2!=MPI_UNEQUAL && size!=0 ) {
        // The intersection matches comm2
        new_comm = comm2;
    } else if ( comm1.isNull() ) {
        // comm1 is null, we can return safely (comm1 is needed for communication)
    } else {
        // The intersection is smaller than comm1 or comm2
        // Check if the new comm is NULL for all processors
        int max_size=0;
        MPI_Allreduce(&size,&max_size,1,MPI_INT,MPI_MAX,comm1.communicator);
        if ( max_size==0 ) {
            // We are dealing with completely disjoint sets
            new_comm = AMP_MPI( AMP_COMM_NULL );
        } else {
            // Create the new comm
            // Note: OpenMPI crashes if the intersection group is EMPTY for any processors
            // We will set it to SELF for the EMPTY processors, then create a NULL comm later
            if ( group12==MPI_GROUP_EMPTY )
                MPI_Comm_group( MPI_COMM_SELF, &group12 );
            MPI_Comm  new_MPI_comm;
            MPI_Comm_create( comm1.communicator, group12, &new_MPI_comm );
            if ( size>0 ) {
                // This is the valid case were we create a new intersection comm
                new_comm = AMP_MPI( new_MPI_comm );
            } else {
                // We actually want a null comm for this communicator
                new_comm = AMP_MPI( AMP_COMM_NULL );
                MPI_Comm_free(&new_MPI_comm);
            }
        }
    }
    if ( group1!=MPI_GROUP_NULL && group1!=MPI_GROUP_EMPTY )
        MPI_Group_free( &group1 );
    if ( group2!=MPI_GROUP_NULL && group2!=MPI_GROUP_EMPTY )
        MPI_Group_free( &group2 );
    if ( group12!=MPI_GROUP_NULL && group12!=MPI_GROUP_EMPTY )
        MPI_Group_free( &group12 );
    return new_comm;
}
#else
AMP_MPI AMP_MPI::intersect( const AMP_MPI &comm1, const AMP_MPI &comm2 ) 
{
    if ( comm1.isNull() || comm2.isNull() )
        return AMP_MPI(AMP_COMM_NULL);
    AMP_ASSERT(comm1.comm_size==1&&comm2.comm_size==1);
    return comm1;
}
#endif


/************************************************************************
*  Split a comm						                                    *
************************************************************************/
AMP_MPI AMP_MPI::split( int color, int key ) const 
{
    MPI_Comm  new_MPI_comm;
    #ifdef USE_MPI
        // USE MPI to split the communicator
        if ( color==-1 ) {
            MPI_Comm_split(communicator,MPI_UNDEFINED,key,&new_MPI_comm);
        } else {
            MPI_Comm_split(communicator,color,key,&new_MPI_comm);
        }
    #else
        if ( color==-1 ) {
            new_MPI_comm = AMP_COMM_NULL;
        } else {
            new_MPI_comm = uniqueGlobalComm;
            uniqueGlobalComm++;
        }
    #endif
    // Create the new object
    AMP_MPI new_comm(new_MPI_comm);
    new_comm.call_abort_in_serial_instead_of_exit = call_abort_in_serial_instead_of_exit;
    // Create the count (Note: we do not need to worry about thread safety)
    #ifdef USE_MPI
        if ( new_comm.communicator != MPI_COMM_NULL ) {
            new_comm.d_count = new int;
            *(new_comm.d_count) = 1;
            ++N_MPI_Comm_created;
        }
    #endif
    return new_comm;
}


/************************************************************************
*  Duplicate an exisiting comm object                                   *
************************************************************************/
AMP_MPI AMP_MPI::dup( ) const 
{
    if ( d_isNull )
        return AMP_MPI(AMP_COMM_NULL);
    MPI_Comm  new_MPI_comm;
    #ifdef USE_MPI
        // USE MPI to duplicate the communicator
        MPI_Comm_dup(communicator,&new_MPI_comm);
    #else
        new_MPI_comm = uniqueGlobalComm;
        uniqueGlobalComm++;
    #endif
    // Create the new comm object
    AMP_MPI new_comm(new_MPI_comm);
    new_comm.d_isNull = d_isNull;
    new_comm.call_abort_in_serial_instead_of_exit = call_abort_in_serial_instead_of_exit;
    // Create the count (Note: we do not need to worry about thread safety)
    #ifdef USE_MPI
        if ( new_comm.communicator != MPI_COMM_NULL ) {
            new_comm.d_count = new int;
            *(new_comm.d_count) = 1;
            ++N_MPI_Comm_created;
        }
    #endif
    return new_comm;
}


/************************************************************************
*  Overload operator ==                                                 *
************************************************************************/
bool AMP_MPI::operator==(const AMP_MPI &comm) const 
{
    return communicator==comm.communicator;
}


/************************************************************************
*  Overload operator !=                                                 *
************************************************************************/
bool AMP_MPI::operator!=(const AMP_MPI &comm) const 
{
    return communicator!=comm.communicator;
}


/************************************************************************
*  Overload operator <                                                  *
************************************************************************/
bool AMP_MPI::operator<(const AMP_MPI &comm) const 
{
    AMP_ASSERT( !this->d_isNull && !comm.d_isNull );
    bool flag = true;
    // First check if either communicator is NULL
    #ifdef USE_MPI
        if ( communicator==MPI_COMM_NULL )
            return false;
        if ( comm.communicator==MPI_COMM_NULL )
            flag = false;
    #else
        if ( communicator==AMP_COMM_NULL )
            return false;
        if ( comm.communicator==AMP_COMM_NULL )
            flag = false;
    #endif
    // Check that the size of the other communicator is > the current communicator size
    if ( comm_size >= comm.comm_size )
        flag = false;
    // Check the union of the communicator groups
    // this is < comm iff this group is a subgroup of comm's group
    #ifdef USE_MPI
        MPI_Group group1, group2, group12;
        MPI_Comm_group ( communicator, &group1 );
        MPI_Comm_group ( comm.communicator, &group2 );
        MPI_Group_union ( group1, group2, &group12 );
        int compare;
        MPI_Group_compare ( group2, group12, &compare );
        if ( compare==MPI_UNEQUAL )
            flag = false;
        MPI_Group_free( &group1 );
        MPI_Group_free( &group2 );
        MPI_Group_free( &group12 );
    #endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}


/************************************************************************
*  Overload operator <=                                                 *
************************************************************************/
bool AMP_MPI::operator<=(const AMP_MPI &comm) const 
{
    AMP_ASSERT( !this->d_isNull && !comm.d_isNull );
    bool flag = true;
    // First check if either communicator is NULL
    #ifdef USE_MPI
        if ( communicator==MPI_COMM_NULL )
            return false;
        if ( comm.communicator==MPI_COMM_NULL )
            flag = false;
    #else
        if ( communicator==AMP_COMM_NULL )
            return false;
        if ( comm.communicator==AMP_COMM_NULL )
            flag = false;
    #endif
    // Check that the size of the other communicator is > the current communicator size
    // this is <= comm iff this group is a subgroup of comm's group
    if ( comm_size > comm.comm_size )
        flag = false;
    // Check the unnion of the communicator groups
    #ifdef USE_MPI
        MPI_Group group1, group2, group12;
        MPI_Comm_group ( communicator, &group1 );
        MPI_Comm_group ( comm.communicator, &group2 );
        MPI_Group_union ( group1, group2, &group12 );
        int compare;
        MPI_Group_compare ( group2, group12, &compare );
        if ( compare==MPI_UNEQUAL )
            flag = false;
        MPI_Group_free( &group1 );
        MPI_Group_free( &group2 );
        MPI_Group_free( &group12 );
    #endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}


/************************************************************************
*  Overload operator >                                                  *
************************************************************************/
bool AMP_MPI::operator>(const AMP_MPI &comm) const 
{
    bool flag = true;
    // First check if either communicator is NULL
    #ifdef USE_MPI
        if ( communicator==MPI_COMM_NULL )
            return false;
        if ( comm.communicator==MPI_COMM_NULL )
            flag = false;
    #else
        if ( communicator==AMP_COMM_NULL )
            return false;
        if ( comm.communicator==AMP_COMM_NULL )
            flag = false;
    #endif
    // Check that the size of the other communicator is > the current communicator size
    if ( comm_size <= comm.comm_size )
        flag = false;
    // Check the unnion of the communicator groups
    // this is > comm iff comm's group is a subgroup of this group
    #ifdef USE_MPI
        MPI_Group group1, group2, group12;
        MPI_Comm_group ( communicator, &group1 );
        MPI_Comm_group ( comm.communicator, &group2 );
        MPI_Group_union ( group1, group2, &group12 );
        int compare;
        MPI_Group_compare ( group1, group12, &compare );
        if ( compare==MPI_UNEQUAL )
            flag = false;
        MPI_Group_free( &group1 );
        MPI_Group_free( &group2 );
        MPI_Group_free( &group12 );
    #endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}


/************************************************************************
*  Overload operator >=                                                 *
************************************************************************/
bool AMP_MPI::operator>=(const AMP_MPI &comm) const 
{
    bool flag = true;
    // First check if either communicator is NULL
    #ifdef USE_MPI
        if ( communicator==MPI_COMM_NULL )
            return false;
        if ( comm.communicator==MPI_COMM_NULL )
            flag = false;
    #else
        if ( communicator==AMP_COMM_NULL )
            return false;
        if ( comm.communicator==AMP_COMM_NULL )
            flag = false;
    #endif
    // Check that the size of the other communicator is > the current communicator size
    if ( comm_size < comm.comm_size )
        flag = false;
    // Check the unnion of the communicator groups
    // this is >= comm iff comm's group is a subgroup of this group
    #ifdef USE_MPI
        MPI_Group group1, group2, group12;
        MPI_Comm_group ( communicator, &group1 );
        MPI_Comm_group ( comm.communicator, &group2 );
        MPI_Group_union ( group1, group2, &group12 );
        int compare;
        MPI_Group_compare ( group1, group12, &compare );
        if ( compare==MPI_UNEQUAL )
            flag = false;
        MPI_Group_free( &group1 );
        MPI_Group_free( &group2 );
        MPI_Group_free( &group12 );
    #endif
    // Perform a global reduce of the flag (equivalent to all operation)
    return allReduce(flag);
}


/************************************************************************
*  Compare two comm objects                                             *
************************************************************************/
int AMP_MPI::compare(const AMP_MPI &comm) const 
{
    if ( communicator==comm.communicator )
        return 1;
    #ifdef USE_MPI
        int result;
        MPI_Comm_compare( communicator, comm.communicator, &result );
        if ( result==MPI_IDENT )
            return 2;
        else if ( result==MPI_CONGRUENT )
            return 3;
        else if ( result==MPI_SIMILAR )
            return 4;
        else if ( result==MPI_UNEQUAL )
            return 0;
        AMP_ERROR("Unknown results from AMP_Comm_compare");
    #else
        if ( comm.communicator==AMP_COMM_NULL || communicator==AMP_COMM_NULL )
            return 0;
        else 
            return 3;
    #endif
    return 0;
}


/************************************************************************
*  Abort the program.                                                   *
************************************************************************/
void AMP_MPI::setCallAbortInSerialInsteadOfExit(bool flag)
{
    call_abort_in_serial_instead_of_exit = flag;
}
void AMP_MPI::abort() const
{
    #ifdef USE_MPI
        int initialized=0, finalized=0;
        MPI_Initialized(&initialized);
        MPI_Finalized(&finalized);
        MPI_Comm comm = communicator;
        if ( comm == MPI_COMM_NULL )
            comm = MPI_COMM_WORLD;
        if ( initialized==0 || finalized!=0 ) {
            // MPI is not availible
            exit(-1);
        } else if ( comm_size > 1 ) {
            MPI_Abort(comm, -1);
        } else if ( call_abort_in_serial_instead_of_exit ) {
            MPI_Abort(comm, -1);
        } else {
            exit(-1);
        }
    #else
        exit(-1);
    #endif
}


/************************************************************************
*  newTag                                                               *
************************************************************************/
int AMP_MPI::newTag() 
{
    // Syncronize the processes to ensure all ranks enter this call 
    // Needed so the count will match
    barrier();
    // Return and increment the tag
    int tag = (*d_currentTag)++;
    AMP_INSIST(tag<=d_maxTag,"Maximum number of tags exceeded\n");
    return tag;
}


/************************************************************************
*  allReduce                                                            *
************************************************************************/
bool AMP_MPI::allReduce(const bool value) const 
{
    bool ret = value;
    if ( comm_size > 1 ) {
        #ifdef USE_MPI
            MPI_Allreduce( (void*) &value, (void*) &ret, 1, MPI_UNSIGNED_CHAR, MPI_MIN, communicator);
        #else
            AMP_ERROR("This shouldn't be possible");
        #endif
    }
    return ret;
}


/************************************************************************
*  anyReduce                                                            *
************************************************************************/
bool AMP_MPI::anyReduce(const bool value) const 
{
    bool ret = value;
    if ( comm_size > 1 ) {
        #ifdef USE_MPI
            MPI_Allreduce( (void*) &value, (void*) &ret, 1, MPI_UNSIGNED_CHAR, MPI_MAX, communicator);
        #else
            AMP_ERROR("This shouldn't be possible");
        #endif
    }
    return ret;
}


/************************************************************************
*  call_sumReduce                                                       *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void AMP_MPI::call_sumReduce<unsigned char>(const unsigned char *send, unsigned char *recv, const int n) const 
{
    PROFILE_START("sumReduce1<unsigned char>",profile_level);
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_SUM, communicator);
    PROFILE_STOP("sumReduce1<unsigned char>",profile_level);
}
template <>
void AMP_MPI::call_sumReduce<unsigned char>(unsigned char *x, const int n) const 
{
    PROFILE_START("sumReduce2<unsigned char>",profile_level);
    unsigned char *send = x;
    unsigned char *recv = new unsigned char[n];
    MPI_Allreduce( send, recv, n, MPI_UNSIGNED_CHAR, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
    PROFILE_STOP("sumReduce2<unsigned char>",profile_level);
}
// char
template <>
void AMP_MPI::call_sumReduce<char>(const char *send, char *recv, const int n) const 
{
    PROFILE_START("sumReduce1<char>",profile_level);
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_SUM, communicator);
    PROFILE_STOP("sumReduce1<char>",profile_level);
}
template <>
void AMP_MPI::call_sumReduce<char>(char *x, const int n) const 
{
    PROFILE_START("sumReduce2<char>",profile_level);
    char *send = x;
    char *recv = new char[n];
    MPI_Allreduce( send, recv, n, MPI_SIGNED_CHAR, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
    PROFILE_STOP("sumReduce2<char>",profile_level);
}
// unsigned int
template <>
void AMP_MPI::call_sumReduce<unsigned int>(const unsigned int *send, unsigned int *recv, const int n) const 
{
    PROFILE_START("sumReduce1<unsigned int>",profile_level);
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_SUM, communicator);
    PROFILE_STOP("sumReduce1<unsigned int>",profile_level);
}
template <>
void AMP_MPI::call_sumReduce<unsigned int>(unsigned int *x, const int n) const 
{
    PROFILE_START("sumReduce2<unsigned int>",profile_level);
    unsigned int *send = x;
    unsigned int *recv = new unsigned int[n];
    MPI_Allreduce( send, recv, n, MPI_UNSIGNED, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
    PROFILE_STOP("sumReduce2<unsigned int>",profile_level);
}
// int
template <>
void AMP_MPI::call_sumReduce<int>(const int *send, int *recv, const int n) const 
{
    PROFILE_START("sumReduce1<int>",profile_level);
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_INT, MPI_SUM, communicator);
    PROFILE_STOP("sumReduce1<int>",profile_level);
}
template <>
void AMP_MPI::call_sumReduce<int>(int *x, const int n) const 
{
    PROFILE_START("sumReduce2<int>",profile_level);
    int *send = x;
    int *recv = new int[n];
    MPI_Allreduce( send, recv, n, MPI_INT, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
    PROFILE_STOP("sumReduce2<int>",profile_level);
}
// long int
template <>
void AMP_MPI::call_sumReduce<long int>(const long int *send, long int *recv, const int n) const 
{
    PROFILE_START("sumReduce1<long int>",profile_level);
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_LONG, MPI_SUM, communicator);
    PROFILE_STOP("sumReduce1<long int>",profile_level);
}
template <>
void AMP_MPI::call_sumReduce<long int>(long int *x, const int n) const 
{
    PROFILE_START("sumReduce2<long int>",profile_level);
    long int *send = x;
    long int *recv = new long int[n];
    MPI_Allreduce( send, recv, n, MPI_LONG, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
    PROFILE_STOP("sumReduce2<long int>",profile_level);
}
// unsigned long int
template <>
void AMP_MPI::call_sumReduce<unsigned long>(const unsigned long *send, unsigned long *recv, const int n) const 
{
    PROFILE_START("sumReduce1<unsigned long>",profile_level);
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_SUM, communicator);
    PROFILE_STOP("sumReduce1<unsigned long>",profile_level);
}
template <>
void AMP_MPI::call_sumReduce<unsigned long>(unsigned long *x, const int n) const 
{
    PROFILE_START("sumReduce2<unsigned long>",profile_level);
    unsigned long int *send = x;
    unsigned long int *recv = new unsigned long int[n];
    MPI_Allreduce( send, recv, n, MPI_UNSIGNED_LONG, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
    PROFILE_STOP("sumReduce2<unsigned long>",profile_level);
}
// size_t
#ifdef USE_WINDOWS
    template <>
    void AMP_MPI::call_sumReduce<size_t>(const size_t *send, size_t *recv, const int n) const 
    {
        PROFILE_START("sumReduce1<size_t>",profile_level);
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIZE_T, MPI_SUM, communicator);
        PROFILE_STOP("sumReduce1<size_t>",profile_level);
    }
    template <>
    void AMP_MPI::call_sumReduce<size_t>(size_t *x, const int n) const 
    {
        PROFILE_START("sumReduce2<size_t>",profile_level);
        size_t *send = x;
        size_t *recv = new size_t[n];
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIZE_T, MPI_SUM, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
        PROFILE_STOP("sumReduce2<size_t>",profile_level);
    }
#endif
// float
template <>
void AMP_MPI::call_sumReduce<float>(const float *send, float *recv, const int n) const 
{
    PROFILE_START("sumReduce1<float>",profile_level);
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_FLOAT, MPI_SUM, communicator);
    PROFILE_STOP("sumReduce1<float>",profile_level);
}
template <>
void AMP_MPI::call_sumReduce<float>(float *x, const int n) const 
{
    PROFILE_START("sumReduce2<float>",profile_level);
    float *send = x;
    float *recv = new float[n];
    MPI_Allreduce( send, recv, n, MPI_FLOAT, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
    PROFILE_STOP("sumReduce2<float>",profile_level);
}
// double
template <>
void AMP_MPI::call_sumReduce<double>(const double *send, double *recv, const int n) const 
{
    PROFILE_START("sumReduce1<double>",profile_level);
    MPI_Allreduce( (void*) send, (void*) recv, n, MPI_DOUBLE, MPI_SUM, communicator);
    PROFILE_STOP("sumReduce1<double>",profile_level);
}
template <>
void AMP_MPI::call_sumReduce<double>(double *x, const int n) const 
{
    PROFILE_START("sumReduce2<double>",profile_level);
    double *send = x;
    double *recv = new double[n];
    MPI_Allreduce( send, recv, n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = recv[i];
    delete [] recv;
    PROFILE_STOP("sumReduce2<double>",profile_level);
}
// std::complex<double>
template <>
void AMP_MPI::call_sumReduce< std::complex<double> >(const std::complex<double> *x, std::complex<double> *y, const int n) const 
{
    PROFILE_START("sumReduce1<complex double>",profile_level);
    double *send = new double[2*n];
    double *recv = new double[2*n];
    for (int i=0; i<n; i++) {
        send[2*i+0] = real(x[i]);
        send[2*i+1] = imag(x[i]);
    }
    MPI_Allreduce( (void*) send, (void*) recv, 2*n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        y[i] = std::complex<double>(recv[2*i+0],recv[2*i+1]);
    delete [] send;
    delete [] recv;
    PROFILE_STOP("sumReduce1<complex double>",profile_level);
}
template <>
void AMP_MPI::call_sumReduce< std::complex<double> >(std::complex<double> *x, const int n) const 
{
    PROFILE_START("sumReduce2<complex double>",profile_level);
    double *send = new double[2*n];
    double *recv = new double[2*n];
    for (int i=0; i<n; i++) {
        send[2*i+0] = real(x[i]);
        send[2*i+1] = imag(x[i]);
    }
    MPI_Allreduce( send, recv, 2*n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        x[i] = std::complex<double>(recv[2*i+0],recv[2*i+1]);
    delete [] send;
    delete [] recv;
    PROFILE_STOP("sumReduce2<complex double>",profile_level);
}
#endif


/************************************************************************
*  call_minReduce                                                       *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void AMP_MPI::call_minReduce<unsigned char>(const unsigned char *send, unsigned char *recv, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("minReduce1<unsigned char>",profile_level);
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_MIN, communicator);
        PROFILE_STOP("minReduce1<unsigned char>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}
template <>
void AMP_MPI::call_minReduce<unsigned char>(unsigned char *x, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("minReduce2<unsigned char>",profile_level);
        unsigned char *send = x;
        unsigned char *recv = new unsigned char[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED_CHAR, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
        PROFILE_STOP("minReduce2<unsigned char>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}
// char
template <>
void AMP_MPI::call_minReduce<char>(const char *send, char *recv, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("minReduce1<char>",profile_level);
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_MIN, communicator);
        PROFILE_STOP("minReduce1<char>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with char is not supported yet");
    }
}
template <>
void AMP_MPI::call_minReduce<char>(char *x, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("minReduce2<char>",profile_level);
        char *send = x;
        char *recv = new char[n];
        MPI_Allreduce( send, recv, n, MPI_SIGNED_CHAR, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
        PROFILE_STOP("minReduce2<char>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with char is not supported yet");
    }
}
// unsigned int
template <>
void AMP_MPI::call_minReduce<unsigned int>(const unsigned int *send, unsigned int *recv, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("minReduce1<unsigned int>",profile_level);
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_MIN, communicator);
        PROFILE_STOP("minReduce1<unsigned int>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}
template <>
void AMP_MPI::call_minReduce<unsigned int>(unsigned int *x, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("minReduce2<unsigned int>",profile_level);
        unsigned int *send = x;
        unsigned int *recv = new unsigned int[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
        PROFILE_STOP("minReduce2<unsigned int>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned int is not supported yet");
    }
}
// int
template <>
void AMP_MPI::call_minReduce<int>(const int *x, int *y, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("minReduce1<int>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_INT, MPI_MIN, communicator);
    } else {
        IntIntStruct *recv = new IntIntStruct[n];
        IntIntStruct *send = new IntIntStruct[n];
        for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
        }
        MPI_Allreduce( send, recv, n, MPI_2INT, MPI_MINLOC, communicator);
        for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
        }
        delete [] recv;
        delete [] send;
    }
    PROFILE_STOP("minReduce1<int>",profile_level);
}
template <>
void AMP_MPI::call_minReduce<int>(int *x, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("minReduce2<int>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        int *send = x;
        int *recv = new int[n];
        MPI_Allreduce( send, recv, n, MPI_INT, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_2INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("minReduce2<int>",profile_level);
}
// unsigned long int
template <>
void AMP_MPI::call_minReduce<unsigned long int>(const unsigned long int *send, unsigned long int *recv, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("minReduce1<unsigned long>",profile_level);
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_MIN, communicator);
        PROFILE_STOP("minReduce1<unsigned long>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}
template <>
void AMP_MPI::call_minReduce<unsigned long int>(unsigned long int *x, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("minReduce2<unsigned long>",profile_level);
        unsigned long int *send = x;
        unsigned long int *recv = new unsigned long int[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED_LONG, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
        PROFILE_STOP("minReduce2<unsigned long>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned long int is not supported yet");
    }
}
// long int
template <>
void AMP_MPI::call_minReduce<long int>(const long int *x, long int *y, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("minReduce1<long int>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_LONG, MPI_MIN, communicator);
    } else {
         LongIntStruct *recv = new LongIntStruct[n];
         LongIntStruct *send = new LongIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_LONG_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("minReduce1<long int>",profile_level);
}
template <>
void AMP_MPI::call_minReduce<long int>(long int *x, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("minReduce2<long int>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        long int *send = x;
        long int *recv = new long int[n];
        MPI_Allreduce( send, recv, n, MPI_LONG, MPI_MIN, communicator);
        for (long int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         LongIntStruct *recv = new LongIntStruct[n];
         LongIntStruct *send = new LongIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_LONG_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("minReduce2<long int>",profile_level);
}
// size_t
#ifdef USE_WINDOWS
    template <>
    void AMP_MPI::call_minReduce<size_t>(const size_t *send, size_t *recv, const int n, int *comm_rank_of_min) const 
    {
        if ( comm_rank_of_min==NULL ) {
            PROFILE_START("minReduce1<size_t>",profile_level);
            MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIZE_T, MPI_MIN, communicator);
            PROFILE_STOP("minReduce1<size_t>",profile_level);
        } else {
             AMP_ERROR("Returning the rank of min with size_t is not supported yet");
        }
    }
    template <>
    void AMP_MPI::call_minReduce<size_t>(size_t *x, const int n, int *comm_rank_of_min) const 
    {
        if ( comm_rank_of_min==NULL ) {
            PROFILE_START("minReduce2<size_t>",profile_level);
            size_t *send = x;
            size_t *recv = new size_t[n];
            MPI_Allreduce( send, recv, n, MPI_SIZE_T, MPI_MIN, communicator);
            for (int i=0; i<n; i++)
                x[i] = recv[i];
            delete [] recv;
            PROFILE_STOP("minReduce2<size_t>",profile_level);
        } else {
             AMP_ERROR("Returning the rank of min with size_t is not supported yet");
        }
    }
#endif
// float
template <>
void AMP_MPI::call_minReduce<float>(const float *x, float *y, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("minReduce1<float>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_INT, MPI_MIN, communicator);
    } else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_FLOAT_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].f;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("minReduce1<float>",profile_level);
}
template <>
void AMP_MPI::call_minReduce<float>(float *x, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("minReduce2<float>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        float *send = x;
        float *recv = new float[n];
        MPI_Allreduce( send, recv, n, MPI_FLOAT, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_FLOAT_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].f;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("minReduce2<float>",profile_level);
}
// double
template <>
void AMP_MPI::call_minReduce<double>(const double *x, double *y, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("minReduce1<double>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_DOUBLE, MPI_MIN, communicator);
    } else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_DOUBLE_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].d;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("minReduce1<double>",profile_level);
}
template <>
void AMP_MPI::call_minReduce<double>(double *x, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("minReduce2<double>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        double *send = x;
        double *recv = new double[n];
        MPI_Allreduce( send, recv, n, MPI_DOUBLE, MPI_MIN, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_DOUBLE_INT, MPI_MINLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].d;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("minReduce2<double>",profile_level);
}
#endif


/************************************************************************
*  call_maxReduce                                                    *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void AMP_MPI::call_maxReduce<unsigned char>(const unsigned char *send, unsigned char *recv, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("maxReduce1<unsigned char>",profile_level);
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_MAX, communicator);
        PROFILE_STOP("maxReduce1<unsigned char>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}template <>
void AMP_MPI::call_maxReduce<unsigned char>(unsigned char *x, const int n, int *comm_rank_of_max) const 
{
    if ( comm_rank_of_max==NULL ) {
        PROFILE_START("maxReduce2<unsigned char>",profile_level);
        unsigned char *send = x;
        unsigned char *recv = new unsigned char[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED_CHAR, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
        PROFILE_STOP("maxReduce2<unsigned char>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}
// char
template <>
void AMP_MPI::call_maxReduce<char>(const char *send, char *recv, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("maxReduce1<char>",profile_level);
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_MAX, communicator);
        PROFILE_STOP("maxReduce1<char>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with char is not supported yet");
    }
}
template <>
void AMP_MPI::call_maxReduce<char>(char *x, const int n, int *comm_rank_of_max) const 
{
    if ( comm_rank_of_max==NULL ) {
        PROFILE_START("maxReduce2<char>",profile_level);
        char *send = x;
        char *recv = new char[n];
        MPI_Allreduce( send, recv, n, MPI_SIGNED_CHAR, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
        PROFILE_STOP("maxReduce2<char>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with char is not supported yet");
    }
}
// unsigned int
template <>
void AMP_MPI::call_maxReduce<unsigned int>(const unsigned int *send, unsigned int *recv, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("maxReduce1<unsigned int>",profile_level);
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_MAX, communicator);
        PROFILE_STOP("maxReduce1<unsigned int>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}
template <>
void AMP_MPI::call_maxReduce<unsigned int>(unsigned int *x, const int n, int *comm_rank_of_max) const 
{
    if ( comm_rank_of_max==NULL ) {
        PROFILE_START("maxReduce2<unsigned int>",profile_level);
        unsigned int *send = x;
        unsigned int *recv = new unsigned int[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
        PROFILE_STOP("maxReduce2<unsigned int>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned int is not supported yet");
    }
}
// int
template <>
void AMP_MPI::call_maxReduce<int>(const int *x, int *y, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("maxReduce1<int>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_INT, MPI_MAX, communicator);
    } else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_2INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("maxReduce1<int>",profile_level);
}
template <>
void AMP_MPI::call_maxReduce<int>(int *x, const int n, int *comm_rank_of_max) const 
{
    PROFILE_START("maxReduce2<int>",profile_level);
    if ( comm_rank_of_max==NULL ) {
        int *send = x;
        int *recv = new int[n];
        MPI_Allreduce( send, recv, n, MPI_INT, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         IntIntStruct *recv = new IntIntStruct[n];
         IntIntStruct *send = new IntIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_2INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            comm_rank_of_max[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("maxReduce2<int>",profile_level);
}
// long int
template <>
void AMP_MPI::call_maxReduce<long int>(const long int *x, long int *y, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("maxReduce1<lond int>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_LONG, MPI_MAX, communicator);
    } else {
         LongIntStruct *recv = new LongIntStruct[n];
         LongIntStruct *send = new LongIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_LONG_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].j;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("maxReduce1<lond int>",profile_level);
}
template <>
void AMP_MPI::call_maxReduce<long int>(long int *x, const int n, int *comm_rank_of_max) const 
{
    PROFILE_START("maxReduce2<lond int>",profile_level);
    if ( comm_rank_of_max==NULL ) {
        long int *send = x;
        long int *recv = new long int[n];
        MPI_Allreduce( send, recv, n, MPI_LONG, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         LongIntStruct *recv = new LongIntStruct[n];
         LongIntStruct *send = new LongIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].j = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_LONG_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].j;
            comm_rank_of_max[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("maxReduce2<lond int>",profile_level);
}
// unsigned long int
template <>
void AMP_MPI::call_maxReduce<unsigned long int>(const unsigned long int *send, unsigned long int *recv, const int n, int *comm_rank_of_min) const 
{
    if ( comm_rank_of_min==NULL ) {
        PROFILE_START("maxReduce1<unsigned long>",profile_level);
        MPI_Allreduce( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_MAX, communicator);
        PROFILE_STOP("maxReduce1<unsigned long>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
    }
}
template <>
void AMP_MPI::call_maxReduce<unsigned long int>(unsigned long int *x, const int n, int *comm_rank_of_max) const 
{
    if ( comm_rank_of_max==NULL ) {
        PROFILE_START("maxReduce2<unsigned long>",profile_level);
        unsigned long int *send = x;
        unsigned long int *recv = new unsigned long int[n];
        MPI_Allreduce( send, recv, n, MPI_UNSIGNED_LONG, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
        PROFILE_STOP("maxReduce2<unsigned long>",profile_level);
    } else {
         AMP_ERROR("Returning the rank of min with unsigned long int is not supported yet");
    }
}
// size_t
#ifdef USE_WINDOWS
    template <>
    void AMP_MPI::call_maxReduce<size_t>(const size_t *send, size_t *recv, const int n, int *comm_rank_of_max) const 
{
        if ( comm_rank_of_max==NULL ) {
            PROFILE_START("minReduce1<size_t>",profile_level);
            MPI_Allreduce( (void*) send, (void*) recv, n, MPI_SIZE_T, MPI_MAX, communicator);
            PROFILE_STOP("minReduce1<size_t>",profile_level);
        } else {
             AMP_ERROR("Returning the rank of min with unsigned char is not supported yet");
        }
    }
    template <>
    void AMP_MPI::call_maxReduce<size_t>(size_t *x, const int n, int *comm_rank_of_max) const 
{
        if ( comm_rank_of_max==NULL ) {
            PROFILE_START("minReduce2<size_t>",profile_level);
            size_t *send = x;
            size_t *recv = new size_t[n];
            MPI_Allreduce( send, recv, n, MPI_SIZE_T, MPI_MAX, communicator);
            for (int i=0; i<n; i++)
                x[i] = recv[i];
            delete [] recv;
            PROFILE_STOP("minReduce2<size_t>",profile_level);
        } else {
             AMP_ERROR("Returning the rank of min with unsigned int is not supported yet");
        }
    }
#endif
// float
template <>
void AMP_MPI::call_maxReduce<float>(const float *x, float *y, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("maxReduce1<float>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_FLOAT, MPI_MAX, communicator);
    } else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_FLOAT_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].f;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("maxReduce1<float>",profile_level);
}
template <>
void AMP_MPI::call_maxReduce<float>(float *x, const int n, int *comm_rank_of_max) const 
{
    PROFILE_START("maxReduce2<float>",profile_level);
    if ( comm_rank_of_max==NULL ) {
        float *send = x;
        float *recv = new float[n];
        MPI_Allreduce( send, recv, n, MPI_FLOAT, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         FloatIntStruct *recv = new FloatIntStruct[n];
         FloatIntStruct *send = new FloatIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].f = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_FLOAT_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].f;
            comm_rank_of_max[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("maxReduce2<float>",profile_level);
}
// double
template <>
void AMP_MPI::call_maxReduce<double>(const double *x, double *y, const int n, int *comm_rank_of_min) const 
{
    PROFILE_START("maxReduce1<double>",profile_level);
    if ( comm_rank_of_min==NULL ) {
        MPI_Allreduce( (void*) x, (void*) y, n, MPI_DOUBLE, MPI_MAX, communicator);
    } else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_DOUBLE_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            y[i] = recv[i].d;
            comm_rank_of_min[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("maxReduce1<double>",profile_level);
}
template <>
void AMP_MPI::call_maxReduce<double>(double *x, const int n, int *comm_rank_of_max) const 
{
    PROFILE_START("maxReduce2<double>",profile_level);
    if ( comm_rank_of_max==NULL ) {
        double *send = x;
        double *recv = new double[n];
        MPI_Allreduce( send, recv, n, MPI_DOUBLE, MPI_MAX, communicator);
        for (int i=0; i<n; i++)
            x[i] = recv[i];
        delete [] recv;
    } else {
         DoubleIntStruct *recv = new DoubleIntStruct[n];
         DoubleIntStruct *send = new DoubleIntStruct[n];
         for ( int i=0; i<n; ++i ) {
            send[i].d = x[i];
            send[i].i = comm_rank;
         }
         MPI_Allreduce( send, recv, n, MPI_DOUBLE_INT, MPI_MAXLOC, communicator);
         for ( int i=0; i<n; ++i ) {
            x[i] = recv[i].d;
            comm_rank_of_max[i] = recv[i].i;
         }
         delete [] recv;
         delete [] send;
    }
    PROFILE_STOP("maxReduce2<double>",profile_level);
}
#endif


/************************************************************************
*  bcast                                                                *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// char
template <>
void AMP_MPI::call_bcast<unsigned char>(unsigned char *x, const int n, const int root) const 
{
    PROFILE_START("bcast<char>",profile_level);
    MPI_Bcast( x, n, MPI_CHAR, root, communicator);
    PROFILE_STOP("bcast<char>",profile_level);
}
template <>
void AMP_MPI::call_bcast<char>(char *x, const int n, const int root) const 
{
    PROFILE_START("bcast<char>",profile_level);
    MPI_Bcast( x, n, MPI_CHAR, root, communicator);
    PROFILE_STOP("bcast<char>",profile_level);
}
// int
template <>
void AMP_MPI::call_bcast<unsigned int>(unsigned int *x, const int n, const int root) const 
{
    PROFILE_START("bcast<int>",profile_level);
    MPI_Bcast( x, n, MPI_INT, root, communicator);
    PROFILE_STOP("bcast<int>",profile_level);
}
template <>
void AMP_MPI::call_bcast<int>(int *x, const int n, const int root) const 
{
    PROFILE_START("bcast<int>",profile_level);
    MPI_Bcast( x, n, MPI_INT, root, communicator);
    PROFILE_STOP("bcast<int>",profile_level);
}
// float
template <>
void AMP_MPI::call_bcast<float>(float *x, const int n, const int root) const 
{
    PROFILE_START("bcast<float>",profile_level);
    MPI_Bcast( x, n, MPI_FLOAT, root, communicator);
    PROFILE_STOP("bcast<float>",profile_level);
}
// double
template <>
void AMP_MPI::call_bcast<double>(double *x, const int n, const int root) const 
{
    PROFILE_START("bcast<double>",profile_level);
    MPI_Bcast( x, n, MPI_DOUBLE, root, communicator);
    PROFILE_STOP("bcast<double>",profile_level);
}
#else
// We need a concrete instantiation of bcast<char>(x,n,root);
template <>
void AMP_MPI::call_bcast<char>(char *x, const int n, const int root) const 
{
    AMP_ERROR("Internal error in AMP_MPI (bcast) ");
}
#endif


/************************************************************************
*  Perform a global barrier across all processors.                      *
************************************************************************/
void AMP_MPI::barrier() const
{
    #ifdef USE_MPI
        MPI_Barrier(communicator);
    #endif
}


/************************************************************************
*  Send data array to another processor.                                *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// char
template <>
void AMP_MPI::send<char>(const char *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Send the data
    PROFILE_START("send<char>",profile_level);
    MPI_Send((void*)buf, length, MPI_CHAR, recv_proc_number, tag, communicator);
    PROFILE_STOP("send<char>",profile_level);
}
// int
template <>
void AMP_MPI::send<int>(const int *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Send the data 
    PROFILE_START("send<int>",profile_level);
    MPI_Send((void*)buf, length, MPI_INT, recv_proc_number, tag, communicator);
    PROFILE_STOP("send<int>",profile_level);
}
// float
template <>
void AMP_MPI::send<float>(const float *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Send the data 
    PROFILE_START("send<float>",profile_level);
    MPI_Send((void*)buf, length, MPI_FLOAT, recv_proc_number, tag, communicator);
    PROFILE_STOP("send<float>",profile_level);
}
// double
template <>
void AMP_MPI::send<double>(const double *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    // Send the data 
    PROFILE_START("send<double>",profile_level);
    MPI_Send((void*)buf, length, MPI_DOUBLE, recv_proc_number, tag, communicator);
    PROFILE_STOP("send<double>",profile_level);
}
#else
// We need a concrete instantiation of send for use without MPI
template <>
void AMP_MPI::send<char>(const char *buf, const int length, 
    const int recv_proc_number, int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    PROFILE_START("send<char>",profile_level);
    MPI_Request id = getRequest( communicator, tag );
    std::map<MPI_Request,Isendrecv_struct>::iterator it = global_isendrecv_list.find(id);
    AMP_INSIST(it==global_isendrecv_list.end(),"send must be paired with a previous call to irecv in serial");
    AMP_ASSERT(it->second.status==2);
    memcpy((char*)it->second.data,buf,length);
    global_isendrecv_list.erase( it );
    PROFILE_START("send<char>",profile_level);
}
#endif


/************************************************************************
*  Non-blocking send data array to another processor.                   *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// char
template <>
MPI_Request AMP_MPI::Isend<char>(const char *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Isend<char>",profile_level);
    MPI_Isend((void*)buf, length, MPI_CHAR, recv_proc, tag, communicator, &request);
    PROFILE_STOP("Isend<char>",profile_level);
    return request;
}
// int
template <>
MPI_Request AMP_MPI::Isend<int>(const int *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Isend<int>",profile_level);
    MPI_Isend((void*)buf, length, MPI_INT, recv_proc, tag, communicator, &request);
    PROFILE_STOP("Isend<int>",profile_level);
    return request;
}
// float
template <>
MPI_Request AMP_MPI::Isend<float>(const float *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Isend<float>",profile_level);
    MPI_Isend((void*)buf, length, MPI_FLOAT, recv_proc, tag, communicator, &request);
    PROFILE_STOP("Isend<float>",profile_level);
    return request;
}
// double
template <>
MPI_Request AMP_MPI::Isend<double>(const double *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Isend<double>",profile_level);
    MPI_Isend((void*)buf, length, MPI_DOUBLE, recv_proc, tag, communicator, &request);
    PROFILE_STOP("Isend<double>",profile_level);
    return request;
}
#else
// We need a concrete instantiation of send for use without mpi
template <>
MPI_Request AMP_MPI::Isend<char>(const char *buf, const int length, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    PROFILE_START("Isend<char>",profile_level);
    MPI_Request id = getRequest( communicator, tag );
    std::map<MPI_Request,Isendrecv_struct>::iterator it = global_isendrecv_list.find(id);
    if ( it==global_isendrecv_list.end() ) {
        // We are calling isend first
        Isendrecv_struct data;
        data.data = buf;
        data.status = 1;
        global_isendrecv_list.insert( std::pair<MPI_Request,Isendrecv_struct>(id,data) );
    } else {
        // We called irecv first
        AMP_ASSERT(it->second.status==2);
        memcpy((char*)it->second.data,buf,length);
        global_isendrecv_list.erase( it );
    }
    PROFILE_STOP("Isend<char>",profile_level);
    return id;
}
#endif


/************************************************************************
*  Send byte array to another processor.                                *
************************************************************************/
void AMP_MPI::sendBytes(const void *buf, const int number_bytes, 
    const int recv_proc_number, int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    send<char>((const char*)buf,number_bytes,recv_proc_number,tag);
}


/************************************************************************
*  Non-blocking send byte array to another processor.                   *
************************************************************************/
MPI_Request AMP_MPI::IsendBytes(const void *buf, const int number_bytes, const int recv_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    return Isend<char>((const char*)buf,number_bytes,recv_proc,tag);
}


/************************************************************************
*  Recieve data array to another processor.                             *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// char
template <>
void AMP_MPI::recv<char>(char *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    PROFILE_START("recv<char>",profile_level);
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe( send_proc_number, tag );
        int recv_length = bytes/sizeof(char);
        AMP_INSIST(length>=recv_length,"Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data 
    MPI_Status status;
    MPI_Recv((void*)buf, length, MPI_CHAR, send_proc_number, tag, communicator, &status);
    PROFILE_STOP("recv<char>",profile_level);
}
// int
template <>
void AMP_MPI::recv<int>(int *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    PROFILE_START("recv<int>",profile_level);
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe( send_proc_number, tag );
        int recv_length = bytes/sizeof(int);
        AMP_INSIST(length>=recv_length,"Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data 
    MPI_Status status;
    MPI_Recv((void*)buf, length, MPI_INT, send_proc_number, tag, communicator, &status);
    PROFILE_STOP("recv<int>",profile_level);
}
// float
template <>
void AMP_MPI::recv<float>(float *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    PROFILE_START("recv<float>",profile_level);
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe( send_proc_number, tag );
        int recv_length = bytes/sizeof(float);
        AMP_INSIST(length>=recv_length,"Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data 
    MPI_Status status;
    MPI_Recv((void*)buf, length, MPI_FLOAT, send_proc_number, tag, communicator, &status);
    PROFILE_STOP("recv<float>",profile_level);
}
// double
template <>
void AMP_MPI::recv<double>(double *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    // Set the tag to 0 if it is < 0
    tag = (tag >= 0) ? tag : 0;
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    PROFILE_START("recv<double>",profile_level);
    // Get the recieve length if necessary
    if (get_length) {
        int bytes = this->probe( send_proc_number, tag );
        int recv_length = bytes/sizeof(double);
        AMP_INSIST(length>=recv_length,"Recived length is larger than allocated array");
        length = recv_length;
    }
    // Send the data 
    MPI_Status status;
    MPI_Recv((void*)buf, length, MPI_DOUBLE, send_proc_number, tag, communicator, &status);
    PROFILE_STOP("recv<double>",profile_level);
}
#else
// We need a concrete instantiation of recv for use without mpi
template <>
void AMP_MPI::recv<char>(char *buf, int &length, 
    const int send_proc_number, const bool get_length, int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    PROFILE_START("recv<char>",profile_level);
    MPI_Request id = getRequest( communicator, tag );
    std::map<MPI_Request,Isendrecv_struct>::iterator it = global_isendrecv_list.find(id);
    AMP_INSIST(it!=global_isendrecv_list.end(),"recv must be paired with a previous call to isend in serial");
    AMP_ASSERT(it->second.status==1);
    memcpy(buf,it->second.data,length);
    global_isendrecv_list.erase( it );
    PROFILE_STOP("recv<char>",profile_level);
}
#endif


/************************************************************************
*  Non-blocking recieve data array to another processor.                *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// char
template <>
MPI_Request AMP_MPI::Irecv<char>(char *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Irecv<char>",profile_level);
    MPI_Irecv((void*)buf, length, MPI_CHAR, send_proc, tag, communicator, &request);
    PROFILE_STOP("Irecv<char>",profile_level);
    return request;
}
// int
template <>
MPI_Request AMP_MPI::Irecv<int>(int *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Irecv<int>",profile_level);
    MPI_Irecv((void*)buf, length, MPI_INT, send_proc, tag, communicator, &request);
    PROFILE_STOP("Irecv<int>",profile_level);
    return request;
}
// float
template <>
MPI_Request AMP_MPI::Irecv<float>(float *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Irecv<float>",profile_level);
    MPI_Irecv((void*)buf, length, MPI_FLOAT, send_proc, tag, communicator, &request);
    PROFILE_STOP("Irecv<float>",profile_level);
    return request;
}
// double
template <>
MPI_Request AMP_MPI::Irecv<double>(double *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Request request;
    PROFILE_START("Irecv<double>",profile_level);
    MPI_Irecv((void*)buf, length, MPI_DOUBLE, send_proc, tag, communicator, &request);
    PROFILE_STOP("Irecv<double>",profile_level);
    return request;
}
#else
// We need a concrete instantiation of irecv for use without mpi
template <>
MPI_Request AMP_MPI::Irecv<char>(char *buf, const int length, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    PROFILE_START("Irecv<char>",profile_level);
    MPI_Request id = getRequest( communicator, tag );
    std::map<MPI_Request,Isendrecv_struct>::iterator it = global_isendrecv_list.find(id);
    if ( it==global_isendrecv_list.end() ) {
        // We are calling Irecv first
        Isendrecv_struct data;
        data.data = buf;
        data.status = 2;
        global_isendrecv_list.insert( std::pair<MPI_Request,Isendrecv_struct>(id,data) );
    } else {
        // We called Isend first
        AMP_ASSERT(it->second.status==1);
        memcpy(buf,it->second.data,length);
        global_isendrecv_list.erase( it );
    }
    PROFILE_STOP("Irecv<char>",profile_level);
    return id;
}
#endif


/************************************************************************
*  Recieve byte array to another processor.                             *
************************************************************************/
void AMP_MPI::recvBytes(void *buf, int &number_bytes, const int send_proc, int tag) const
{
    recv<char>((char*)buf,number_bytes,send_proc,false,tag);
}


/************************************************************************
*  Recieve byte array to another processor.                             *
************************************************************************/
MPI_Request AMP_MPI::IrecvBytes(void *buf, const int number_bytes, const int send_proc, const int tag) const
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    return Irecv<char>((char*)buf,number_bytes,send_proc,tag);
}


/************************************************************************
*  allGather                                                            *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void AMP_MPI::call_allGather<unsigned char>(const unsigned char x_in, unsigned char *x_out) const 
{
    PROFILE_START("allGather<unsigned char>",profile_level);
    MPI_Allgather( (void*) &x_in, 1, MPI_UNSIGNED_CHAR, (void*) x_out, 1, MPI_UNSIGNED_CHAR, communicator );
    PROFILE_STOP("allGather<unsigned char>",profile_level);
}
template <>
void AMP_MPI::call_allGather<unsigned char>(const unsigned char *x_in, int size_in, unsigned char *x_out, int *size_out, int *disp_out) const 
{
    PROFILE_START("allGatherv<unsigned char>",profile_level);
    MPI_Allgatherv( (void*) x_in, size_in, MPI_CHAR, (void*) x_out, size_out, disp_out, MPI_CHAR, communicator );
    PROFILE_STOP("allGatherv<unsigned char>",profile_level);
}
// char
template <>
void AMP_MPI::call_allGather<char>(const char x_in, char *x_out) const 
{
    PROFILE_START("allGather<char>",profile_level);
    MPI_Allgather( (void*) &x_in, 1, MPI_CHAR, (void*) x_out, 1, MPI_CHAR, communicator );
    PROFILE_STOP("allGather<char>",profile_level);
}
template <>
void AMP_MPI::call_allGather<char>(const char *x_in, int size_in, char *x_out, int *size_out, int *disp_out) const 
{
    PROFILE_START("allGatherv<char>",profile_level);
    MPI_Allgatherv( (void*) x_in, size_in, MPI_CHAR, (void*) x_out, size_out, disp_out, MPI_CHAR, communicator );
    PROFILE_STOP("allGatherv<char>",profile_level);
}
// unsigned int
template <>
void AMP_MPI::call_allGather<unsigned int>(const unsigned int x_in, unsigned int *x_out) const 
{
    PROFILE_START("allGather<unsigned int>",profile_level);
    MPI_Allgather( (void*) &x_in, 1, MPI_UNSIGNED, (void*) x_out, 1, MPI_UNSIGNED, communicator );
    PROFILE_STOP("allGather<unsigned int>",profile_level);
}
template <>
void AMP_MPI::call_allGather<unsigned int>(const unsigned int *x_in, int size_in, unsigned int *x_out, int *size_out, int *disp_out) const 
{
    PROFILE_START("allGatherv<unsigned int>",profile_level);
    MPI_Allgatherv( (void*) x_in, size_in, MPI_UNSIGNED, (void*) x_out, size_out, disp_out, MPI_UNSIGNED, communicator );
    PROFILE_STOP("allGatherv<unsigned int>",profile_level);
}
// int
template <>
void AMP_MPI::call_allGather<int>(const int x_in, int *x_out) const 
{
    PROFILE_START("allGather<int>",profile_level);
    MPI_Allgather( (void*) &x_in, 1, MPI_INT, (void*) x_out, 1, MPI_INT, communicator );
    PROFILE_STOP("allGather<int>",profile_level);
}
template <>
void AMP_MPI::call_allGather<int>(const int *x_in, int size_in, int *x_out, int *size_out, int *disp_out) const 
{
    PROFILE_START("allGatherv<int>",profile_level);
    MPI_Allgatherv( (void*) x_in, size_in, MPI_INT, (void*) x_out, size_out, disp_out, MPI_INT, communicator );
    PROFILE_STOP("allGatherv<int>",profile_level);
}
// unsigned long int
template <>
void AMP_MPI::call_allGather<unsigned long int>(const unsigned long int x_in, unsigned long int *x_out) const 
{
    PROFILE_START("allGather<unsigned long>",profile_level);
    MPI_Allgather( (void*) &x_in, 1, MPI_UNSIGNED_LONG, (void*) x_out, 1, MPI_UNSIGNED_LONG, communicator );
    PROFILE_STOP("allGather<unsigned long>",profile_level);
}
template <>
void AMP_MPI::call_allGather<unsigned long int>(const unsigned long int *x_in, int size_in, unsigned long int *x_out, int *size_out, int *disp_out) const 
{
    PROFILE_START("allGatherv<unsigned long>",profile_level);
    MPI_Allgatherv( (void*) x_in, size_in, MPI_UNSIGNED_LONG, (void*) x_out, size_out, disp_out, MPI_UNSIGNED_LONG, communicator );
    PROFILE_STOP("allGatherv<unsigned long>",profile_level);
}
// long int
template <>
void AMP_MPI::call_allGather<long int>(const long int x_in, long int *x_out) const 
{
    PROFILE_START("allGather<long int>",profile_level);
    MPI_Allgather( (void*) &x_in, 1, MPI_LONG, (void*) x_out, 1, MPI_LONG, communicator );
    PROFILE_STOP("allGather<long int>",profile_level);
}
template <>
void AMP_MPI::call_allGather<long int>(const long int *x_in, int size_in, long int *x_out, int *size_out, int *disp_out) const 
{
    PROFILE_START("allGatherv<long int>",profile_level);
    MPI_Allgatherv( (void*) x_in, size_in, MPI_LONG, (void*) x_out, size_out, disp_out, MPI_LONG, communicator );
    PROFILE_STOP("allGatherv<long int>",profile_level);
}
// float
template <>
void AMP_MPI::call_allGather<float>(const float x_in, float *x_out) const 
{
    PROFILE_START("allGather<float>",profile_level);
    MPI_Allgather( (void*) &x_in, 1, MPI_FLOAT, (void*) x_out, 1, MPI_FLOAT, communicator );
    PROFILE_STOP("allGather<float>",profile_level);
}
template <>
void AMP_MPI::call_allGather<float>(const float *x_in, int size_in, float *x_out, int *size_out, int *disp_out) const 
{
    PROFILE_START("allGatherv<float>",profile_level);
    MPI_Allgatherv( (void*) x_in, size_in, MPI_FLOAT, (void*) x_out, size_out, disp_out, MPI_FLOAT, communicator );
    PROFILE_STOP("allGatherv<float>",profile_level);
}
// double
template <>
void AMP_MPI::call_allGather<double>(const double x_in, double *x_out) const 
{
    PROFILE_START("allGather<double>",profile_level);
    MPI_Allgather( (void*) &x_in, 1, MPI_DOUBLE, (void*) x_out, 1, MPI_DOUBLE, communicator );
    PROFILE_STOP("allGather<double>",profile_level);
}
template <>
void AMP_MPI::call_allGather<double>(const double *x_in, int size_in, double *x_out, int *size_out, int *disp_out) const 
{
    PROFILE_START("allGatherv<double>",profile_level);
    MPI_Allgatherv( (void*) x_in, size_in, MPI_DOUBLE, (void*) x_out, size_out, disp_out, MPI_DOUBLE, communicator );
    PROFILE_STOP("allGatherv<double>",profile_level);
}
#else
// We need a concrete instantiation of call_allGather<char>(x_in,size_in,x_out,size_out)
template <>
void AMP_MPI::call_allGather<char>(const char *x_in, int size_in, char *x_out, int *size_out, int *disp_out) const 
{
    AMP_ERROR("Internal error in AMP_MPI (allGather) ");
}
#endif


/************************************************************************
*  allToAll                                                             *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
template <> void AMP_MPI::allToAll<unsigned char>(const int n, const unsigned char *send, unsigned char *recv ) const 
{
    PROFILE_START("allToAll<unsigned char>",profile_level);
    MPI_Alltoall( (void*) send, n, MPI_UNSIGNED_CHAR, (void*) recv, n, MPI_UNSIGNED_CHAR, communicator);
    PROFILE_STOP("allToAll<unsigned char>",profile_level);
}
template <> void AMP_MPI::allToAll<char>(const int n, const char *send, char *recv ) const 
{
    PROFILE_START("allToAll<char>",profile_level);
    MPI_Alltoall( (void*) send, n, MPI_CHAR, (void*) recv, n, MPI_CHAR, communicator);
    PROFILE_STOP("allToAll<char>",profile_level);
}
template <> void AMP_MPI::allToAll<unsigned int>(const int n, const unsigned int *send, unsigned int *recv ) const 
{
    PROFILE_START("allToAll<unsigned int>",profile_level);
    MPI_Alltoall( (void*) send, n, MPI_UNSIGNED, (void*) recv, n, MPI_UNSIGNED, communicator);
    PROFILE_STOP("allToAll<unsigned int>",profile_level);
}
template <> void AMP_MPI::allToAll<int>(const int n, const int *send, int *recv ) const 
{
    PROFILE_START("allToAll<int>",profile_level);
    MPI_Alltoall( (void*) send, n, MPI_INT, (void*) recv, n, MPI_INT, communicator);
    PROFILE_STOP("allToAll<int>",profile_level);
}
template <> void AMP_MPI::allToAll<unsigned long int>(const int n, const unsigned long int *send, unsigned long int *recv ) const 
{
    PROFILE_START("allToAll<unsigned long>",profile_level);
    MPI_Alltoall( (void*) send, n, MPI_UNSIGNED_LONG, (void*) recv, n, MPI_UNSIGNED_LONG, communicator);
    PROFILE_STOP("allToAll<unsigned long>",profile_level);
}
template <> void AMP_MPI::allToAll<long int>(const int n, const long int *send, long int *recv ) const 
{
    PROFILE_START("allToAll<long int>",profile_level);
    MPI_Alltoall( (void*) send, n, MPI_LONG, (void*) recv, n, MPI_LONG, communicator);
    PROFILE_STOP("allToAll<long int>",profile_level);
}
template <> void AMP_MPI::allToAll<float>(const int n, const float *send, float *recv ) const 
{
    PROFILE_START("allToAll<float>",profile_level);
    MPI_Alltoall( (void*) send, n, MPI_FLOAT, (void*) recv, n, MPI_FLOAT, communicator);
    PROFILE_STOP("allToAll<float>",profile_level);
}
template <> void AMP_MPI::allToAll<double>(const int n, const double *send, double *recv ) const 
{
    PROFILE_START("allToAll<double>",profile_level);
    MPI_Alltoall( (void*) send, n, MPI_DOUBLE, (void*) recv, n, MPI_DOUBLE, communicator);
    PROFILE_STOP("allToAll<double>",profile_level);
}
#endif


/************************************************************************
*  call_allToAll                                                        *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void AMP_MPI::call_allToAll<unsigned char>(const unsigned char *send_data, const int send_cnt[], 
        const int send_disp[], unsigned char *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    PROFILE_START("allToAllv<unsigned char>",profile_level);
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_UNSIGNED_CHAR, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_UNSIGNED_CHAR, communicator ); 
    PROFILE_STOP("allToAllv<unsigned char>",profile_level);
}
// char
template <>
void AMP_MPI::call_allToAll<char>(const char *send_data, const int send_cnt[], 
        const int send_disp[], char *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    PROFILE_START("allToAllv<char>",profile_level);
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_CHAR, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_CHAR, communicator );
    PROFILE_STOP("allToAllv<char>",profile_level); 
}
// unsigned int
template <>
void AMP_MPI::call_allToAll<unsigned int>(const unsigned int *send_data, const int send_cnt[], 
        const int send_disp[], unsigned int *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    PROFILE_START("allToAllv<unsigned int>",profile_level);
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_UNSIGNED, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_UNSIGNED, communicator ); 
    PROFILE_STOP("allToAllv<unsigned int>",profile_level);
}
// int
template <>
void AMP_MPI::call_allToAll<int>(const int *send_data, const int send_cnt[], 
        const int send_disp[], int *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    PROFILE_START("allToAllv<int>",profile_level);
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_INT, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_INT, communicator ); 
    PROFILE_STOP("allToAllv<int>",profile_level);
}
// unsigned long int
template <>
void AMP_MPI::call_allToAll<unsigned long int>(const unsigned long int *send_data, const int send_cnt[], 
        const int send_disp[], unsigned long int *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    PROFILE_START("allToAllv<unsigned long>",profile_level);
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_UNSIGNED_LONG, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_UNSIGNED_LONG, communicator ); 
    PROFILE_STOP("allToAllv<unsigned long>",profile_level);
}
// long int
template <>
void AMP_MPI::call_allToAll<long int>(const long int *send_data, const int send_cnt[], 
        const int send_disp[], long int *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    PROFILE_START("allToAllv<long int>",profile_level);
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_LONG, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_LONG, communicator ); 
    PROFILE_STOP("allToAllv<long int>",profile_level);
}
// float
template <>
void AMP_MPI::call_allToAll<float>(const float *send_data, const int send_cnt[], 
        const int send_disp[], float *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    PROFILE_START("allToAllv<float>",profile_level);
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_FLOAT, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_FLOAT, communicator ); 
    PROFILE_STOP("allToAllv<float>",profile_level);
}
// double
template <>
void AMP_MPI::call_allToAll<double>(const double *send_data, const int send_cnt[], 
        const int send_disp[], double *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    PROFILE_START("allToAllv<double>",profile_level);
    MPI_Alltoallv( (void*) send_data, (int*) send_cnt, (int*) send_disp, MPI_DOUBLE, 
        (void*) recv_data, (int*) recv_cnt, (int*) recv_disp, MPI_DOUBLE, communicator ); 
    PROFILE_STOP("allToAllv<double>",profile_level);
}
#else
// Default instatiation of unsigned char
template <>
void AMP_MPI::call_allToAll<char>(const char *send_data, const int send_cnt[], 
        const int send_disp[], char *recv_data, const int *recv_cnt, const int *recv_disp) const
{
    AMP_ERROR("Should not reach this point");
}
#endif


/************************************************************************
*  call_sumScan                                                         *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void AMP_MPI::call_sumScan<unsigned char>(const unsigned char *send, unsigned char *recv, int n) const 
{
    PROFILE_START("sumScan<unsigned char>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<unsigned char>",profile_level);
}
// char
template <>
void AMP_MPI::call_sumScan<char>(const char *send, char *recv, int n) const 
{
    PROFILE_START("sumScan<char>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<char>",profile_level);
}
// unsigned int
template <>
void AMP_MPI::call_sumScan<unsigned int>(const unsigned int *send, unsigned int *recv, int n) const 
{
    PROFILE_START("sumScan<unsigned int>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<unsigned int>",profile_level);
}
// int
template <>
void AMP_MPI::call_sumScan<int>(const int *send, int *recv, int n) const 
{
    PROFILE_START("sumScan<int>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_INT, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<int>",profile_level);
}
// long int
template <>
void AMP_MPI::call_sumScan<long int>(const long int *send, long int *recv, int n) const 
{
    PROFILE_START("sumScan<long int>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_LONG, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<long int>",profile_level);
}
// unsigned long int
template <>
void AMP_MPI::call_sumScan<unsigned long>(const unsigned long *send, unsigned long *recv, int n) const 
{
    PROFILE_START("sumScan<unsigned long>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<unsigned long>",profile_level);
}
// size_t
#ifdef USE_WINDOWS
template <>
void AMP_MPI::call_sumScan<size_t>(const size_t *send, size_t *recv, int n) const 
{
    PROFILE_START("sumScan<size_t>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_SIZE_T, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<size_t>",profile_level);
}
#endif
// float
template <>
void AMP_MPI::call_sumScan<float>(const float *send, float *recv, int n) const 
{
    PROFILE_START("sumScan<float>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_FLOAT, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<float>",profile_level);
}
// double
template <>
void AMP_MPI::call_sumScan<double>(const double *send, double *recv, int n) const 
{
    PROFILE_START("sumScan<double>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_DOUBLE, MPI_SUM, communicator);
    PROFILE_STOP("sumScan<double>",profile_level);
}
// std::complex<double>
template <>
void AMP_MPI::call_sumScan< std::complex<double> >(const std::complex<double> *x, std::complex<double> *y, int n) const 
{
    double *send = new double[2*n];
    double *recv = new double[2*n];
    for (int i=0; i<n; i++) {
        send[2*i+0] = real(x[i]);
        send[2*i+1] = imag(x[i]);
    }
    MPI_Scan( (void*) send, (void*) recv, 2*n, MPI_DOUBLE, MPI_SUM, communicator);
    for (int i=0; i<n; i++)
        y[i] = std::complex<double>(recv[2*i+0],recv[2*i+1]);
    delete [] send;
    delete [] recv;
}
#endif


/************************************************************************
*  call_minScan                                                         *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void AMP_MPI::call_minScan<unsigned char>(const unsigned char *send, unsigned char *recv, int n) const 
{
    PROFILE_START("minScan<unsigned char>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_MIN, communicator);
    PROFILE_STOP("minScan<unsigned char>",profile_level);
}
// char
template <>
void AMP_MPI::call_minScan<char>(const char *send, char *recv, int n) const 
{
    PROFILE_START("minScan<char>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_MIN, communicator);
    PROFILE_STOP("minScan<char>",profile_level);
}
// unsigned int
template <>
void AMP_MPI::call_minScan<unsigned int>(const unsigned int *send, unsigned int *recv, int n) const 
{
    PROFILE_START("minScan<unsigned int>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_MIN, communicator);
    PROFILE_STOP("minScan<unsigned int>",profile_level);
}
// int
template <>
void AMP_MPI::call_minScan<int>(const int *send, int *recv, int n) const 
{
    PROFILE_START("minScan<int>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_INT, MPI_MIN, communicator);
    PROFILE_STOP("minScan<int>",profile_level);
}
// unsigned long int
template <>
void AMP_MPI::call_minScan<unsigned long int>(const unsigned long int *send, unsigned long int *recv, int n) const 
{
    PROFILE_START("minScan<unsigned long>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_MIN, communicator);
    PROFILE_STOP("minScan<unsigned long>",profile_level);
}
// long int
template <>
void AMP_MPI::call_minScan<long int>(const long int *send, long int *recv, int n) const 
{
    PROFILE_START("minScan<long int>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_LONG, MPI_MIN, communicator);
    PROFILE_STOP("minScan<long int>",profile_level);
}
// size_t
#ifdef USE_WINDOWS
template <>
void AMP_MPI::call_minScan<size_t>(const size_t *send, size_t *recv, int n) const 
{
    PROFILE_START("minScan<size_t>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_SIZE_T, MPI_MIN, communicator);
    PROFILE_STOP("minScan<size_t>",profile_level);
}
#endif
// float
template <>
void AMP_MPI::call_minScan<float>(const float *send, float *recv, int n) const 
{
    PROFILE_START("minScan<float>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_FLOAT, MPI_MIN, communicator);
    PROFILE_STOP("minScan<float>",profile_level);
}
// double
template <>
void AMP_MPI::call_minScan<double>(const double *send, double *recv, int n) const 
{
    PROFILE_START("minScan<double>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_DOUBLE, MPI_MIN, communicator);
    PROFILE_STOP("minScan<double>",profile_level);
}
#endif


/************************************************************************
*  call_maxScan                                                         *
*  Note: these specializations are only called when using MPI.          *
************************************************************************/
#ifdef USE_MPI
// unsigned char
template <>
void AMP_MPI::call_maxScan<unsigned char>(const unsigned char *send, unsigned char *recv, int n) const 
{
    PROFILE_START("maxScan<unsigned char>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_CHAR, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<unsigned char>",profile_level);
}
// char
template <>
void AMP_MPI::call_maxScan<char>(const char *send, char *recv, int n) const 
{
    PROFILE_START("maxScan<char>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_SIGNED_CHAR, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<char>",profile_level);
}
// unsigned int
template <>
void AMP_MPI::call_maxScan<unsigned int>(const unsigned int *send, unsigned int *recv, int n) const 
{
    PROFILE_START("maxScan<unsigned int>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<unsigned int>",profile_level);
}
// int
template <>
void AMP_MPI::call_maxScan<int>(const int *send, int *recv, int n) const 
{
    PROFILE_START("maxScan<int>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_INT, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<int>",profile_level);
}
// long int
template <>
void AMP_MPI::call_maxScan<long int>(const long int *send, long int *recv, int n) const 
{
    PROFILE_START("maxScan<long int>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_LONG, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<long int>",profile_level);
}
// unsigned long int
template <>
void AMP_MPI::call_maxScan<unsigned long int>(const unsigned long int *send, unsigned long int *recv, int n) const 
{
    PROFILE_START("maxScan<unsigned long>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_UNSIGNED_LONG, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<unsigned long>",profile_level);
}
// size_t
#ifdef USE_WINDOWS
template <>
void AMP_MPI::call_maxScan<size_t>(const size_t *send, size_t *recv, int n) const 
{
    PROFILE_START("maxScan<size_t>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_SIZE_T, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<size_t>",profile_level);
}
#endif
// float
template <>
void AMP_MPI::call_maxScan<float>(const float *send, float *recv, int n) const 
{
    PROFILE_START("maxScan<float>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_INT, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<float>",profile_level);
}
// double
template <>
void AMP_MPI::call_maxScan<double>(const double *send, double *recv, int n) const 
{
    PROFILE_START("maxScan<double>",profile_level);
    MPI_Scan( (void*) send, (void*) recv, n, MPI_DOUBLE, MPI_MAX, communicator);
    PROFILE_STOP("maxScan<double>",profile_level);
}
#endif


/************************************************************************
*  Wait functions                                                       *
************************************************************************/
#ifdef USE_MPI
void AMP_MPI::wait( MPI_Request request) 
{
    PROFILE_START("wait",profile_level);
    MPI_Status  status;
    int flag = 0;
    int err = MPI_Test( &request, &flag, &status );
    AMP_ASSERT(err==MPI_SUCCESS);   // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Test( &request, &flag, &status );
    }
    PROFILE_STOP("wait",profile_level);
}
int AMP_MPI::waitAny( int count, MPI_Request *request) 
{
    if ( count==0 ) 
        return -1;
    PROFILE_START("waitAny",profile_level);
    int index = -1;
    int flag = 0;
    MPI_Status *status = new MPI_Status[count];
    int err = MPI_Testany( count, request, &index, &flag, status );
    AMP_ASSERT(err==MPI_SUCCESS);   // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testany( count, request, &index, &flag, status );
    }
    AMP_ASSERT(index>=0);   // Check that the index is valid
    delete [] status;
    PROFILE_STOP("waitAny",profile_level);
    return index;
}
void AMP_MPI::waitAll( int count, MPI_Request *request) 
{
    if ( count==0 ) 
        return;
    PROFILE_START("waitAll",profile_level);
    int flag = 0;
    MPI_Status *status = new MPI_Status[count];
    int err = MPI_Testall( count, request, &flag, status );
    AMP_ASSERT(err==MPI_SUCCESS);   // Check that the first call is valid
    while ( !flag ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testall( count, request, &flag, status );
    }
    PROFILE_STOP("waitAll",profile_level);
    delete [] status;
}
std::vector<int> AMP_MPI::waitSome( int count, MPI_Request *request )
{
    if ( count==0 ) 
        return std::vector<int>();
    PROFILE_START("waitSome",profile_level);
    std::vector<int> indicies(count,-1);
    MPI_Status *status = new MPI_Status[count];
    int outcount=0;
    int err = MPI_Testsome( count, request, &outcount, &indicies[0], status );
    AMP_ASSERT(err==MPI_SUCCESS);           // Check that the first call is valid
    AMP_ASSERT(outcount!=MPI_UNDEFINED);    // Check that the first call is valid
    while ( outcount==0 ) {
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
        // Check if the request has finished
        MPI_Testsome( count, request, &outcount, &indicies[0], status );
    }
    indicies.resize(outcount);
    delete [] status;
    PROFILE_STOP("waitSome",profile_level);
    return indicies;
}
#else
void AMP_MPI::wait( MPI_Request request) 
{
    PROFILE_START("wait",profile_level);
    while ( 1 ) {
        // Check if the request is in our list
        if ( global_isendrecv_list.find(request)==global_isendrecv_list.end() )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP("wait",profile_level);
}
int AMP_MPI::waitAny( int count, MPI_Request *request) 
{
    if ( count==0 ) 
        return -1;
    PROFILE_START("waitAny",profile_level);
    int index = 0;
    while ( 1 ) {
        // Check if the request is in our list
        bool found_any = false;
        for (int i=0; i<count; i++) {
            if ( global_isendrecv_list.find(request[i])==global_isendrecv_list.end() ) {
                found_any = true;
                index = i;
            }
        }
        if ( found_any )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP("waitAny",profile_level);
    return index;
}
void AMP_MPI::waitAll( int count, MPI_Request *request) 
{
    if ( count==0 ) 
        return;
    PROFILE_START("waitAll",profile_level);
    while ( 1 ) {
        // Check if the request is in our list
        bool found_all = true;
        for (int i=0; i<count; i++) {
            if ( global_isendrecv_list.find(request[i])!=global_isendrecv_list.end() )
                found_all = false;
        }
        if ( found_all )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP("waitAll",profile_level);
}
std::vector<int> AMP_MPI::waitSome( int count, MPI_Request *request )
{
    if ( count==0 ) 
        return std::vector<int>();
    PROFILE_START("waitSome",profile_level);
    std::vector<int> indicies();
    while ( 1 ) {
        // Check if the request is in our list
        for (int i=0; i<count; i++) {
            if ( global_isendrecv_list.find(request[i])==global_isendrecv_list.end() ) {
                found_any = true;
                indicies.push_back(i);
            }
        }
        if ( !indicies.empty() )
            break;
        // Put the current thread to sleep to allow other threads to run
        sched_yield();
    }
    PROFILE_STOP("waitSome",profile_level);
    return indicies;
}
#endif


/************************************************************************
*  Probe functions                                                      *
************************************************************************/
#ifdef USE_MPI
int AMP_MPI::Iprobe( int source, int tag) const 
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Status status;
    int flag = 0;
    MPI_Iprobe(source,tag,communicator,&flag,&status);
    if ( flag==0 )
        return -1;
    int count;
    MPI_Get_count(&status,MPI_BYTE,&count);
    AMP_ASSERT(count>=0);
    return count;
}
int AMP_MPI::probe( int source, int tag) const 
{
    AMP_INSIST(tag<=d_maxTag,"Maximum tag value exceeded");
    AMP_INSIST(tag>=0,"tag must be >= 0");
    MPI_Status status;
    MPI_Probe(source,tag,communicator,&status);
    int count;
    MPI_Get_count(&status,MPI_BYTE,&count);
    AMP_ASSERT(count>=0);
    return count;
}
#else
int AMP_MPI::Iprobe( int source, int tag) const 
{
    AMP_ERROR("Not implimented for serial codes (Iprobe)");
    return 0;
}
int AMP_MPI::probe( int source, int tag) const 
{
    AMP_ERROR("Not implimented for serial codes (probe)");
    return 0;
}
#endif



/************************************************************************
*  Timer functions                                                      *
************************************************************************/
#ifdef USE_MPI
    double AMP_MPI::time() 
    { 
        return MPI_Wtime();
    }
    double AMP_MPI::tick() 
    { 
        return MPI_Wtick();
    }
#else
    #if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        double AMP_MPI::time() 
        { 
            LARGE_INTEGER end, f;
            QueryPerformanceFrequency(&f);
            QueryPerformanceCounter(&end);       
            double time = ((double)end.QuadPart)/((double)f.QuadPart);
            return time;
        }
        double AMP_MPI::tick() 
        { 
            LARGE_INTEGER f;
            QueryPerformanceFrequency(&f);
            double resolution = ((double)1.0)/((double)f.QuadPart);
            return resolution;
        }
    #else
        double AMP_MPI::time() 
        { 
            timeval current_time;
            gettimeofday(&current_time,NULL);
            double time = ((double)current_time.tv_sec)+1e-6*((double)current_time.tv_usec);
            return time;
        }
        double AMP_MPI::tick() 
        { 
            timeval start, end;
            gettimeofday(&start,NULL);
            gettimeofday(&end,NULL);
            while ( end.tv_sec==start.tv_sec &&  end.tv_usec==start.tv_usec )
                gettimeofday(&end,NULL);
            double resolution = ((double)(end.tv_sec-start.tv_sec))+1e-6*((double)(end.tv_usec-start.tv_usec));
            return resolution;
        }
    #endif
#endif


} // namespace

