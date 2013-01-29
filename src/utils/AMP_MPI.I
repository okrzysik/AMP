// This file contains the default instantiations for templated operations
// Note: Intel compilers need definitions before all default instantions to compile correctly


/************************************************************************
*  sumReduce                                                            *
************************************************************************/
template <class type>
inline type AMP_MPI::sumReduce(const type value) const {
    if ( comm_size > 1 ) {
        type tmp = value;
        call_sumReduce(&tmp,1);
        return tmp;
    } else {
        return value;
    }
}
template <class type>
inline void AMP_MPI::sumReduce(type *x, const int n) const {
    if ( comm_size > 1 ) 
        call_sumReduce(x,n);
}
template <class type>
inline void AMP_MPI::sumReduce(const type *x, type *y, const int n) const {
    if ( comm_size > 1 ) {
        call_sumReduce(x,y,n);
    } else {
        for (int i=0; i<n; i++)
            y[i] = x[i];
    }
}
// Define specializations of call_sumReduce(type*, const int)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_sumReduce<unsigned char>(unsigned char*, const int) const;
template <> void AMP_MPI::call_sumReduce<char>(char*, const int) const;
template <> void AMP_MPI::call_sumReduce<unsigned int>(unsigned int*, const int) const;
template <> void AMP_MPI::call_sumReduce<int>(int*, const int) const;
template <> void AMP_MPI::call_sumReduce<unsigned long int>(unsigned long int*, const int) const;
template <> void AMP_MPI::call_sumReduce<long int>(long int*, const int) const;
template <> void AMP_MPI::call_sumReduce<float>(float*, const int) const;
template <> void AMP_MPI::call_sumReduce<double>(double*, const int) const;
template <> void AMP_MPI::call_sumReduce< std::complex<double> >(std::complex<double>*, const int) const;
#endif
// Default instantiations of call_sumReduce(type*, const int)
template <class type>
void AMP_MPI::call_sumReduce(type *x, const int n) const {
    char message[200];
    sprintf(message,"Default instantion of sumReduce in parallel is not supported (%s)",typeid(type).name());
    AMP_ERROR(message);
}
// Define specializations of call_sumReduce(const type*, type*, const int)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_sumReduce<unsigned char>(const unsigned char*, unsigned char*, const int) const;
template <> void AMP_MPI::call_sumReduce<char>(const char*,  char*, const int) const;
template <> void AMP_MPI::call_sumReduce<unsigned int>(const unsigned int*, unsigned int*, const int) const;
template <> void AMP_MPI::call_sumReduce<int>(const int*, int*, const int) const;
template <> void AMP_MPI::call_sumReduce<unsigned long int>(const unsigned long int*, unsigned long int*, const int) const;
template <> void AMP_MPI::call_sumReduce<long int>(const long int*, long int*, const int) const;
template <> void AMP_MPI::call_sumReduce<float>(const float*, float*, const int) const;
template <> void AMP_MPI::call_sumReduce<double>(const double*, double*, const int) const;
template <> void AMP_MPI::call_sumReduce< std::complex<double> >(const std::complex<double>*, std::complex<double>*, const int) const;
#endif
// Default instantiations of call_sumReduce(const type*, type*, const int)
template <class type>
void AMP_MPI::call_sumReduce(const type *x, type *y, const int n) const {
    char message[200];
    sprintf(message,"Default instantion of sumReduce in parallel is not supported (%s)",typeid(type).name());
    AMP_ERROR(message);
}


/************************************************************************
*  minReduce                                                            *
************************************************************************/
template <class type>
inline type AMP_MPI::minReduce(const type value) const {
    if ( comm_size > 1 ) {
        type tmp = value;
        call_minReduce(&tmp,1,NULL);
        return tmp;
    } else {
        return value;
    }
}
template <class type>
inline void AMP_MPI::minReduce(type *x, const int n, int *rank_of_min) const {
    if ( comm_size > 1 ) {
        call_minReduce(x,n,rank_of_min);
    } else {
        if ( rank_of_min != NULL ) {
            for (int i=0; i<n; i++)
                rank_of_min[i] = 0;
        }
    }
}
template <class type>
inline void AMP_MPI::minReduce(const type *x, type *y, const int n, int *rank_of_min) const {
    if ( comm_size > 1 ) {
        call_minReduce(x,y,n,rank_of_min);
    } else {
        for (int i=0; i<n; i++) {
            y[i] = x[i];
            if ( rank_of_min != NULL )
                rank_of_min[i] = 0;
        }
    }
}
// Define specializations of call_minReduce(type*, const int, int*)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_minReduce<unsigned char>(unsigned char*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<char>(char*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<unsigned int>(unsigned int*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<int>(int*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<unsigned long int>(unsigned long int*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<long int>(long int*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<float>(float*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<double>(double*, const int, int*) const;
#endif
// Default instantiations of call_minReduce(type*, const int, int*)
template <class type>
void AMP_MPI::call_minReduce(type *x, const int n, int *rank_of_min) const {
    char message[200];
    sprintf(message,"Default instantion of minReduce in parallel is not supported (%s)",typeid(type).name());
    AMP_ERROR(message);
}
// Define specializations of call_minReduce(const type*, type*, const int, int*)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_minReduce<unsigned char>(const unsigned char*, unsigned char*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<char>(const char*, char*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<unsigned int>(const unsigned int*, unsigned int*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<int>(const int*, int*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<unsigned long int>(const unsigned long int*, unsigned long int*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<long int>(const long int*, long int*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<float>(const float*, float*, const int, int*) const;
template <> void AMP_MPI::call_minReduce<double>(const double*, double*, const int, int*) const;
#endif
// Default instantiations of call_minReduce(const type*, type*, const int, int*)
template <class type>
void AMP_MPI::call_minReduce(const type *x, type *y, const int n, int *rank_of_min) const {
    char message[200];
    sprintf(message,"Default instantion of minReduce in parallel is not supported (%s)",typeid(type).name());
    AMP_ERROR(message);
}


/************************************************************************
*  maxReduce                                                         *
************************************************************************/
template <class type>
inline type AMP_MPI::maxReduce(const type value) const {
    if ( comm_size > 1 ) {
        type tmp = value;
        call_maxReduce(&tmp,1,NULL);
        return tmp;
    } else {
        return value;
    }
}
template <class type>
inline void AMP_MPI::maxReduce(type *x, const int n, int *rank_of_max) const {
    if ( comm_size > 1 ) {
        call_maxReduce(x,n,rank_of_max);
    } else {
        if ( rank_of_max != NULL ) {
            for (int i=0; i<n; i++)
                rank_of_max[i] = 0;
        }
    }
}
template <class type>
inline void AMP_MPI::maxReduce(const type *x, type *y, const int n, int *rank_of_max) const {
    if ( comm_size > 1 ) {
        call_maxReduce(x,y,n,rank_of_max);
    } else {
        for (int i=0; i<n; i++) {
            y[i] = x[i];
            if ( rank_of_max != NULL )
                rank_of_max[i] = 0;
        }
    }
}
// Define specializations of call_maxReduce(type*, const int, int*)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_maxReduce<unsigned char>(unsigned char*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<char>(char*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<unsigned int>(unsigned int*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<int>(int*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<unsigned long int>(unsigned long int*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<long int>(long int*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<float>(float*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<double>(double*, const int, int*) const;
#endif
// Default instantiations of call_maxReduce(type*, const int, int*)
template <class type>
void AMP_MPI::call_maxReduce(type *x, const int n, int *rank_of_max) const {
    char message[200];
    sprintf(message,"Default instantion of maxReduce in parallel is not supported (%s)",typeid(type).name());
    AMP_ERROR(message);
}
// Define specializations of call_maxReduce(const type*, type*, const int, int*)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_maxReduce<unsigned char>(const unsigned char*, unsigned char*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<char>(const char*, char*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<unsigned int>(const unsigned int*, unsigned int*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<int>(const int*, int*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<unsigned long int>(const unsigned long int*, unsigned long int*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<long int>(const long int*, long int*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<float>(const float*, float*, const int, int*) const;
template <> void AMP_MPI::call_maxReduce<double>(const double*, double*, const int, int*) const;
#endif
// Default instantiations of call_maxReduce(const type*, type*, const int, int*)
template <class type>
void AMP_MPI::call_maxReduce(const type *x, type *y, const int n, int *rank_of_max) const {
    char message[200];
    sprintf(message,"Default instantion of maxReduce in parallel is not supported (%s)",typeid(type).name());
    AMP_ERROR(message);
}


/************************************************************************
*  bcast                                                                *
************************************************************************/
template <class type>
inline type AMP_MPI::bcast(const type value, const int root) const {
    if ( root >= comm_size )
        AMP_ERROR("root cannot be >= size in bcast");
    if ( comm_size > 1 ) {
        type tmp = value;
        call_bcast(&tmp,1,root);
        return tmp;
    } else {
        return value;
    }
}
template <class type>
inline void AMP_MPI::bcast(type *x, const int n, const int root) const {
    if ( root >= comm_size )
        AMP_ERROR("root cannot be >= size in bcast");
    if ( comm_size > 1 )
        call_bcast(x,n,root);
}
// Define specializations of bcast(type*, const int, const int)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_bcast<unsigned char>(unsigned char*, const int, const int) const;
template <> void AMP_MPI::call_bcast<char>(char*, const int, const int) const;
template <> void AMP_MPI::call_bcast<unsigned int>(unsigned int*, const int, const int) const;
template <> void AMP_MPI::call_bcast<int>(int*, const int, const int) const;
template <> void AMP_MPI::call_bcast<float>(float*, const int, const int) const;
template <> void AMP_MPI::call_bcast<double>(double*, const int, const int) const;
#else
template <> void AMP_MPI::call_bcast<char>(char*, const int, const int) const;
#endif
// Default instantiations of bcast(type*, const int, const int)
template <class type>
void AMP_MPI::call_bcast(type *x, const int n, const int root) const {
    call_bcast<char>((char*)x,(int)n*sizeof(type),root);
}


/************************************************************************
*  send                                                                 *
************************************************************************/
// Define specializations of send(const type*, const int, const int, int)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::send<char>(const char*, const int, const int, int) const;
template <> void AMP_MPI::send<int>(const int*, int, const int, int) const;
template <> void AMP_MPI::send<float>(const float*, const int, const int, int) const;
template <> void AMP_MPI::send<double>(const double*, const int, const int, int) const;
#else
template <> void AMP_MPI::send<char>(const char*, const int, const int, int) const;
#endif
// Default instantiations of send(const type*, const int, const int, int)
template <class type>
inline void AMP_MPI::send(const type *buf, const int length, const int recv_proc_number, int tag) const
{
    send<char>((const char*)buf,length*sizeof(type),recv_proc_number,tag);
}


/************************************************************************
*  Isend                                                                *
************************************************************************/
// Define specializations of Isend(const type*, const int, const int, const int)
#ifdef USE_EXT_MPI
template <> MPI_Request AMP_MPI::Isend<char>(const char*, const int, const int, const int) const;
template <> MPI_Request AMP_MPI::Isend<int>(const int*, int, const int, const int) const;
template <> MPI_Request AMP_MPI::Isend<float>(const float*, const int, const int, const int) const;
template <> MPI_Request AMP_MPI::Isend<double>(const double*, const int, const int, const int) const;
#else
template <> MPI_Request AMP_MPI::Isend<char>(const char*, const int, const int, const int) const;
#endif
// Default instantiations of Isend(const type*, const int, const int, const int)
template <class type>
inline MPI_Request AMP_MPI::Isend(const type *buf, const int length, const int recv_proc_number, const int tag) const
{
    return Isend<char>((const char*)buf,length*sizeof(type),recv_proc_number,tag);
}


/************************************************************************
*  recv                                                                 *
************************************************************************/
// Define specializations of recv(type*, int&, const int, const bool, int)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::recv<char>(char*, int&, const int, const bool, int) const;
template <> void AMP_MPI::recv<int>(int*, int&, const int, const bool, int) const;
template <> void AMP_MPI::recv<float>(float*, int&, const int, const bool, int) const;
template <> void AMP_MPI::recv<double>(double*, int&, const int, const bool, int) const;
#else
template <> void AMP_MPI::recv<char>(char*, int&, const int, const bool, int) const;
#endif
// Default instantiations of recv(type*, int&, const int, const bool, int)
template <class type>
inline void AMP_MPI::recv(type *buf, int &length, const int send_proc_number, const bool get_length, int tag) const
{
    int size = length*sizeof(type);
    recv<char>((char*)buf,size,send_proc_number,get_length,tag);
    if ( get_length ) {
        AMP_ASSERT(size%sizeof(type)==0);
        length = size/sizeof(type);
    }
}


/************************************************************************
*  Irecv                                                                *
************************************************************************/
// Define specializations of recv(type*, int&, const int, const bool, int)
#ifdef USE_EXT_MPI
template <> MPI_Request AMP_MPI::Irecv<char>(char*, const int, const int, const int) const;
template <> MPI_Request AMP_MPI::Irecv<int>(int*, const int, const int, const int) const;
template <> MPI_Request AMP_MPI::Irecv<float>(float*, const int, const int, const int) const;
template <> MPI_Request AMP_MPI::Irecv<double>(double*, const int, const int, const int) const;
#else
template <> MPI_Request AMP_MPI::Irecv<char>(char*, const int, const int, const int) const;
#endif
// Default instantiations of recv(type*, int&, const int, const bool, int)
template <class type>
inline MPI_Request AMP_MPI::Irecv(type *buf, const int length, const int send_proc, const int tag) const
{
    return Irecv<char>((char*)buf,length*sizeof(type),send_proc,tag);
}


/************************************************************************
*  allGather                                                            *
************************************************************************/
template <class type>
inline void AMP_MPI::allGather(const type x_in, type *x_out) const
{
    if ( comm_size > 1 ) {
        // We can use the vector form of allGather with a char array to ge the data we want
        call_allGather( x_in, x_out );
    } else {
        // Single processor case
        x_out[0] = x_in;
    }
}
template <class type>
int AMP_MPI::allGather(const type *send_data, const int send_cnt, type *recv_data, 
    int *recv_cnt, int* recv_disp, bool known_recv) const
{
    // Check the inputs
    if ( known_recv && ( recv_cnt==NULL || recv_disp==NULL ) )
        AMP_ERROR("Error calling allGather");
    // Check if we are dealing with a single processor
    if ( comm_size==1 ) {
        if ( !known_recv ) {
            // We do not know the recieved sizes
            for (int i=0; i<send_cnt; i++)
                recv_data[i] = send_data[i];
            if ( recv_cnt!=NULL )
                recv_cnt[0] = send_cnt;
            if ( recv_disp!=NULL )
                recv_disp[0] = 0;
        } else {
            // We know the recieved sizes
            for (int i=0; i<send_cnt; i++)
                recv_data[i+recv_disp[0]] = send_data[i];
        }
        return send_cnt;
    }
    // Get the sizes of the recieved data (if necessary)
    int *recv_cnt2 = recv_cnt;
    int *recv_disp2 = recv_disp;
    if ( !known_recv ) {
        if ( recv_cnt==NULL )
            recv_cnt2 = new int[comm_size];
        if ( recv_disp==NULL )
            recv_disp2 = new int[comm_size];
        call_allGather( send_cnt, recv_cnt2 );
        recv_disp2[0] = 0;
        for (int i=1; i<comm_size; i++)
            recv_disp2[i] = recv_disp2[i-1]+recv_cnt2[i-1];
    }
    int N_recv = 0;
    for (int i=0; i<comm_size; i++)
        N_recv += recv_cnt2[i];
    // Send/recv the data
    call_allGather( send_data, send_cnt, recv_data, recv_cnt2, recv_disp2 );
    // Delete any temporary memory
    if ( recv_cnt==NULL )
        delete [] recv_cnt2;
    if ( recv_disp==NULL )
        delete [] recv_disp2;
    return N_recv;
}
// Define specializations of call_allGather(const type, type*) 
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_allGather<unsigned char>(const unsigned char, unsigned char*) const;
template <> void AMP_MPI::call_allGather<char>(const char, char*) const;
template <> void AMP_MPI::call_allGather<unsigned int>(const unsigned int, unsigned int*) const;
template <> void AMP_MPI::call_allGather<int>(const int, int*) const;
template <> void AMP_MPI::call_allGather<unsigned long int>(const unsigned long int, unsigned long int*) const;
template <> void AMP_MPI::call_allGather<long int>(const long int, long int*) const;
template <> void AMP_MPI::call_allGather<float>(const float, float*) const;
template <> void AMP_MPI::call_allGather<double>(const double, double*) const;
#endif
// Default instantiations of call_allGather(const type, type*) 
template <class type>
void AMP_MPI::call_allGather(const type x_in, type *x_out)  const
{
    allGather<char>( (const char*) &x_in, (int) sizeof(type), (char*) x_out );
}
// Define specializations of call_allGather(const type*, int, type*, int*, int*)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_allGather<unsigned char>(const unsigned char*, int, unsigned char*, int*, int*) const;
template <> void AMP_MPI::call_allGather<char>(const char*, int, char*, int*, int*) const;
template <> void AMP_MPI::call_allGather<unsigned int>(const unsigned int*, int, unsigned int*, int*, int*) const;
template <> void AMP_MPI::call_allGather<int>(const int*, int, int*, int*, int*) const;
template <> void AMP_MPI::call_allGather<unsigned long int>(const unsigned long int*, int, unsigned long int*, int*, int*) const;
template <> void AMP_MPI::call_allGather<long int>(const long int*, int, long int*, int*, int*) const;
template <> void AMP_MPI::call_allGather<float>(const float*, int, float*, int*, int*) const;
template <> void AMP_MPI::call_allGather<double>(const double*, int, double*, int*, int*) const;
#else
template <> void AMP_MPI::call_allGather<char>(const char*, int, char*, int*, int*) const;
#endif
// Default instantiations of int call_allGather(const type*, int, type*, int*)
template <class type>
void AMP_MPI::call_allGather(const type *x_in, int size_in, type *x_out, int *size_out, int *disp_out) const
{
    int *size2 = new int[comm_size];
    int *disp2 = new int[comm_size];
    for (int i=0; i<comm_size; i++) {
        size2[i] = size_out[i]*sizeof(type);
        disp2[i] = disp_out[i]*sizeof(type);
    }
    call_allGather<char>( (const char*) x_in, (int)size_in*sizeof(type), (char*) x_out, size2, disp2 );
    delete [] size2;
    delete [] disp2;
}


/************************************************************************
*  setGather                                                            *
************************************************************************/
template <class type>
inline void AMP_MPI::setGather( std::set<type> &set ) const
{
    std::vector<type> send_buf(set.begin(),set.end());
    std::vector<int> recv_cnt(this->comm_size,0);
    this->allGather<int>( (int)send_buf.size(), &recv_cnt[0] );
    std::vector<int> recv_disp(this->comm_size,0);
    for (int i=1; i<this->comm_size; i++)
        recv_disp[i] = recv_disp[i-1]+recv_cnt[i-1];
    size_t N_recv_tot = 0;
    for (int i=0; i<this->comm_size; i++)
        N_recv_tot += recv_cnt[i];
    if ( N_recv_tot==0 )
        return;
    std::vector<type> recv_buf(N_recv_tot);
    type *send_data = NULL;
    if ( send_buf.size()>0 ) { send_data = &send_buf[0]; }
    type *recv_data = &recv_buf[0];
    this->allGather<type>( send_data, (int)send_buf.size(), recv_data, &recv_cnt[0], &recv_disp[0], true);
    for (size_t i=0; i<recv_buf.size(); i++)
        set.insert(recv_buf[i]);
}


/************************************************************************
*  mapGather                                                            *
************************************************************************/
template <class KEY, class DATA>
inline void AMP_MPI::mapGather( std::map<KEY,DATA> &map ) const
{
    std::vector<KEY> send_id;
    std::vector<DATA> send_data;
    send_id.reserve(map.size());
    send_data.reserve(map.size());
    std::map<KEY,DATA> map2;
    typename std::map<KEY,DATA>::iterator it;
    for (it=map.begin(); it!=map.end(); ++it) {
        send_id.push_back(it->first);
        send_data.push_back(it->second);
    }
    int send_size = (int)send_id.size();
    std::vector<int> recv_cnt(this->comm_size,0);
    this->allGather<int>( send_size, &recv_cnt[0] );
    std::vector<int> recv_disp(this->comm_size,0);
    for (int i=1; i<this->comm_size; i++)
        recv_disp[i] = recv_disp[i-1]+recv_cnt[i-1];
    size_t N_recv_tot = 0;
    for (int i=0; i<this->comm_size; i++)
        N_recv_tot += recv_cnt[i];
    if ( N_recv_tot==0 )
        return;
    std::vector<KEY> recv_id(N_recv_tot);
    std::vector<DATA> recv_data(N_recv_tot);
    KEY *send_data1 = NULL;
    DATA *send_data2 = NULL;
    if ( send_id.size()>0 ) { 
        send_data1 = &send_id[0]; 
        send_data2 = &send_data[0]; 
    }
    this->allGather<KEY>(   send_data1, send_size, &recv_id[0],   &recv_cnt[0], &recv_disp[0], true);
    this->allGather<DATA>( send_data2, send_size, &recv_data[0], &recv_cnt[0], &recv_disp[0], true);
    map = std::map<KEY,DATA>();
    for (size_t i=0; i<N_recv_tot; i++)
        map.insert( std::pair<KEY,DATA>(recv_id[i],recv_data[i]) );
}


/************************************************************************
*  sumScan                                                              *
************************************************************************/
template <class type>
inline void AMP_MPI::sumScan(const type *x, type *y, const int n) const {
    if ( comm_size > 1 ) {
        call_sumScan(x,y,n);
    } else {
        for (int i=0; i<n; i++)
            y[i] = x[i];
    }
}
// Define specializations of call_sumScan(const type*, type*, int)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_sumScan<unsigned char>(const unsigned char*, unsigned char*, int) const;
template <> void AMP_MPI::call_sumScan<char>(const char*, char*, int) const;
template <> void AMP_MPI::call_sumScan<unsigned int>(const unsigned int*, unsigned int*, int) const;
template <> void AMP_MPI::call_sumScan<int>(const int*, int*, int) const;
template <> void AMP_MPI::call_sumScan<unsigned long int>(const unsigned long int*, unsigned long int*, int) const;
template <> void AMP_MPI::call_sumScan<long int>(const long int*, long int*, int) const;
template <> void AMP_MPI::call_sumScan<float>(const float*, float*, int) const;
template <> void AMP_MPI::call_sumScan<double>(const double*, double*, int) const;
template <> void AMP_MPI::call_sumScan< std::complex<double> >(const std::complex<double>*, std::complex<double>*, int) const;
#endif
// Default instantiations of call_sumScan(const type*, type*, int)
template <class type>
void AMP_MPI::call_sumScan(const type *x, type *y, int n) const {
    char message[200];
    sprintf(message,"Default instantion of sumScan in parallel is not supported (%s)",typeid(type).name());
    AMP_ERROR(message);
}


/************************************************************************
*  minScan                                                              *
************************************************************************/
template <class type>
inline void AMP_MPI::minScan(const type *x, type *y, const int n) const {
    if ( comm_size > 1 ) {
        call_minScan(x,y,n);
    } else {
        for (int i=0; i<n; i++)
            y[i] = x[i];
    }
}
// Define specializations of call_minScan(const type*, type*, int)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_minScan<unsigned char>(const unsigned char*, unsigned char*, int) const;
template <> void AMP_MPI::call_minScan<char>(const char*, char*, int) const;
template <> void AMP_MPI::call_minScan<unsigned int>(const unsigned int*, unsigned int*, int) const;
template <> void AMP_MPI::call_minScan<int>(const int*, int*, int) const;
template <> void AMP_MPI::call_minScan<unsigned long int>(const unsigned long int*, unsigned long int*, int) const;
template <> void AMP_MPI::call_minScan<long int>(const long int*, long int*, int) const;
template <> void AMP_MPI::call_minScan<float>(const float*, float*, int) const;
template <> void AMP_MPI::call_minScan<double>(const double*, double*, int) const;
#endif
// Default instantiations of call_minScan(const type*, type*, int)
template <class type>
void AMP_MPI::call_minScan(const type *x, type *y, int n) const {
    char message[200];
    sprintf(message,"Default instantion of minScan in parallel is not supported (%s)",typeid(type).name());
    AMP_ERROR(message);
}


/************************************************************************
*  maxScan                                                              *
************************************************************************/
template <class type>
inline void AMP_MPI::maxScan(const type *x, type *y, const int n) const {
    if ( comm_size > 1 ) {
        call_maxScan(x,y,n);
    } else {
        for (int i=0; i<n; i++)
            y[i] = x[i];
    }
}
// Define specializations of call_maxScan(const type*, type*, int)
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_maxScan<unsigned char>(const unsigned char*, unsigned char*, int) const;
template <> void AMP_MPI::call_maxScan<char>(const char*, char*, int) const;
template <> void AMP_MPI::call_maxScan<unsigned int>(const unsigned int*, unsigned int*, int) const;
template <> void AMP_MPI::call_maxScan<int>(const int*, int*, int) const;
template <> void AMP_MPI::call_maxScan<unsigned long int>(const unsigned long int*, unsigned long int*, int) const;
template <> void AMP_MPI::call_maxScan<long int>(const long int*, long int*, int) const;
template <> void AMP_MPI::call_maxScan<float>(const float*, float*, int) const;
template <> void AMP_MPI::call_maxScan<double>(const double*, double*, int) const;
#endif
// Default instantiations of call_maxScan(const type*, type*, int)
template <class type>
void AMP_MPI::call_maxScan(const type *x, type *y, int n) const {
    char message[200];
    sprintf(message,"Default instantion of maxReduce in parallel is not supported (%s)",typeid(type).name());
    AMP_ERROR(message);
}


/************************************************************************
*  allToAll                                                             *
************************************************************************/
// Define specializations of allToAll(const int n, const char*, char* )
#ifdef USE_EXT_MPI
template <> void AMP_MPI::allToAll<unsigned char>(const int n, const unsigned char*, unsigned char* ) const;
template <> void AMP_MPI::allToAll<char>(const int n, const char*, char* ) const;
template <> void AMP_MPI::allToAll<unsigned int>(const int n, const unsigned int*, unsigned int* ) const;
template <> void AMP_MPI::allToAll<int>(const int n, const int*, int* ) const;
template <> void AMP_MPI::allToAll<unsigned long int>(const int n, const unsigned long int*, unsigned long int* ) const;
template <> void AMP_MPI::allToAll<long int>(const int n, const long int*, long int* ) const;
template <> void AMP_MPI::allToAll<float>(const int n, const float*, float* ) const;
template <> void AMP_MPI::allToAll<double>(const int n, const double*, double* ) const;
#endif
// Default instantiations of allToAll(const int n, const char*, char* )
#ifdef USE_EXT_MPI
template <class type>
void AMP_MPI::allToAll(const int n, const type *send_data, type *recv_data ) const {
    allToAll<char>(n*sizeof(type),(char*)send_data,(char*)recv_data);
}
#else
template <class type>
void AMP_MPI::allToAll(const int n, const type *send_data, type *recv_data ) const {
    if ( comm_size != 1 )
        AMP_ERROR("Invalid size for allToAll");
    for (int i=0; i<n; i++)
        recv_data[i] = send_data[i];
}
#endif


/************************************************************************
*  allToAll                                                             *
************************************************************************/
template <class type>
int AMP_MPI::allToAll(const type *send_data, const int send_cnt[], const int send_disp[], 
        type *recv_data, int *recv_cnt, int* recv_disp, bool known_recv)  const
{
    int N_recieved = 0;
    if ( comm_size == 1 ) {
        // Special case for single-processor communicators
        if ( known_recv ) {
            if ( recv_cnt[0]!=send_cnt[0] && send_cnt[0]>0 )
                AMP_ERROR("Single processor send/recv are different sizes");
        } else {
            if ( recv_cnt != NULL )
                recv_cnt[0] = send_cnt[0];
            if ( recv_disp != NULL )
                recv_disp[0] = send_disp[0];
        }
        for (int i=0; i<send_cnt[0]; i++)
            recv_data[i+recv_disp[0]] = send_data[i+send_disp[0]];
        N_recieved = send_cnt[0];
    } else if ( known_recv ) {
        // The recieve sizes are known
        AMP_ASSERT( recv_cnt!=NULL && recv_disp!=NULL);
        call_allToAll(send_data,send_cnt,send_disp,recv_data,recv_cnt,recv_disp);
        for (int i=0; i<comm_size; i++)
            N_recieved += recv_cnt[i];
    } else {
        // The recieve sizes are not known, we need to communicate that information first
        int *recv_cnt2 = recv_cnt;
        int *recv_disp2 = recv_disp;
        if ( recv_cnt==NULL ) 
            recv_cnt2 = new int[comm_size];
        if ( recv_disp==NULL ) 
            recv_disp2 = new int[comm_size];
        // Communicate the size we will be recieving from each processor
        allToAll<int>( 1, send_cnt, recv_cnt2 );
        recv_disp2[0] = 0;
        for (int i=1; i<comm_size; i++)
            recv_disp2[i] = recv_disp2[i-1]+recv_cnt2[i-1];
        // Send the data
        call_allToAll(send_data,send_cnt,send_disp,recv_data,recv_cnt2,recv_disp2);
        for (int i=0; i<comm_size; i++)
            N_recieved += recv_cnt2[i];
        if ( recv_cnt==NULL ) 
            delete [] recv_cnt2;
        if ( recv_disp==NULL ) 
            delete [] recv_disp2;
    }
    return N_recieved;
}
// Define specializations of call_allToAll(const type*, const int*, const int*, type*, const int*, const int*);
#ifdef USE_EXT_MPI
template <> void AMP_MPI::call_allToAll<unsigned char>(const unsigned char*, const int*, const int*, unsigned char*, const int*, const int*) const;
template <> void AMP_MPI::call_allToAll<char>(const char*, const int*, const int*, char*, const int*, const int*) const;
template <> void AMP_MPI::call_allToAll<unsigned int>(const unsigned int*, const int*, const int*, unsigned int*, const int*, const int*) const;
template <> void AMP_MPI::call_allToAll<int>(const int*, const int*, const int*, int*, const int*, const int*) const;
template <> void AMP_MPI::call_allToAll<unsigned long int>(const unsigned long int*, const int*, const int*, unsigned long int*, const int*, const int*) const;
template <> void AMP_MPI::call_allToAll<long int>(const long int*, const int*, const int*, long int*, const int*, const int*) const;
template <> void AMP_MPI::call_allToAll<float>(const float*, const int*, const int*, float*, const int*, const int*) const;
template <> void AMP_MPI::call_allToAll<double>(const double*, const int*, const int*, double*, const int*, const int*) const;
#else
template <> void AMP_MPI::call_allToAll<char>(const char*, const int*, const int*, char*, const int*, const int*) const;
#endif
// Default instantiations of call_allToAll(const type*, const int*, const int*, type*, const int*, const int*);
template <class type>
void AMP_MPI::call_allToAll(const type *send_data, const int send_cnt[], const int send_disp[], type *recv_data, const int *recv_cnt, const int* recv_disp)  const
{
    int *send_cnt2  = new int[comm_size];
    int *recv_cnt2  = new int[comm_size];
    int *send_disp2 = new int[comm_size];
    int *recv_disp2 = new int[comm_size];
    for (int i=0; i<comm_size; i++) {
        send_cnt2[i] = send_cnt[i]*sizeof(type);
        send_disp2[i] = send_disp[i]*sizeof(type);
        recv_cnt2[i] = recv_cnt[i]*sizeof(type);
        recv_disp2[i] = recv_disp[i]*sizeof(type);
    }
    call_allToAll<char>( (char*) send_data, send_cnt2, send_disp2, (char*) recv_data, recv_cnt2, recv_disp2 );
    delete [] send_cnt2;
    delete [] recv_cnt2;
    delete [] send_disp2;
    delete [] recv_disp2;
}



