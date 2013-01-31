#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <complex>

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/AMP_MPI.h"
#include "utils/ProfilerApp.h"
#include "utils/PIO.h"


struct mytype{
    int a;
    double b;
    mytype( ) {
        a = -1;
        b = -1.0;
    }
    mytype( int i ) {
        a = i;
        b = -1.0;
    }
    mytype( int i, double d ) {
        a = i;
        b = d;
    }
    bool operator==( mytype &other ) {
        if ( a==other.a && b==other.b )
            return true;
        return false;
    }
    bool operator!=( mytype &other ) {
        if ( a!=other.a || b!=other.b )
            return true;
        return false;
    }
};


// Routines to test Reduce with known data types
// flag - 0: all tests should pass
//        1: basic reduce should pass, reduce with rank should fail with error message
template <class type>
int testReduce(AMP::AMP_MPI comm, AMP::UnitTest *ut, int flag);
template <>
int testReduce<std::complex<double> >(AMP::AMP_MPI comm, AMP::UnitTest *ut, int flag) {
    PROFILE_START("testReduce<complex double>");
    char message[128];
    std::complex<double> rank = comm.getRank()+1;
    std::complex<double> N = ((comm.getSize()*(comm.getSize()+1))/2);
    // Test sumReduce
    sprintf(message,"sumReduce (%s)",typeid(std::complex<double>).name());
    if ( comm.sumReduce< std::complex<double> >(rank) == N )
        ut->passes(message);
    else
        ut->failure(message);
    sprintf(message,"sumReduce (%s) (x,y)",typeid(std::complex<double>).name());
    std::complex<double> y;
    comm.sumReduce< std::complex<double> >(&rank,&y,1);
    if ( y == N )
        ut->passes(message);
    else
        ut->failure(message);
    PROFILE_STOP("testReduce<complex double>");
    return 2;   // Return the number of tests
}
template <class type>
int testReduce(AMP::AMP_MPI comm, AMP::UnitTest *ut, int flag) {
    PROFILE_START("testReduce");
    char message[128];
    type rank = (type) comm.getRank();
    type size = (type) comm.getSize();
    if ( (int)(size) != comm.getSize() ) {
        sprintf(message,"Reduce (%s) cannot represent the number of processors",typeid(type).name());
        ut->expected_failure(message);
        PROFILE_STOP2("testReduce<class type>");
        return 0;
    }
    type x, y;
    int N = ((comm.getSize()*(comm.getSize()+1))/2);
    // Test sumReduce
    sprintf(message,"sumReduce (%s)",typeid(type).name());
    if ( ((int)((type)N)) != N )       
        ut->expected_failure(message);      // type cannot represent N
    else if ( comm.sumReduce<type>(rank+1) == (type) N )
        ut->passes(message);
    else
        ut->failure(message);
    sprintf(message,"sumReduce (%s) (x,y)",typeid(type).name());
    x = rank+1;
    comm.sumReduce<type>(&x,&y,1);
    if ( ((int)((type)N)) != N )    
        ut->expected_failure(message);
    else if ( y == (type) N )
        ut->passes(message);
    else
        ut->failure(message);
    // Test minReduce
    sprintf(message,"minReduce (%s)",typeid(type).name());    
    if ( comm.minReduce<type>(rank+1) == 1 )
        ut->passes(message);
    else
        ut->failure(message);
    sprintf(message,"minReduce (%s) (x,y)",typeid(type).name());    
    comm.minReduce<type>(&x,&y,1,NULL);
    if ( y == 1 )
        ut->passes(message);
    else
        ut->failure(message);
    // Test minReduce
    sprintf(message,"maxReduce (%s)",typeid(type).name());    
    if ( comm.maxReduce<type>(rank+1) == size )
        ut->passes(message);
    else
        ut->failure(message);
    sprintf(message,"maxReduce (%s) (x,y)",typeid(type).name());    
    comm.maxReduce<type>(&x,&y,1,NULL);
    if ( y == size )
        ut->passes(message);
    else
        ut->failure(message);
    // Test minReduce with rank
    int rank_of_min=-1;
    int rank_of_max=-1;
    type rank_min = rank+1;
    type rank_max = rank+1;
    sprintf(message,"minReduce-rank (%s)",typeid(type).name());
    try {
        comm.minReduce<type>(&rank_min,1,&rank_of_min);
        if ( rank_min==1 && rank_of_min==0 )
            ut->passes(message);
        else
            ut->failure(message);
        if ( flag==1 && comm.getSize()>1 )
            ut->failure(message);
    } catch (...) {
        if ( flag==1 && comm.getSize()>1 )
            ut->expected_failure(message);
        else
            ut->failure(message);
    }
    sprintf(message,"minReduce-rank (%s) (x,y)",typeid(type).name());
    try {
        comm.minReduce<type>(&x,&rank_min,1,&rank_of_min);
        if ( rank_min==1 && rank_of_min==0 )
            ut->passes(message);
        else
            ut->failure(message);
        if ( flag==1 && comm.getSize()>1 )
            ut->failure(message);
    } catch (...) {
        if ( flag==1 && comm.getSize()>1 )
            ut->expected_failure(message);
        else
            ut->failure(message);
    }
    // Test maxReduce with rank
    sprintf(message,"maxReduce-rank (%s)",typeid(type).name());    
    try {
        comm.maxReduce<type>(&rank_max,1,&rank_of_max);
        if ( rank_max==size && rank_of_max==comm.getSize()-1 )
            ut->passes(message);
        else
            ut->failure(message);
        if ( flag==1 && comm.getSize()>1 )
            ut->failure(message);
    } catch (...) {
        if ( flag==1 && comm.getSize()>1 )
            ut->expected_failure(message);
        else
            ut->failure(message);
    }
    sprintf(message,"maxReduce-rank (%s) (x,y)",typeid(type).name());    
    try {
        comm.maxReduce<type>(&x,&rank_max,1,&rank_of_max);
        if ( rank_max==size && rank_of_max==comm.getSize()-1 )
            ut->passes(message);
        else
            ut->failure(message);
        if ( flag==1 && comm.getSize()>1 )
            ut->failure(message);
    } catch (...) {
        if ( flag==1 && comm.getSize()>1 )
            ut->expected_failure(message);
        else
            ut->failure(message);
    }
    PROFILE_STOP("testReduce");
    return 10;   // Return the number of tests
}


// Routine to test Scan with known data types
// flag - 0: all tests should pass
//        1: only sumScan is valid (complex<double>)
template <class type>
int testScan(AMP::AMP_MPI comm, AMP::UnitTest *ut, int flag=0) {
    PROFILE_START("testScan");
    char message[500];
    type x = (type) (comm.getRank()+1);
    type y;
    sprintf(message,"sumScan (%s)",typeid(type).name());
    comm.sumScan<type>(&x,&y,1);
    type N = (type) (((comm.getRank()+1)*(comm.getRank()+2))/2);
    if ( y == N )
        ut->passes(message);
    else
        ut->failure(message);
    if ( flag==1 ) {
        PROFILE_STOP2("testScan");
        return 1;
    }
    sprintf(message,"minScan (%s)",typeid(type).name());    
    comm.minScan<type>(&x,&y,1);
    if ( y == (type) 1 )
        ut->passes(message);
    else
        ut->failure(message);
    sprintf(message,"maxScan (%s)",typeid(type).name());    
    comm.maxScan<type>(&x,&y,1);
    if ( y == x )
        ut->passes(message);
    else
        ut->failure(message);
    PROFILE_STOP("testScan");
    return 3;   // Return the number of tests
}


// Routine to test bcast
template <class type>
int testBcast(AMP::AMP_MPI comm, AMP::UnitTest *ut, type default_val, type new_val) {
    PROFILE_START("testBcast");
    char message[128];
    for (int i=0; i<comm.getSize(); i++) {
        type tmp1 = default_val;
        if ( comm.getRank() == i )
            tmp1 = new_val;
        sprintf(message,"bcast scalar (%s) from rank %i",typeid(type).name(),i);
        if ( comm.bcast(tmp1,i)==new_val )
            ut->passes(message);
        else
            ut->failure(message);
        type tmp2[2];
        tmp2[0] = default_val;
        tmp2[1] = default_val;
        if ( comm.getRank() == i ) {
            tmp2[0] = new_val;
            tmp2[1] = new_val;
        }
        sprintf(message,"bcast vector (%s) from rank %i",typeid(type).name(),i);
        comm.bcast(tmp2,2,i);
        if ( tmp2[0]==new_val && tmp2[1]==new_val )
            ut->passes(message);
        else
            ut->failure(message);
    }
    PROFILE_STOP("testBcast");
    return 2*comm.getSize();   // Return the number of tests
}


// Routine to test allGather
template <class type>
int testAllGather(AMP::AMP_MPI comm, AMP::UnitTest *ut) {
    PROFILE_START("testAllGather");
    char message[128];
    // Test scalar allGather
    type x1 = (type) comm.getRank();
    type *x2 = new type[comm.getSize()];
    comm.allGather(x1,x2);
    bool pass = true;
    for (int i=0; i<comm.getSize(); i++) {
        type test = i;
        if ( x2[i] != test )
            pass = false;
    }
    sprintf(message,"allGather scalar (%s)",typeid(type).name());
    if ( pass )
        ut->passes(message);
    else
        ut->failure(message);
    // Test vector allGather
    int N = (comm.getSize()*(comm.getSize()+1))/2;
    type *x3 = new type[comm.getRank()+1];
    type *x4 = new type[N];
    type *x5 = new type[N];
    int *size = new int[comm.getSize()];
    for (int i=0; i<=comm.getRank(); i++)
        x3[i] = (type) comm.getRank();
    int tot1 = comm.allGather(x3,comm.getRank()+1,x4);
    int tot2 = comm.allGather(x3,comm.getRank()+1,x5,size);
    pass = true;
    if ( tot1!=N || tot2!=N )
        pass = false;
    int k = 0;
    for (int i=0; i<comm.getSize(); i++) {
        if ( size[i] != i+1 )
            pass = false;
        if ( !pass ) 
            break;
        for (int j=0; j<=i; j++) {
            type test = i;
            if ( x4[k]!=test || x5[k]!=test )
                pass = false;
            k++;
        }
    }
    sprintf(message,"allGather vector (%s)",typeid(type).name());
    if ( pass )
        ut->passes(message);
    else
        ut->failure(message);
    delete [] x2;
    delete [] x3;
    delete [] x4;
    delete [] x5;
    delete [] size;
    // Test vector allGather with know recive sizes and non-zero displacements
    type *send = new type[comm.getRank()+1];
    type *recv = new type[comm.getSize()*comm.getSize()+1];
    int *recv_size = new int[comm.getSize()];
    int *recv_disp = new int[comm.getSize()];
    for (int i=0; i<=comm.getRank(); i++)
        send[i] = i;
    for (int i=0; i<comm.getSize(); i++)
        recv_size[i] = i+1;
    for (int i=0; i<comm.getSize(); i++)
        recv_disp[i] = 1 + i*comm.getSize() + comm.getSize()-i-1;
    for (int i=0; i<=comm.getSize()*comm.getSize(); i++)
        recv[i] = (type) -1;
    int tot = comm.allGather(send,comm.getRank()+1,recv,recv_size,recv_disp,true);
    pass = true;
    if ( tot != N )
        pass = false;
    type test = (type) -1;
    if ( recv[0] != test )
        pass = false;
    for (int i=0; i<comm.getSize(); i++) {
        for (int j=0; j<comm.getSize(); j++) {
            int k = j+i*comm.getSize()+1 - recv_disp[i];
            if ( k>=0 )
                test = k;
            else
                test = (type) -1;
            if ( recv[j+i*comm.getSize()+1] != test ) 
                pass = false;
        }
    }
    sprintf(message,"allGather vector with known recv and non-zero displacements (%s)",typeid(type).name());
    if ( pass )
        ut->passes(message);
    else
        ut->failure(message);
    delete [] send;
    delete [] recv;
    delete [] recv_size;
    delete [] recv_disp;
    // Test vector allGather with no elements
    size = new int[comm.getSize()];
    sprintf(message,"allGather scalar (%s)",typeid(type).name());
    try {
        comm.allGather(&x1,0,(type*)NULL,size);
        ut->passes(message);
    } catch (...) {
        ut->failure(message);
    }
    delete [] size;
    PROFILE_STOP("testAllGather");
    return 4;   // Return the number of tests
}


// Routine to test setGather
template <class type>
int testSetGather(AMP::AMP_MPI comm, AMP::UnitTest *ut) {
    PROFILE_START("testSetGather");
    char message[500];
    type x1 = (type) comm.getRank();
    std::set<type> set;
    set.insert( x1 );
    comm.setGather( set );
    bool pass = true;
    for (int i=0; i<comm.getSize(); i++) {
        type x2 = i;
        if ( set.find(x2)==set.end() )
            pass = false;
    }
    sprintf(message,"setGather (%s)",typeid(type).name());
    if ( pass )
        ut->passes(message);
    else
        ut->failure(message);
    PROFILE_STOP("testSetGather");
    return 1;   // Return the number of tests
}


// Routine to test mapGather
template <class type>
int testMapGather(AMP::AMP_MPI comm, AMP::UnitTest *ut) {
    PROFILE_START("testMapGather");
    char message[128];
    type x1 = (type) comm.getRank();
    std::map<int,type> map;
    map.insert( std::pair<int,type>( comm.getRank(), x1 ) );
    comm.mapGather( map );
    bool pass = true;
    for (int i=0; i<comm.getSize(); i++) {
        type x2 = i;
        typename std::map<int,type>::iterator it;
        it = map.find(i);
        if ( it==map.end() )
            pass = false;
        else if ( it->second != x2 )
            pass = false;
    }
    sprintf(message,"mapGather (%s)",typeid(type).name());
    if ( pass )
        ut->passes(message);
    else
        ut->failure(message);
    PROFILE_STOP("testMapGather");
    return 1;   // Return the number of tests
}


// Routine to test allToAll
template <class type>
int testAllToAll(AMP::AMP_MPI comm, AMP::UnitTest *ut) {
    PROFILE_START("testAllToAll");
    bool pass;
    char message[128];
    int size = 0;    
    type *send_data, *recv_data;
    int *send_cnt = new int[comm.getSize()];
    int *recv_cnt = new int[comm.getSize()];
    int *send_disp = new int[comm.getSize()];
    int *recv_disp = new int[comm.getSize()];
    // Test allToAll with a scalar value to each processor
    send_data = new type[comm.getSize()];
    recv_data = new type[comm.getSize()];
    for (int i=0; i<comm.getSize(); i++)
        send_data[i] = comm.getSize();
    comm.allToAll(1,send_data,recv_data);
    pass = true;
    for (int i=0; i<comm.getSize(); i++) {
        type test = comm.getSize();
        if ( recv_data[i] != test ) 
            pass = false;
    }
    delete [] send_data;
    delete [] recv_data;
    sprintf(message,"allToAll with scalar (%s)",typeid(type).name());
    if ( pass )
        ut->passes(message);
    else
        ut->failure(message);
    // Test allToAll vector with a scalar value to each processor
    send_data = new type[comm.getSize()];
    recv_data = new type[comm.getSize()];
    for (int i=0; i<comm.getSize(); i++) {
        send_cnt[i] = 1;
        recv_cnt[i] = 1;
        send_disp[i] = i;
        recv_disp[i] = i;
        send_data[i] = comm.getSize();
        recv_data[i] = 0;
    }
    size = comm.allToAll(send_data,send_cnt,send_disp,recv_data,recv_cnt,recv_disp,true);
    pass = true;
    if ( size!=comm.getSize() )
        pass = false;
    for (int i=0; i<comm.getSize(); i++) {
        type test = comm.getSize();
        if ( recv_data[i] != test ) 
            pass = false;
    }
    delete [] send_data;
    delete [] recv_data;
    sprintf(message,"allToAll vector with scalar (%s)",typeid(type).name());
    if ( pass )
        ut->passes(message);
    else
        ut->failure(message);
    // Test allToAll with a variable number of values per processor and spacing
    send_data = new type[comm.getSize()*comm.getSize()];
    recv_data = new type[2*comm.getRank()*comm.getSize()];
    for (int i=0; i<comm.getSize(); i++) {
        send_cnt[i] = i;
        recv_cnt[i] = comm.getRank();
        send_disp[i] = i*comm.getSize();
        recv_disp[i] = 2*i*comm.getRank();
        for (int j=0; j<comm.getSize(); j++) {
            if ( j < i )
                send_data[j+send_disp[i]] = i;
            else
                send_data[j+send_disp[i]] = (type) -1;
        }
    }
    for (int i=0; i<2*comm.getRank()*comm.getSize(); i++)
        recv_data[i] = (type) -2;
    size = comm.allToAll(send_data,send_cnt,send_disp,recv_data,recv_cnt,recv_disp,true);
    pass = true;
    if ( size!=comm.getRank()*comm.getSize() )
        pass = false;
    for (int i=0; i<comm.getSize(); i++) {
        for (int j=0; j<2*comm.getRank(); j++) {
            if ( j < comm.getRank() ) {
                type test = comm.getRank();
                if ( recv_data[j+recv_disp[i]] != test )
                    pass = false;
            } else {
                type test = (type) -2;
                if ( recv_data[j+recv_disp[i]] != test )
                    pass = false;
            }
        }
    }
    delete [] send_data;
    delete [] recv_data;
    sprintf(message,"allToAll with vector of known size and displacements (%s)",typeid(type).name());
    if ( pass )
        ut->passes(message);
    else
        ut->failure(message);
    // Test allToAll with a unknown recieve length
    send_data = new type[comm.getSize()*comm.getSize()];
    type *recv_data1 = new type[comm.getSize()*comm.getSize()];
    type *recv_data2 = new type[comm.getSize()*comm.getSize()];
    for (int i=0; i<comm.getSize(); i++) {
        send_cnt[i] = i;
        recv_cnt[i] = -1;
        send_disp[i] = i*comm.getSize();
        recv_disp[i] = -1;
        for (int j=0; j<comm.getSize(); j++) {
            if ( j < i ) 
                send_data[j+send_disp[i]] = i;
            else
                send_data[j+send_disp[i]] = (type) -1;
        }
    }
    for (int i=0; i<comm.getSize()*comm.getSize(); i++) {
        recv_data1[i] = (type) -2;
        recv_data2[i] = (type) -2;
    }
    int size1 = comm.allToAll(send_data,send_cnt,send_disp,recv_data1,recv_cnt,recv_disp,false);
    int size2 = comm.allToAll(send_data,send_cnt,send_disp,recv_data2);
    bool pass1 = true;
    bool pass2 = true;
    if ( size1!=comm.getRank()*comm.getSize() )
        pass1 = false;
    if ( size2!=comm.getRank()*comm.getSize() )
        pass2 = false;
    for (int i=0; i<comm.getSize(); i++) {
        if ( recv_cnt[i]!=comm.getRank() || recv_disp[i]!=i*comm.getRank() )
            pass1 = false;
    }
    for (int i=0; i<comm.getRank()*comm.getSize(); i++) {
        type test = comm.getRank();
        if ( recv_data1[i] != test )
            pass1 = false;
        if ( recv_data2[i] != test )
            pass2 = false;
    }
    delete [] send_data;
    delete [] recv_data1;
    delete [] recv_data2;
    sprintf(message,"allToAll with vector of unknown size (%s)",typeid(type).name());
    if ( pass1 )
        ut->passes(message);
    else
        ut->failure(message);
    sprintf(message,"allToAll with vector of unknown size with NULL recv(%s)",typeid(type).name());
    if ( pass2 )
        ut->passes(message);
    else
        ut->failure(message);
    // Free temporary variables
    delete [] send_cnt;
    delete [] recv_cnt;
    delete [] send_disp;
    delete [] recv_disp;
    PROFILE_STOP("testAllToAll");
    return 5;   // Return the number of tests

}


// Routine to test send/recv
template <class type>
int testSendRecv(AMP::AMP_MPI comm, AMP::UnitTest *ut, type v1, type v2) {
    PROFILE_START("testSendRecv");
    char message[128];
    // Test send-recv with a known length
    for (int i=0; i<comm.getSize(); i++) {
        for (int j=0; j<comm.getSize(); j++) {
            type x = v1;
            int tag = i+j*comm.getSize();
            sprintf(message,"send-recv %i-%i known length (%s)",i,j,typeid(type).name());
            if ( i==j ) {
                // We are not allowed to send/recieve from the same processor
                continue;
            } else if ( i==comm.getRank() ) {
                // We are sending
                x = v2;
                comm.send(&x,1,j,tag);
            } else if ( j==comm.getRank() ) {
                // We are recieving
                int size=1;
                comm.recv(&x,size,i,false,tag);
                if ( size==1 && x==v2 )
                    ut->passes(message);
                else
                    ut->failure(message);
            }
        }
    }
    // Test send-recv with an unknown length
    for (int i=0; i<comm.getSize(); i++) {
        for (int j=0; j<comm.getSize(); j++) {
            type x = v1;
            int tag = i+j*comm.getSize();
            sprintf(message,"send-recv %i-%i unknown length (%s)",i,j,typeid(type).name());
            if ( i==j ) {
                // We are not allowed to send/recieve from the same processor
                continue;
            } else if ( i==comm.getRank() ) {
                // We are sending
                x = v2;
                comm.send(&x,1,j,tag);
            } else if ( j==comm.getRank() ) {
                // We are recieving
                int size=1;
                comm.recv(&x,size,i,true,tag);
                if ( size==1 && x==v2 )
                    ut->passes(message);
                else
                    ut->failure(message);
            }
        }
    }
    // Test send-recv with an empty length
    for (int i=0; i<comm.getSize(); i++) {
        for (int j=0; j<comm.getSize(); j++) {
            type x = v1;
            int tag = i+j*comm.getSize();
            sprintf(message,"send-recv %i-%i empty length (%s)",i,j,typeid(type).name());
            if ( i==j ) {
                // We are not allowed to send/recieve from the same processor
                continue;
            } else if ( i==comm.getRank() ) {
                // We are sending
                x = v2;
                comm.send(&x,0,j,tag);
            } else if ( j==comm.getRank() ) {
                // We are recieving
                int size = comm.probe(i,tag);
                comm.recv(&x,size,i,false,tag);
                if ( size==0 )
                    ut->passes(message);
                else
                    ut->failure(message);
            }
        }
    }
    PROFILE_STOP("testSendRecv");
    return 3*comm.getSize()*comm.getSize();   // Return the number of tests
}


// Routine to test Isend/Irecv
template <class type>
int testIsendIrecv(AMP::AMP_MPI comm, AMP::UnitTest *ut, type v1, type v2) {
    PROFILE_START("testIsendIrecv");
    char message[128];
    std::vector<MPI_Request> sendRequest;
    std::vector<MPI_Request> recvRequest;
    // Send all messages
    for (int i=0; i<comm.getSize(); i++) {
        // Check if the current rank is sending
        if ( i!=comm.getRank() )
            continue;
        for (int j=0; j<comm.getSize(); j++) {
            // Start a non-blocking send
            int tag = i+j*comm.getSize();
            MPI_Request request = comm.Isend(&v1,1,j,tag);
            sendRequest.insert(sendRequest.begin(),request);
        }
    }
    // Recv all messages
    type *recv_buffer = new type[comm.getSize()];
    for (int i=0; i<comm.getSize(); i++)
        recv_buffer[i] = v2;
    recv_buffer[comm.getRank()] = v1;
    for (int j=0; j<comm.getSize(); j++) {
        // Check if the current rank is recieving
        if ( j!=comm.getRank() )
            continue;
        for (int i=0; i<comm.getSize(); i++) {
            // Start a non-blocking recv
            int tag = i+j*comm.getSize();
            MPI_Request request = comm.Irecv(&recv_buffer[i],1,i,tag);
            recvRequest.insert(recvRequest.begin(),request);
        }
    }
    // Wait for all communications to finish
    AMP::AMP_MPI::wait(sendRequest[0]);
    sendRequest.erase(sendRequest.begin()+0);
    while ( sendRequest.size() > 0 ) {
        int index = comm.waitAny(sendRequest.size(),&(sendRequest[0]));
        sendRequest.erase(sendRequest.begin()+index);
    }
    AMP::AMP_MPI::waitAll(recvRequest.size(),&(recvRequest[0]));
    // Check the recieved values
    bool pass = true;
    for (int i=0; i<comm.getSize(); i++) {
        if ( recv_buffer[i] != v1 )
            pass = false;
    }
    sprintf(message,"Isend-Irecv (%s)",typeid(type).name());
    if ( pass )
        ut->passes(message);
    else
        ut->failure(message);
    delete [] recv_buffer;
    PROFILE_STOP("testIsendIrecv");
    return comm.getSize()*comm.getSize();   // Return the number of tests
}


// Structure to contain timer results
struct testCommTimerResults {
    int N_reduce;
    int N_scan;
    int N_bcast;
    int N_allGather;
    int N_setGather;
    int N_mapGather;
    int N_allToAll;
    int N_sendRecv;
    int N_IsendIrecv;
    double t_reduce;
    double t_scan;
    double t_bcast;
    double t_allGather;
    double t_setGather;
    double t_mapGather;
    double t_allToAll;
    double t_sendRecv;
    double t_IsendIrecv;
    // Constructor
    testCommTimerResults() {
        N_reduce = 0;
        N_scan = 0;
        N_bcast = 0;
        N_allGather = 0;
        N_setGather = 0;
        N_mapGather = 0;
        N_allToAll = 0;
        N_sendRecv = 0;
        N_IsendIrecv = 0;
        t_reduce = 0.0;
        t_scan = 0.0;
        t_bcast = 0.0;
        t_allGather = 0.0;
        t_setGather = 0.0;
        t_mapGather = 0.0;
        t_allToAll = 0.0;
        t_sendRecv = 0.0;
        t_IsendIrecv = 0.0;
    }
    // Print the results
    void print() {
        printf("   Reduce:      N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",N_reduce,t_reduce,1e6*t_reduce/N_reduce);
        printf("   Scan:        N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",N_scan,t_scan,1e6*t_scan/N_scan);
        printf("   Bcast:       N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",N_bcast,t_bcast,1e6*t_bcast/N_bcast);
        printf("   allGather:   N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",N_allGather,t_allGather,1e6*t_allGather/N_allGather);
        printf("   allToAll:    N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",N_allToAll,t_allToAll,1e6*t_allToAll/N_allToAll);
        printf("   send-recv:   N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",N_sendRecv,t_sendRecv,1e6*t_sendRecv/N_sendRecv);
        printf("   Isend-Irecv: N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",N_IsendIrecv,t_IsendIrecv,1e6*t_IsendIrecv/N_IsendIrecv);
        printf("   setGather:   N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",N_setGather,t_setGather,1e6*t_setGather/N_setGather);
        printf("   mapGather:   N = %5i, t_tot = %0.5e, t_avg = %6.1f us\n",N_mapGather,t_mapGather,1e6*t_mapGather/N_mapGather);
    }
};


// This routine will test a single MPI communicator
testCommTimerResults testComm(AMP::AMP_MPI comm, AMP::UnitTest *ut) {
    PROFILE_START("testComm");
    testCommTimerResults timer;
    double start_time;
    // Test the tag
    int tag0 = comm.newTag();
    bool pass = tag0>0 && tag0<comm.maxTag();
    for (int i=1; i<255; i++) {
        if ( comm.newTag()!=tag0+i )
            pass = false;
    }
    if ( pass ) 
        ut->passes("newTag");
    else
        ut->failure("newTag");
    // Test all and any reduce
    bool test1 = !comm.allReduce(comm.getRank()!=0);
    bool test2 = comm.allReduce(true);
    if ( test1 && test2 )
        ut->passes("allReduce");
    else
        ut->failure("allReduce");
    test1 = comm.anyReduce(comm.getRank()==0);
    test2 = !comm.anyReduce(false);
    if ( test1 && test2 )
        ut->passes("anyReduce");
    else
        ut->failure("anyReduce");
    // Test min, max, and sum reduce
    start_time = AMP::AMP_MPI::time();
    timer.N_reduce += testReduce<unsigned char>(comm,ut,1);         // does not support rank of min/max
    timer.N_reduce += testReduce<char>(comm,ut,1);                  // does not support rank of min/max
    timer.N_reduce += testReduce<unsigned int>(comm,ut,1);          // does not support rank of min/max
    timer.N_reduce += testReduce<int>(comm,ut,0);
    timer.N_reduce += testReduce<unsigned long int>(comm,ut,1);     // does not support rank of min/max
    timer.N_reduce += testReduce<long int>(comm,ut,0);
    timer.N_reduce += testReduce<size_t>(comm,ut,1);                // does not support rank of min/max
    timer.N_reduce += testReduce<float>(comm,ut,0);
    timer.N_reduce += testReduce<double>(comm,ut,0);
    timer.N_reduce += testReduce<std::complex<double> >(comm,ut,2); // only sumreduce is valid for complex numbers
    mytype tmp1(1,-1.0);
    mytype tmp2;
    if ( comm.getSize() > 1 ) {
        // We can't perform a reduce on an unknown data type (this should throw an error)
        try {
            // This should fail
            tmp2 = comm.sumReduce<mytype>(tmp1);
            ut->failure("sumReduce should give an error with an unknown type");
        } catch (...) {
            ut->passes("sumReduce should give an error with an unknown type");
        }
        try {
            // This should fail
            tmp2 = comm.minReduce<mytype>(tmp1);
            ut->failure("minReduce should give an error with an unknown type");
        } catch (...) {
            ut->passes("minReduce should give an error with an unknown type");
        }
        try {
            // This should fail
            tmp2 = comm.maxReduce<mytype>(tmp1);
            ut->failure("maxReduce should give an error with an unknown type");
        } catch (...) {
            ut->passes("maxReduce should give an error with an unknown type");
        }
        timer.N_reduce += 3;
    }
    timer.t_reduce = AMP::AMP_MPI::time()-start_time;
    // Test min, max, and sum scan
    start_time = AMP::AMP_MPI::time();
    timer.N_scan += testScan<unsigned char>(comm,ut);
    timer.N_scan += testScan<char>(comm,ut);
    timer.N_scan += testScan<unsigned int>(comm,ut);
    timer.N_scan += testScan<int>(comm,ut);
    timer.N_scan += testScan<unsigned long int>(comm,ut);
    timer.N_scan += testScan<long int>(comm,ut);
    timer.N_scan += testScan<size_t>(comm,ut);
    timer.N_scan += testScan<float>(comm,ut);
    timer.N_scan += testScan<double>(comm,ut);
    timer.N_scan += testScan< std::complex<double> >(comm,ut,1);    // Only sumScan is valid with complex data
    if ( comm.getSize() > 1 ) {
        // We can't perform a reduce on an unknown data type (this should throw an error)
        try {
            // This should fail
            comm.sumScan<mytype>(&tmp1,&tmp2,1);
            ut->failure("sumReduce should give an error with an unknown type");
        } catch (...) {
            ut->passes("sumReduce should give an error with an unknown type");
        }
        try {
            // This should fail
            comm.minScan<mytype>(&tmp1,&tmp2,1);
            ut->failure("minReduce should give an error with an unknown type");
        } catch (...) {
            ut->passes("minReduce should give an error with an unknown type");
        }
        try {
            // This should fail
            comm.maxScan<mytype>(&tmp1,&tmp2,1);
            ut->failure("maxReduce should give an error with an unknown type");
        } catch (...) {
            ut->passes("maxReduce should give an error with an unknown type");
        }
        timer.N_scan += 3;
    }
    timer.t_scan = AMP::AMP_MPI::time()-start_time;
    // Test bcast
    start_time = AMP::AMP_MPI::time();
    timer.N_bcast += testBcast<unsigned char>(comm,ut,0,1);
    timer.N_bcast += testBcast<char>(comm,ut,-1,1);
    timer.N_bcast += testBcast<unsigned int>(comm,ut,0,1);
    timer.N_bcast += testBcast<int>(comm,ut,-1,1);
    timer.N_bcast += testBcast<unsigned long int>(comm,ut,0,1);
    timer.N_bcast += testBcast<long int>(comm,ut,-1,1);
    timer.N_bcast += testBcast<size_t>(comm,ut,-1,1);
    timer.N_bcast += testBcast<float>(comm,ut,-1.0,1.0);
    timer.N_bcast += testBcast<double>(comm,ut,-1.0,1.0);
    mytype tmp3(-1,-1.0);
    mytype tmp4(1,1.0);
    timer.N_bcast += testBcast<mytype>(comm,ut,tmp3,tmp4);
    timer.t_bcast = AMP::AMP_MPI::time()-start_time;
    // Test barrier
    comm.barrier();
    // Test gather
    start_time = AMP::AMP_MPI::time();
    timer.N_allGather += testAllGather<unsigned char>(comm,ut);
    timer.N_allGather += testAllGather<char>(comm,ut);
    timer.N_allGather += testAllGather<unsigned int>(comm,ut);
    timer.N_allGather += testAllGather<int>(comm,ut);
    timer.N_allGather += testAllGather<unsigned long int>(comm,ut);
    timer.N_allGather += testAllGather<long int>(comm,ut);
    timer.N_allGather += testAllGather<size_t>(comm,ut);
    timer.N_allGather += testAllGather<float>(comm,ut);
    timer.N_allGather += testAllGather<double>(comm,ut);
    timer.N_allGather += testAllGather< std::complex<double> >(comm,ut);
    timer.N_allGather += testAllGather<mytype>(comm,ut);
    timer.t_allGather = AMP::AMP_MPI::time()-start_time;
    // Test std::set gather
    start_time = AMP::AMP_MPI::time();
    timer.N_setGather += testSetGather<unsigned char>(comm,ut);
    timer.N_setGather += testSetGather<char>(comm,ut);
    timer.N_setGather += testSetGather<unsigned int>(comm,ut);
    timer.N_setGather += testSetGather<int>(comm,ut);
    timer.N_setGather += testSetGather<unsigned long int>(comm,ut);
    timer.N_setGather += testSetGather<long int>(comm,ut);
    timer.N_setGather += testSetGather<size_t>(comm,ut);
    timer.N_setGather += testSetGather<float>(comm,ut);
    timer.N_setGather += testSetGather<double>(comm,ut);
    timer.t_setGather = AMP::AMP_MPI::time()-start_time;
    // Test std::map gather
    start_time = AMP::AMP_MPI::time();
    timer.N_mapGather += testMapGather<unsigned char>(comm,ut);
    timer.N_mapGather += testMapGather<char>(comm,ut);
    timer.N_mapGather += testMapGather<unsigned int>(comm,ut);
    timer.N_mapGather += testMapGather<int>(comm,ut);
    timer.N_mapGather += testMapGather<unsigned long int>(comm,ut);
    timer.N_mapGather += testMapGather<long int>(comm,ut);
    timer.N_mapGather += testMapGather<size_t>(comm,ut);
    timer.N_mapGather += testMapGather<float>(comm,ut);
    timer.N_mapGather += testMapGather<double>(comm,ut);
    timer.t_mapGather = AMP::AMP_MPI::time()-start_time;
    // Test allToAlll
    start_time = AMP::AMP_MPI::time();
    timer.N_allToAll += testAllToAll<unsigned char>(comm,ut);
    timer.N_allToAll += testAllToAll<char>(comm,ut);
    timer.N_allToAll += testAllToAll<unsigned int>(comm,ut);
    timer.N_allToAll += testAllToAll<int>(comm,ut);
    timer.N_allToAll += testAllToAll<unsigned long int>(comm,ut);
    timer.N_allToAll += testAllToAll<long int>(comm,ut);
    timer.N_allToAll += testAllToAll<size_t>(comm,ut);
    timer.N_allToAll += testAllToAll<float>(comm,ut);
    timer.N_allToAll += testAllToAll<double>(comm,ut);
    timer.N_allToAll += testAllToAll< std::complex<double> >(comm,ut);
    timer.N_allToAll += testAllToAll<mytype>(comm,ut);
    timer.t_allToAll = AMP::AMP_MPI::time()-start_time;
    // Test send/recv
    start_time = AMP::AMP_MPI::time();
    timer.N_sendRecv += testSendRecv<unsigned char>(comm,ut,0,1);
    timer.N_sendRecv += testSendRecv<char>(comm,ut,-1,1);
    timer.N_sendRecv += testSendRecv<unsigned int>(comm,ut,0,1);
    timer.N_sendRecv += testSendRecv<int>(comm,ut,-1,1);
    timer.N_sendRecv += testSendRecv<unsigned long int>(comm,ut,0,1);
    timer.N_sendRecv += testSendRecv<long int>(comm,ut,-1,1);
    timer.N_sendRecv += testSendRecv<size_t>(comm,ut,0,1);
    timer.N_sendRecv += testSendRecv<float>(comm,ut,-1.0,1.0);
    timer.N_sendRecv += testSendRecv<double>(comm,ut,-1.0,1.0);
    timer.N_sendRecv += testSendRecv<mytype>(comm,ut,tmp3,tmp4);
    timer.t_sendRecv = AMP::AMP_MPI::time()-start_time;
    // Test Isend/Irecv
    start_time = AMP::AMP_MPI::time();
    timer.N_IsendIrecv += testIsendIrecv<unsigned char>(comm,ut,0,1);
    timer.N_IsendIrecv += testIsendIrecv<char>(comm,ut,-1,1);
    timer.N_IsendIrecv += testIsendIrecv<unsigned int>(comm,ut,0,1);
    timer.N_IsendIrecv += testIsendIrecv<int>(comm,ut,-1,1);
    timer.N_IsendIrecv += testIsendIrecv<unsigned long int>(comm,ut,0,1);
    timer.N_IsendIrecv += testIsendIrecv<long int>(comm,ut,-1,1);
    timer.N_IsendIrecv += testIsendIrecv<size_t>(comm,ut,0,1);
    timer.N_IsendIrecv += testIsendIrecv<float>(comm,ut,-1.0,1.0);
    timer.N_IsendIrecv += testIsendIrecv<double>(comm,ut,-1.0,1.0);
    timer.N_IsendIrecv += testIsendIrecv<mytype>(comm,ut,tmp3,tmp4);
    timer.t_IsendIrecv = AMP::AMP_MPI::time()-start_time;
    PROFILE_STOP("testComm");
    return timer;
}


// Test comm dup and the number of communicators that can be created
void testCommDup(AMP::UnitTest *ut) {
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    AMP::AMP_MPI dupComm = globalComm.dup();
    if ( globalComm.getCommunicator()!=dupComm.getCommunicator() &&
        dupComm.getSize()==globalComm.getSize() && dupComm.getRank()==globalComm.getRank() ) 
    {
        ut->passes("dup comm");
    } else {
        ut->failure("dup comm");
        return;
    }
    size_t N_comm_try = 10000;  // Maximum number of comms to try and create
    std::vector<AMP::AMP_MPI> comms;
    comms.reserve(N_comm_try);
    try {
        for (size_t i=0; i<N_comm_try; i++) {
            comms.push_back( globalComm.dup() );
            AMP_ASSERT(globalComm.getCommunicator()!=comms[i].getCommunicator());
        }
        ut->passes("Created an unlimited number of comms");
    } catch (...) {
        if ( comms.size() < 252 ) {
            ut->failure("Could not create 252 different communicators");
        } else {
            char message[128];
            sprintf(message,"Failed to create an unlimited number of comms (%llu)",
                static_cast<unsigned long long int>(comms.size()));
            ut->expected_failure(message);
        }
        AMP::pout << "Maximum number of concurrent communicators: " << comms.size() << std::endl;
    }
    comms = std::vector<AMP::AMP_MPI>();
    size_t N_dup = 0;
    globalComm.barrier();
    try {
        double start = AMP::AMP_MPI::time();
        for (size_t i=0; i<N_comm_try; i++) {
            AMP::AMP_MPI tmp_comm = globalComm.dup();
            AMP_ASSERT(globalComm.getCommunicator()!=tmp_comm.getCommunicator());
            N_dup++;
        }
        double stop = AMP::AMP_MPI::time();
        ut->passes("Created/Destroyed an unlimited number of comms");
        char message[128];
        sprintf(message,"Time to create/destroy comm using AMP::AMP_MPI::dup() is: %0.1f us",1e6*(stop-start)/N_dup);
        AMP::pout << message << std::endl;
    } catch (...) {
        ut->failure("Failed to create/destroy an unlimited number of comms");
        AMP::pout << "Maximum number of communicators created with destruction: " << N_dup << std::endl;
    }
}


//  This test will test the AMP_MPI routines
int main(int argc, char *argv[])
{
    // Startup
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    int num_failed = 0;
    PROFILE_ENABLE(0);

    // Limit the scope so objects are destroyed
    {
        PROFILE_START("Main");
        // Create the unit test
        AMP::UnitTest ut;

        // Get the start time for the tests
        double start_time = AMP::AMP_MPI::time();

        // Print the global size (if we are using MPI)
        int global_size = 0;
        #ifdef USE_EXT_MPI
            MPI_Comm_size(MPI_COMM_WORLD,&global_size);
        #else
            global_size = 1;
        #endif

        // Test the global communicator (AMP_COMM_WORLD)
        AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
        if ( !globalComm.isNull() )
            ut.passes("Global communicator created");
        else
            ut.failure("Global communicator created");
        if ( globalComm.getSize()==global_size )
            ut.passes("Global communicator size");
        else
            ut.failure("Global communicator size");
        if ( globalComm.getRank()==0 )  {
            std::cout << "MPI_COMM_WORLD = " << global_size << " processors" << std::endl;
            std::cout << "   Largest tag value = " << globalComm.maxTag() << std::endl << std::endl;
        }
        #ifdef USE_EXT_MPI
            if ( globalComm.getCommunicator() == MPI_COMM_WORLD )
                ut.passes("Communicator == MPI_COMM_WORLD");
            else
                ut.failure("Communicator == MPI_COMM_WORLD");
        #endif
        testCommTimerResults commTimer = testComm(globalComm,&ut);
        if ( globalComm.getRank()==0 ) {
            std::cout << "Results for global timer (rank 0)" << std::endl;
            commTimer.print();
            std::cout << std::endl;
        }

        // Test AMP_COMM_SELF
        AMP::AMP_MPI selfComm = AMP::AMP_MPI(AMP_COMM_SELF);
        if ( !selfComm.isNull() )
            ut.passes("Self communicator created");
        else
            ut.failure("Self communicator created");
        #ifdef USE_EXT_MPI
            if ( selfComm.getCommunicator() == MPI_COMM_SELF )
                ut.passes("Communicator == MPI_COMM_SELF");
            else
                ut.failure("Communicator == MPI_COMM_SELF");
        #endif
        testComm(selfComm,&ut);
        
        // Test == and !=
        if ( globalComm==globalComm && !(selfComm==globalComm) )
            ut.passes("==");
        else
            ut.failure("==");
        if ( selfComm!=globalComm && !(globalComm!=globalComm) )
            ut.passes("!=");
        else
            ut.failure("!=");

        // Test AMP_COMM_NULL
        AMP::AMP_MPI nullComm = AMP::AMP_MPI(AMP_COMM_NULL);
        if ( nullComm.isNull() )
            ut.passes("Null communicator created");
        else
            ut.failure("Null communicator created");
        if ( nullComm.getSize()==0 )
            ut.passes("Null communicator has zero size");
        else
            ut.failure("Null communicator has zero size");
        #ifdef USE_EXT_MPI
            if ( nullComm.getCommunicator() == MPI_COMM_NULL )
                ut.passes("Communicator == MPI_COMM_NULL");
            else
                ut.failure("Communicator == MPI_COMM_NULL");
        #endif

        // Test dup
        testCommDup(&ut);
        AMP::AMP_MPI dupComm = globalComm.dup();
        
        // Test compare
        if ( globalComm.compare(globalComm)==1 )
            ut.passes("compare comm global==global");
        else
            ut.failure("compare comm global==global");
        if ( globalComm.compare(dupComm)==3 )
            ut.passes("compare comm global~=dup");
        else
            ut.failure("compare comm global~=dup");
        if ( global_size==1 ) {
            if ( globalComm.compare(selfComm)==3 )
                ut.passes("compare comm global~=self (global size=1)");
            else
                ut.failure("compare comm global~=self (global size=1)");
        } else {
            if ( globalComm.compare(selfComm)==0 )
                ut.passes("compare comm global!=self");
            else
                ut.failure("compare comm global!=self");
        }

        // Split the global comm and test
        PROFILE_START("Split");
        int color;
        if ( globalComm.getRank()==0 )
            color = 0;
        else if ( globalComm.getRank()<3 )
            color = 1;
        else 
            color = 2+(globalComm.getRank()-2)/4;
        std::vector<AMP::AMP_MPI> splitComms(4);
        splitComms[0] = globalComm.split( color );
        splitComms[1] = globalComm.split( color, globalComm.getRank() );
        if ( splitComms[0].getCommunicator()!=globalComm.getCommunicator() && 
             splitComms[1].getCommunicator()!=globalComm.getCommunicator() && 
             splitComms[0].getCommunicator()!=splitComms[1].getCommunicator() )
            ut.passes("split comm has different communicator");
        else
            ut.failure("split comm has different communicator");
        if ( globalComm.getSize()>1 ) {
            if ( splitComms[0].getSize()<globalComm.getSize() )
                ut.passes("split comm is smaller");
            else
                ut.failure("split comm is smaller");
        }
        if ( splitComms[0].getRank()==splitComms[1].getRank() )
            ut.passes("split sort by rank");
        else
            ut.failure("split sort by rank");
        testComm(splitComms[0],&ut);
        splitComms[2] = globalComm.split( -1 );
        if ( splitComms[2].isNull() )
            ut.passes("split with color=-1 returns NULL communicator");
        else
            ut.failure("split with color=-1 returns NULL communicator");
        if ( globalComm.getRank()==0 )
            color = -1;
        splitComms[3] = splitComms[0];  // Make a copy to ensure there are no memory leaks
        splitComms[3] = splitComms[2];  // Perform assignement to check memory leaks
        AMP_ASSERT(splitComms[3]==splitComms[2]);
        PROFILE_STOP("Split");

        // Test  <  <=  >  >=
        if ( globalComm.getSize()>1 ) {
            if ( splitComms[0]<globalComm && splitComms[1]<globalComm && !(globalComm<globalComm) && !(globalComm<splitComms[0]) ) 
                ut.passes(" < comm");
            else
                ut.failure(" < comm");
            if ( splitComms[0]<=globalComm && splitComms[1]<=globalComm && globalComm<=globalComm && !(globalComm<=splitComms[0]) ) 
                ut.passes(" <= comm");
            else
                ut.failure(" <= comm");
            if ( globalComm>splitComms[0] && globalComm>splitComms[1] && !(globalComm>globalComm) && !(splitComms[0]>globalComm) ) 
                ut.passes(" > comm");
            else
                ut.failure(" > comm");
            if ( globalComm>=splitComms[0] && globalComm>=splitComms[1] && globalComm>=globalComm && !(splitComms[0]>=globalComm) ) 
                ut.passes(" >= comm");
            else
                ut.failure(" >= comm");
        }

        // Test intersection
        // Test globalComm with selfComm
        if ( globalComm.getSize() > 1 ) {
            AMP::AMP_MPI comm1 = AMP::AMP_MPI::intersect( globalComm, selfComm );
            AMP::AMP_MPI comm2 = AMP::AMP_MPI::intersect( selfComm, globalComm );
            AMP::AMP_MPI comm3 = AMP::AMP_MPI::intersect( globalComm, globalComm );
            if ( comm1.compare(globalComm)==0 && comm1.compare(selfComm)!=0 &&
                 comm2.compare(globalComm)==0 && comm2.compare(selfComm)!=0 &&  
                 comm3.compare(globalComm)!=0 && comm3.compare(selfComm)==0 )
                ut.passes("intersection of globalComm and selfComm");
            else
                ut.failure("intersection of globalComm and selfComm");
        }
        // Test case where we have disjoint sets (this can only happen of one of the comms is null)
        AMP::AMP_MPI intersection = AMP::AMP_MPI::intersect( globalComm, nullComm );
        if ( intersection.isNull() )
            ut.passes("intersection of non-overlapping comms");
        else
            ut.failure("intersection of non-overlapping comms");
        // Test case where the comms partially overlap
        if ( globalComm.getSize() > 2 ) {
            int n = globalComm.getSize()-1;
            // Intersect 2 comms (all other ranks will be null)
            AMP::AMP_MPI split1 = globalComm.split(globalComm.getRank()==0?-1:0);
            AMP::AMP_MPI split2 = globalComm.split(globalComm.getRank()==n?-1:0);
            AMP::AMP_MPI intersection = AMP::AMP_MPI::intersect( split1, split2 );
            bool pass = true;
            if ( globalComm.getRank()==0 || globalComm.getRank()==n ) {
                if ( !intersection.isNull() )
                    pass = false;
            } else {
                if ( intersection.compare(split1)!=0 || intersection.compare(split2)!=0 ||
                     intersection.getSize()!=globalComm.getSize()-2 )
                    pass = false;
            }
            // Intersect 2 sets for ranks (3 groups should result)
            /*split1 = globalComm.split(globalComm.getRank()==0?1:2);
            split2 = globalComm.split(globalComm.getRank()==n?1:2);
            intersection = AMP::AMP_MPI::intersect( split1, split2 );
            bool pass = true;
            if ( globalComm.getRank()==0 || globalComm.getRank()==n ) {
                if ( intersection.compare(selfComm)==0 )
                    pass = false;
            } else {
                if ( intersection.compare(split1)!=0 || intersection.compare(split2)!=0 ||
                     intersection.getSize()!=globalComm.getSize()-2 )
                    pass = false;
            }*/
            if ( pass )
                ut.passes("intersection of partially overlapping comms");
            else
                ut.failure("intersection of partially overlapping comms");
        }

        // Test time and tick
        double end_time = AMP::AMP_MPI::time();
        double time_res = AMP::AMP_MPI::tick();
        if ( globalComm.getRank() == 0 ) {
            std::cout << "Time to run tests: " << end_time-start_time << std::endl;
            std::cout << "Timer resolution: " << time_res << std::endl;
            if ( time_res>0 && time_res<1 && (end_time-start_time)>=time_res )
                ut.passes("time and tick");
            else
                ut.failure("time and tick");
            std::cout << std::endl;
        }
        
        // Finished testing, report the results
        PROFILE_START("Report");
        start_time = AMP::AMP_MPI::time();
        ut.report();
        num_failed = ut.NumFailGlobal();
        end_time = AMP::AMP_MPI::time();
        if ( globalComm.getRank() == 0 )
            std::cout << "Time to report: " << end_time-start_time << std::endl << std::endl;
        PROFILE_STOP("Report");

        PROFILE_STOP("Main");
    } // Limit the scope so objects are detroyed

    // Shutdown
    PROFILE_SAVE("test_AMP_MPI");
    AMP::AMPManager::shutdown();
    return num_failed;
}   


