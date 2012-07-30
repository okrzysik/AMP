#include "Utilities.h"

#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/Logger.h"
#include "utils/PIO.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdexcept>

// Detect the OS and include system dependent headers
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    // Note: windows has not been testeds
    #define USE_WINDOWS
    #include <windows.h>
    #include <stdio.h>   
    #include <tchar.h>
    #include <psapi.h>
    //#pragma comment(lib, psapi.lib) //added
    //#pragma comment(linker, /DEFAULTLIB:psapi.lib)
#elif defined(__APPLE__)
    #define USE_MAC
    #include <signal.h>
    #include <execinfo.h>
    #include <cxxabi.h>
    #include <dlfcn.h>
    #include<mach/mach.h>
#elif defined(__linux) || defined(__unix) || defined(__posix)
    #define USE_LINUX
    #include <signal.h>
    #include <execinfo.h>
    #include <cxxabi.h>
    #include <dlfcn.h>
#else
    #error Unknown OS
#endif


namespace AMP{

// Routine to check if a file exists
bool Utilities::fileExists( const std::string& filename )
{
   std::ifstream ifile(filename.c_str());
   return ifile!=NULL;
}

// Routine to rename a file.
void Utilities::renameFile(const std::string& old_filename, const std::string& new_filename)
{
   AMP_ASSERT(!old_filename.empty());
   AMP_ASSERT(!new_filename.empty());
   rename(old_filename.c_str(), new_filename.c_str());
}

// Routine to delete a file.
void Utilities::deleteFile(const std::string& filename)
{
   AMP_ASSERT(!filename.empty());
   if ( fileExists( filename ) ) {
      int error = remove(filename.c_str());
      AMP_INSIST(error==0,"Error deleting file");
   }
}

/*
 * Routine to recursively construct directories based on a relative path name.
 */
void Utilities::recursiveMkdir(
   const std::string& path, 
   mode_t mode,
   bool only_node_zero_creates)
{

#ifdef _MSC_VER
   const char seperator = '/';
   #define mkdir(path, mode) _mkdir(path)
#else
   const char seperator = '/';
#endif

   AMP_MPI comm = AMP_MPI(AMP_COMM_WORLD);
   if ( (!only_node_zero_creates) || (comm.getRank() == 0)) {
      int length = (int) path.length();
      char *path_buf= new char[length+1];
      sprintf(path_buf,"%s",path.c_str());
      struct stat status;
      int pos = length - 1;
   
      /* find part of path that has not yet been created */
      while ( (stat(path_buf,&status) != 0) && (pos >= 0) ) {
   
         /* slide backwards in string until next slash found */
         bool slash_found = false;
         while ( (!slash_found) && (pos >= 0) ) {
           if (path_buf[pos] == seperator) {
              slash_found = true;
              if (pos >= 0) path_buf[pos] = '\0';
           } else pos--;
         }
      } 

      /* 
       * if there is a part of the path that already exists make sure
       * it is really a directory
       */
      if (pos >= 0) {
         if ( !S_ISDIR(status.st_mode) ) {
            AMP_ERROR("Error in Utilities::recursiveMkdir...\n"
               << "    Cannot create directories in path = " << path
               << "\n    because some intermediate item in path exists and"
               << "is NOT a directory" << std::endl);
         }
      }
   
      /* make all directories that do not already exist */
   
      /* 
       * if (pos < 0), then there is no part of the path that
       * already exists.  Need to make the first part of the 
       * path before sliding along path_buf.
       */
      if (pos < 0) {
	 if(mkdir(path_buf,mode) != 0) {
	    AMP_ERROR("Error in Utilities::recursiveMkdir...\n"
		       << "    Cannot create directory  = " 
		       << path_buf << std::endl);
	 }
	 pos = 0;
      }
   
      /* make rest of directories */
      do {
   
         /* slide forward in string until next '\0' found */
         bool null_found = false;
         while ( (!null_found) && (pos < length) ) {
           if (path_buf[pos] == '\0') {
              null_found = true;
              path_buf[pos] = seperator;
           }
           pos++;
         }
   
         /* make directory if not at end of path */
	 if (pos < length) {
	    if(mkdir(path_buf,mode) != 0) {
	       AMP_ERROR("Error in Utilities::recursiveMkdir...\n"
			  << "    Cannot create directory  = " 
			  << path_buf << std::endl);
	    }
	 }
      } while (pos < length);

      delete [] path_buf;
   }

   /* 
    * Make sure all processors wait until node zero creates 
    * the directory structure.
    */
   if (only_node_zero_creates) {
     comm.barrier();
   }
}

/*
 * Routine to convert an integer to a string. 
 */
std::string Utilities::intToString(int num, int min_width)
{
   int tmp_width = ( min_width > 0 ? min_width : 1 );
   std::ostringstream os;
   if ( num < 0 ) {
      os << '-' << std::setw(tmp_width-1) << std::setfill('0') << -num;
   } else {
      os << std::setw(tmp_width) << std::setfill('0') << num;
   }
   os << std::flush;

   return(os.str()); //returns the string form of the stringstream object
}

std::string Utilities::nodeToString(int num) 
{
   return intToString(num, 5);
}

std::string Utilities::processorToString(int num) 
{
   return intToString(num, 5);
}


std::string Utilities::patchToString(int num) {
   return intToString(num, 4);
}


std::string Utilities::levelToString(int num) {
   return intToString(num, 4);
}


std::string Utilities::blockToString(int num) {
   return intToString(num, 4);
}


/*
 * Routine that calls abort and prints calling location to error stream.
 */
void Utilities::abort(const std::string &message, 
	              const std::string &filename, 
	              const int line) 
{
    if ( AMP::AMPManager::use_MPI_Abort==true) {
        // Print the call stack and memory usage
        long long unsigned int N_bytes = getMemoryUsage();
        printf("Bytes used = %llu\n",N_bytes);
        std::vector<std::string> stack = getCallStack();
        printf("Stack Trace:\n");
        for (size_t i=0; i<stack.size(); i++)
            printf("   %s",stack[i].c_str());
        printf("\n");
        // Log the abort message
        Logger::getInstance() -> logAbort(message, filename, line);
        // Use MPI_abort (will terminate all processes)
        AMP_MPI comm = AMP_MPI(AMP_COMM_WORLD);
        comm.abort();
    } else {
        // Throw and standard exception (allows the use of try, catch)
        // std::stringstream  stream;
        // stream << message << std::endl << "  " << filename << ":  " << line;
        // std::cout << stream.str() << std::endl;
        throw std::logic_error(message);
    }
}


// Function to create a 32-bit hash key from a character array
unsigned int Utilities::hash_char(const char* name)
{
    AMP_INSIST(sizeof(unsigned int)==4,"Need unsigned 32-bit int");
    unsigned int hash = 5381;
    unsigned char c;
    while((c = *name++)) {
        // hash = hash * 33 ^ c
        hash = ((hash << 5) + hash) ^ c;
    }
    return hash;
}


// Function to get the memory usage
size_t Utilities::getMemoryUsage()
{
    size_t N_bytes = 0;
    #if defined(USE_LINUX)
        std::string mem;
        std::ifstream proc("/proc/self/status");
        std::string s;
        while(getline(proc, s), !proc.fail()) {
            if(s.substr(0, 6) == "VmSize") {
                mem = s.substr(7);
                break;
            }
        }
        for (size_t i=0; i<mem.size(); i++) {
            if ( mem[i]==9 || mem[i]==10 || mem[i]==12 || mem[i]==13 )
                mem[i] = 32;
        }
        size_t mult = 1;
        if ( mem.find("kB")!=std::string::npos ) {
            mult = 0x400;
            mem.erase(mem.find("kB"),2);
        } else if ( mem.find("MB")!=std::string::npos ) {
            mult = 0x100000;
            mem.erase(mem.find("MB"),2);
        } else if ( mem.find("GB")!=std::string::npos ) {
            mult = 0x40000000;
            mem.erase(mem.find("GB"),2);
        }
        for (size_t i=0; i<mem.size(); i++) {
            if ( ( mem[i]<48 || mem[i]>57 ) && mem[i]!=32 ) {
                printf("Unable to get size from string: %s\n",s.c_str());
                return 0;
            }
        }
        long int N = atol(mem.c_str());
        N_bytes = ((size_t)N)*mult;
    #elif defined(USE_MAC)
        struct task_basic_info t_info;
        mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
        if (KERN_SUCCESS != task_info(mach_task_self(),
                              TASK_BASIC_INFO, (task_info_t)&t_info, 
                              &t_info_count)) {
            return 0;
        }
        N_bytes = t_info.resident_size;
    #elif defined(USE_WINDOWS)
        PROCESS_MEMORY_COUNTERS memCounter;
        GetProcessMemoryInfo( GetCurrentProcess(), &memCounter, sizeof(memCounter) );
        N_bytes = memCounter.WorkingSetSize;
    #endif
    return N_bytes;
}


//! Function to print the current call stack
std::vector<std::string>  Utilities::getCallStack()
{
    std::vector<std::string>  stack_list;
    #if defined(USE_LINUX) || defined(USE_MAC)
        void *trace[100];
        Dl_info dlinfo;
        int status;
        const char *symname;
        char *demangled=NULL;
        int trace_size = backtrace(trace,100);
        for (int i=0; i<trace_size; ++i) {  
            if(!dladdr(trace[i], &dlinfo))
                continue;
            symname = dlinfo.dli_sname;
            demangled = abi::__cxa_demangle(symname, NULL, 0, &status);
            if(status == 0 && demangled)
                symname = demangled;
            std::string object = std::string(dlinfo.dli_fname);
            std::string function = "";
            if ( symname!=NULL )
                function = std::string(symname);
            if ( i!=0 ) {  // Skip the current function
                std::string stack_item = object + ":   " + function + "\n";
                //stack_item = "object: " + object + "\n";
                //stack_item += "function: " + function + "\n";
                stack_list.push_back(stack_item);
            }
            if ( demangled==NULL ) {
                free(demangled);
                demangled=NULL;
            }
        } 
    #elif defined(USE_WINDOWS)


    #endif
    return stack_list;
}


// Print AMP Banner
void Utilities::printBanner()
{
    std::ostringstream banner;
    banner << std::endl;
    banner << "            _____                    _____                    _____"           << std::endl;
    banner << "           /\\    \\                  /\\    \\                  /\\    \\ "         << std::endl;
    banner << "          /::\\    \\                /::\\____\\                /::\\    \\"         << std::endl;        
    banner << "         /::::\\    \\              /::::|   |               /::::\\    \\"        << std::endl;       
    banner << "        /::::::\\    \\            /:::::|   |              /::::::\\    \\"       << std::endl;      
    banner << "       /:::/\\:::\\    \\          /::::::|   |             /:::/\\:::\\    \\"      << std::endl;     
    banner << "      /:::/__\\:::\\    \\        /:::/|::|   |            /:::/__\\:::\\    \\"     << std::endl;    
    banner << "     /::::\\   \\:::\\    \\      /:::/ |::|   |           /::::\\   \\:::\\    \\"    << std::endl;   
    banner << "    /::::::\\   \\:::\\    \\    /:::/  |::|___|______    /::::::\\   \\:::\\    \\"   << std::endl;  
    banner << "   /:::/\\:::\\   \\:::\\    \\  /:::/   |::::::::\\    \\  /:::/\\:::\\   \\:::\\____\\"  << std::endl; 
    banner << "  /:::/  \\:::\\   \\:::\\____\\/:::/    |:::::::::\\____\\/:::/  \\:::\\   \\:::|    |" << std::endl;
    banner << "  \\::/    \\:::\\  /:::/    /\\::/    / ~~~~~/:::/    /\\::/    \\:::\\  /:::|____|" << std::endl;
    banner << "   \\/____/ \\:::\\/:::/    /  \\/____/      /:::/    /  \\/_____/\\:::\\/:::/    /"  << std::endl; 
    banner << "            \\::::::/    /               /:::/    /            \\::::::/    /"   << std::endl;  
    banner << "             \\::::/    /               /:::/    /              \\::::/    /"    << std::endl;   
    banner << "             /:::/    /               /:::/    /                \\::/____/"     << std::endl;    
    banner << "            /:::/    /               /:::/    /"                               << std::endl;                            
    banner << "           /:::/    /               /:::/    /"                                << std::endl;                               
    banner << "          /:::/    /               /:::/    /"                                 << std::endl;                                
    banner << "          \\::/    /                \\::/    /"                                  << std::endl;                                 
    banner << "           \\/____/                  \\/____/"                                   << std::endl;                                  
    banner << std::endl << std::endl;

    AMP::pout << banner.str();

}

// Factor a number into it's prime factors
std::vector<int> Utilities::factor(size_t number)
{
    if ( number<=3 ) 
        return std::vector<int>(1,(int)number);
    size_t i, n, n_max;
    bool factor_found;
    // Compute the maximum number of factors
    int N_primes_max = 1;
    n = number;
    while (n >>= 1) ++N_primes_max;
    // Initialize n, factors 
    n = number;
    std::vector<int> factors;
    factors.reserve(N_primes_max);
    while ( 1 ) {
        // Check if n is a trivial prime number
        if ( n==2 || n==3 || n==5 ) {
            factors.push_back( (int) n );
            break;
        } 
        // Check if n is divisible by 2
        if ( n%2 == 0 ) {
            factors.push_back( 2 );
            n/=2;
            continue;
        } 
        // Check each odd number until a factor is reached
        n_max = (size_t) floor(sqrt((double) n));
        factor_found = false;
        for (i=3; i<=n_max; i+=2) {
            if ( n%i == 0 ) {
                factors.push_back( i );
                n/=i;
                factor_found = true;
                break;
            } 
        }
        if ( factor_found )
            continue;
        // No factors were found, the number must be prime
        factors.push_back( (int) n );
        break;
    }
    // Sort the factors
    AMP::Utilities::quicksort(factors);
    return factors;
}



}

