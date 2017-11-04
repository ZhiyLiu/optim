/* srepFolder: The directory where s-reps stored.
 * For example: const char* srepFolder = "/NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/";
 * command like: ./srepProcrustes /NIRAL/work/ltu/WorkSpace/data_pablo/Pablo2/lib/vtksrep/Correspondence/models/alignment/
 * Liyun Tu
 * Mar 20, 2014
*/



#include "srepprocrustes.h"


int main( int argc, char* argv[] ){

    if( argc != 2 ) {
        std::cerr << "Usage: "<< std::endl;
        std::cerr << argv[0];
        std::cerr << " <srepFolder>";
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    const char * srepFolder = argv[1];

    srepprocrustes proc(srepFolder);

    proc.procrustesAlignment();


    return 0;
}

