/* Doing thin plate spline using multiple template.
 * The target s-reps are the avarage of those from different template.
 *
 * Liyun Tu
 * Aug 4, 2014
*/


#include "thinplatesplinepdmtosrep.h"



int main( int argc, char* argv[] ){

    if( argc != 5 ) {
        std::cerr << "Usage: "<< std::endl;
        std::cerr << argv[0];
        std::cerr << " <tamplatePDMFileList targetPDM tamplateSrepFileList targetSrepFolder>"<<std::endl;
        std::cerr << " tamplatePDMFileList: where the file listing all the filenames of tamplates surface PDM"<<std::endl;
        std::cerr << " targetPDM: where the target surface PDM stored"<<std::endl;
        std::cerr << " tamplateSrepFileList: where the file listing all the filename of tamplates' s-reps"<<std::endl;
        std::cerr << " targetSrepFolder: where the generated s-reps stored"<<std::endl;
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    thinplatesplinepdmtosrep obj;
    obj.multiple_template_tps(argv[1], argv[2], argv[3], argv[4]);

    return 0;
}


