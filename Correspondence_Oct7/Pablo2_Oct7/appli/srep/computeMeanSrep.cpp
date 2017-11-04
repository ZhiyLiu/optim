/* Given a list of sreps, this function compute the mean srep of them.
 * Liyun Tu
 * Aug 4, 2014
 *
 * ./computeMeanSrep /NIRAL/work/ltu/WorkSpace/Fitting-LV-sreps/Step_4_multi_template_TPS/test/srepFileList.txt /NIRAL/work/ltu/WorkSpace/Fitting-LV-sreps/Step_4_multi_template_TPS/test/N0147-2-01_segmented_vent_pp_surfSPHARM_scale_avarage.m3d

*/


#include "thinplatesplinepdmtosrep.h"



int main( int argc, char* argv[] ){

    if( argc != 3 ) {
        std::cerr << "Usage: "<< std::endl;
        std::cerr << argv[0];
        std::cerr << " <srepFileList outFileName>"<<std::endl;
        std::cerr << " srepFileList: where the file listing all the filenames of sreps"<<std::endl;
        std::cerr << " outFileName: the file name where the output mean stored"<<std::endl;
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    thinplatesplinepdmtosrep obj;

    // get srep pathname
    vector< std::string > srepPathNames = obj.readTempalte(argv[1]);


    //obj.avarage_srep(srepPathNames, argv[2]); //NOT CORRECT? THE SPOKE TIP SHOUDN'T BE USE TO COMPUTE THE AVARAGE.


    //compuate the avarage of the skeletal point, avarage of the spoke direction and avarage of the spoke radius.
    obj.avarage_srep_new(srepPathNames, argv[2]);


    obj.cpns_avarage_srep(srepPathNames, argv[2]);


    return 0;
}



