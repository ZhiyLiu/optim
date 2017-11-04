/* srepFolder: The directory where s-reps stored.
 * command: ./tps_srep template.vtk test_target.vtk template.m3d ../target_srep_folder/
 * Liyun Tu
 * Mar 20, 2014
*/



//#include "thinplatesplinesrep.h"

#include "thinplatesplinepdmtosrep.h"



int main( int argc, char* argv[] ){

    if( argc != 5 ) {
        std::cerr << "Usage: "<< std::endl;
        std::cerr << argv[0];
        std::cerr << " <tamplate PDM, target PDM, tamplate s-rep and target s-rep folders>";
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    const char* srepFolder = argv[1];

    //string sourcefilename = string(srepFolder) + "/template_PP_SPHARMSPHARM.vtk";
    //string targetfilename = string(srepFolder) + "/N0113-1-01_segmented_vent_pp_surfSPHARM_scale.vtk";
    //string sourcesrepfilename = string(srepFolder) + "/template.m3d";

    // Method one
    //thinplatesplinesrep obj;
    //obj.tps_srep(sourcefilename.c_str(), targetfilename.c_str(), sourcesrepfilename.c_str());


    // Method two
    thinplatesplinepdmtosrep obj;
    //obj.tps_to_srep(sourcefilename.c_str(), targetfilename.c_str(), sourcesrepfilename.c_str());
    obj.tps_to_srep(argv[1], argv[2], argv[3], argv[4]);


    return 0;
}


