
#include <stdio.h>
#include <unistd.h>
#include "plan_file_io.h"
#include "gen.h"
#include "extbm.h"
#include "libplanio.h"

#include <fcntl.h>
#include <fstream>

int
get_accessory(
    char	*phys_dat_dir,
    int		accessory_id,
    ACCESSORY	*accessory)
{
    int		fdes,
		loop;
    int		count;
    char	file_name[150];
    char	*malloc();

    sprintf(file_name, "%s/accessory", phys_dat_dir);
    fdes = open(file_name, O_RDONLY, 0);
    if (fdes < 0)
	return(-1);

    if (read_int(fdes, &count)) {
	close(fdes);
	return(-1);
    }

    for (loop = 0; loop < count; loop++) {
	if (read_accessory(fdes, accessory)) {
	    close(fdes);
	    return(-1);
	}

	if (accessory->accessory_id == accessory_id) {
	    close(fdes);
	    return(0);
	}
    }

    close(fdes);
    return(-1);
}

