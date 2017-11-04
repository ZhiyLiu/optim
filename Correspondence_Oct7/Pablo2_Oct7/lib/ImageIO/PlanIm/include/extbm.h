/** @file extbm.h
    This header file contains data type and record declarations specific
    to external beam calculations.
*/

// HISTORY
// 29-july-2004 zsx: Add doxygen-style comments.

#ifndef _PLUNC_HAVE_EXTBM_H
#define _PLUNC_HAVE_EXTBM_H

/** Defines the maximum number of filters per beam. */
#define MAX_FILTER_COUNT	(2)
/** The maximum number of parts in a filter assembly, such as the
base plate and the wedge. */
#define MAX_CHUNK_COUNT		(4)

/** @name TMR table sizes */
//@{
/** Defines the maximum number of depths in a TMR table. */
#define MAX_SAR_DEPTHS		(40)
/** Defines the maximum number of radii in a TMR table. */
#define MAX_SAR_RADII		(30)
/** Defines the maximum number of flatness depths in a OAR table. */
#define MAX_FLATNESS_DEPTHS	(5)
/** Defines the maximum number of flatness radii in a OAR table. */
#define MAX_FLATNESS_RADII	(40)
/** Defines the maximum number of field sizes for output factors. */
#define MAX_Sc_POSITIONS	(25)
//@}

/*
  Some constants for dose engine
*/
/** @name Differential TMR table sizes */
//@{
/** Defines the maximum number of radii in a differential TMR table. */
#define dSAR_MAX_RADII		(30)
/** Defines the reasonable number of radii in a differential TMR
table. */
#define REASONABLE_RADIUS_COUNT	(15)
/** Defines the maximum number of sectors in a differential TMR table. */
#define dSAR_MAX_ANGLES		(60)
/** Defines the reasonable number of sectors in a differential TMR
table. */
#define REASONABLE_SECTOR_COUNT	(15)
//@}

/** @name Types of patient placement */
//@{
/** Type 1 of patient placement methods - SSD, which stands for
    source-to-surface distance. */
#define SSD_PLACEMENT		(1)
/** Type 2 of patient placement methods - isocentric. */
#define ISOCENTRIC		(2)
//@}

/** @name Types of inhomogeneity correction */
//@{
/** Type 1 inhomogeneity correction in dose calculation.
    Hyperlink to related reference here. */
#define BATHO			(1)
/** Type 2 inhomogeneity correction in dose calculation.
    Hyperlink to related reference here. */
#define EQ_TAR			(2)
//@}

/** @name Types of beam movement during treatment */
//@{
/** Type 1 of beam movement during treatment - dose is not delivered
    during beam movement. */
#define STATIC			(1)
/** Type 2 of beam movement during treatment - dose is delivered
    during beam movement. Not currently used. */
#define ARC			(2)
//@}

/** @name Types of treatment units */
//@{
/** Type 1 of units: photon external beam */
#define PHOTON			(1)
/** Type 2 of units: electron external beam */
#define ELECTRON		(2)
/** Type 3 of units: simulator */
#define SIMULATOR		(3)
/** Type 4 of units: ??? */
#define GEOM			(4)
//@}

/** @name Types of scatter tables */
//@{
/** Type 1 of scatter tables: TAR - Tissue-Air Ratio. */
#define TAR			(1)
/** Type 2 of scatter tables: TPR - Tissue-Phantom Ratio. */
#define TPR			(2)
/** Type 3 of scatter tables: TMR - Tissue-Maximum Ratio. */
#define TMR			(3)
//@}

/** This is not currently used.*/
#define S_THING(type) (((type) == TAR) ? "SAR" : \
		       (((type) == TMR) ? "SMR" : \
			(((type) == TPR) ? "SPR" : \
			 "S*R")))
/** This is not currently used.*/
#define T_THING(type) (((type) == TAR) ? "TAR" : \
		       (((type) == TMR) ? "TMR" : \
			(((type) == TPR) ? "TPR" : \
			 "T*R")))

/** @name Types of treatment times */
//@{
/** Type 1 of treatment times - monitor units. */
#define MU			(1)
/** Type 2 of treatment times - example, 3:15 (3 min and 15 sec) */
#define ENGLISH			(2)
/** Type 3 of treatment times - example, 3.25 (3.25 min) */
#define DECIMAL_UNIT		(3)
//@}

/** Macro to relate the name of a treatment time unit to the
    treatment time type. */
#define TIME_UNITS(units) (((units) == MU) ? "monitor unit(s)" : \
			   (((units) == ENGLISH) ? "minute(s)" : \
			    (((units) == DECIMAL_UNIT) ? "minute(s)" : \
			     "mystery unit(s)")))

/** Defines the maximum number of jaws (leaf pairs) allowed on a collimator. */
#define MAX_JAW_COUNT		(100)
/** Defines the maximum number of dose delivery segments for each
    beam.
*/
#define MAX_SEGMENT_COUNT	(100)

/** JAW_struct is a data structure defining the general
    parameters of the collimator jaw. */
typedef struct JAW_struct
{
/** Index of jaw movement type.
    0: left and right jaws move sysmetrically;\n
    1: left and right jaws move independently;\n
    2: MLC - multiple divergent jaws;\n
    3: MLC - multiple non-divergent jaws.
*/
    int independent;

/** @name Jaw Limits

    The central axis is 0.0; and the following parameters are relative
    to that.
    1 is for left or bottom jaw, and 2 is for right or top jaw.\n
    Negative values are to the left (or below) the central axis,
    positive values to the right (or above) the central axis.\n
    For example, a typical traditional accelerator may have values:\n
	independent = 0 \n
	min_1 = -20.0 \n
	max_1 = -2.0 \n
	min_2 = 2.0 \n
	max_2 = 20.0 \n
	min_w = -20.0 \n
	max_w = 20.0
*/
//@{
    float min_1;	/**< Leftmost or bottommost travel of the
			         left or bottom jaw. */
    float max_1;	/**< Rightmost or topmost travel of the
    			     left or bottom jaw. */
    float min_2;	/**< Leftmost or bottommost travel of the
    			     right or top jaw. */
    float max_2;	/**< Rightmost or topmost travel of the
    			     right or top jaw. */
    float min_w;	/**< Start of jaw in the width dimension. */
    float max_w;	/**< End of jaw in the width dimension. */
//@}

    float SDD;		/**< Primary source to jaw (diaphram) distance. */
    char  name_1[NAME_LENGTH]; /**< Name of the left or bottom jaw. */
    char  name_2[NAME_LENGTH]; /**< Name of the right or top jaw. */
} JAW;

/** JAW_POSITION_struct is a data structure for the position of a
    pair of jaws.
*/
typedef struct JAW_POSITION_struct
{
/** Index of jaw movement type.
    0: left and right (or top and bottom) jaws move sysmetrically;\n
    1: left and right (or top and bottom) jaws move independently;\n
    2: MLC - multiple leaf collimator, divergent leafs;\n
    3: MLC - multiple leaf collimator, non-divergent leafs.
*/
    int independent;
    float pos_1; //!< Position of the left/bottom jaw.\ 0.0 at c-axis.
    float pos_2; //!< Position of the right/top jaw.\ 0.0 at c-axis.
    float min_w; //!< Start of jaw orthogonal to jaw travel dir.
    float max_w; //!< End of jaw orthogonal to jaw travel dir.
} JAW_POSITION;

/** @name Types of jaws */
//@{
#define DOUBLE_DIVERGENT			0
#define DIVERGENT_SIDE_ROUNDED_END		1
#define NON_DIVERGENT_SIDE_DIVERGENT_END	2
#define NON_DIVERGENT_SIDE_ROUNDED_END		3
//@}

/** UNIT_struct is a data structure describing the characteristics of
    a treatment unit, which is defined as one modality at one energy
    level of a physical treatment machine.
*/
typedef struct UNIT_struct
{
    int version;        //!< Version of this unit data structure.
    int unit_id;        //!< Unique unit id number.
    char name[NAME_LENGTH]; //!< Name of this unit.
    int modality;       //!< PHOTON, ELECTRON, or SIMULATOR.
    float energy;	    //!< Energy in MV or Mev.
    float dmax;	        //!< The nominal depth of maximum dose.
    float SAD;          //!< Source to Axis Distance.
    float PSD;          //!< Primary src to Secondary src Distance.
    float source_diam1;	//!< Of primary, for penumbra calculations.
    float source_diam2;	//!< Of secondary, for penumbra calculations.
    float source_ratio;	//!< Fraction of source due to primary.
    float e_eq_diam;	//!< Electronic equilibrium blurring diameter.

    int x_jaw_count;    //!< Number of jaws (leaf pairs) in x-direction.
    int y_jaw_count;    //!< Number of jaws (leaf pairs) in y-direction.
    int x_jaw_type;     //!< DOUBLE_DIVERGENT, etc.
    int y_jaw_type;     //!< DOUBLE_DIVERGENT, etc.
    JAW jaw[MAX_JAW_COUNT]; //!< Array of jaws of this unit.

    /** @name Unit angles
    The gantry, table and collimator move in a unit-specific
	manner which is dependent on the manufacturer's markings, and
	these control the GUI and printouts that the user sees.
	See also the document on PLUNC's coordination systems. */
    //@{
    int gantry_sense;	//!< +1 => +angle=ccw; -1 => +angle=cw rotation
    float gantry_offset;  //!< Currently unused.
    float gantry_minimum; //!< Minimum limit of gantry rotation.
    float gantry_maximum; //!< Maximum limit of gantry rotation.

    int table_sense;	//!< +1 => +angle=ccw; -1 => +angle=cw rotation
    float table_offset;  //!< Currently unused.
    float table_minimum; //!< Minimum limit of table rotation.
    float table_maximum; //!< Maximum limit of table rotation.

    int collimator_sense; //!< +1 => +angle=ccw; -1 => +angle=cw rotation
    float collimator_offset;  //!< Currently unused.
    float collimator_minimum; //!< Minimum limit of collimator rotation.
    float collimator_maximum; //!< Maximum limit of collimator rotation.
    //@}
} UNIT;

/** TIME_CALC_PARAMETERS_struct is a data structure describing the
    time calculation parameters of a unit. */
typedef struct TIME_CALC_PARAMETERS_struct
{
    int version; //!< Version of this time_cal data structure.
    int unit_id; //!< Unit on which it can be used.
    float dose_rate; //!< Calibrated dose rate at cal_dis and cal_depth.
    PDATE calibration_date; //!< Date when calibratation was performed.
    float decay_constant; //!< In inverse hours, -1.0 if no decay.
    float cal_dist; //!< SSD + depth for calibration.
    float cal_depth; //!< Depth of calibration; negative if "in air".
    float end_effect; //!< actual time = set time + end_effect.
    int time_units;	//!< One of MU, ENGLISH, DECIMAL.
/** @name Output factor
  Output factor - this corrects for dose rate variation due strictly
  to changes in collimator size.

  The numbers below are used in a way consistent with Khan's Sc.  The
  value for a given field is looked for the rectangular collimator
  setting without regard to any beam shaping. \n

  The output factors, or Sc values if you prefer, are tabulated for
  your choice of rectangular fields. \n

  NOTE that the Sc_*_positions below are the collimator position.
  For instance, the x and y positions for a 10 by 10 would be
  5 and 5. \n
*/
//@{
    int Sc_x_count;	//!< Number of positions in x.
    int Sc_y_count;	//!< Number of positions in y.
    float Sc_x_positions[MAX_Sc_POSITIONS]; //!< Array of x positions.
    float Sc_y_positions[MAX_Sc_POSITIONS]; //!< Array of y positions.
    float Sc[MAX_Sc_POSITIONS][MAX_Sc_POSITIONS];	//!< Sc[y][x]
//@}
} TIME_CALC_PARAMETERS;

/** @name Types of filters */
//@{
#define REAL_FILTER		0
#define DYNAMIC_FILTER		1
#define BOLUS_FILTER		2
#define COMPENSATOR_FILTER	3
#define CUSTOM_FILTER_ID2	-3
//@}

/** @name Types of filters (Pre-version 6 style) */
//@{
#define CUSTOM_FILTER_ID	-1
#define ANASTRUCT_FILTER_ID	-2
//@}

/** FILTER_struct is a data structure for beam filters.
*/
typedef struct FILTER_struct
{
    int	version; //!< Version of this filter data structure.
    int unit_id; //!< ID of the unit to which the filter belongs.
    int type;    //!< Filter type constant.
    int filter_id;	//!< Unique filter id number.
/** @name Mirrorness of wedges

    Mirrorness is with respect to the y axis of the beam.  Conventional
    wedges which have their point toward or away from the gantry have
    themselves as mirrors.  Conventional wedges which have their point
    to the right or left usually have mirror images (often the same
    physical wedge inserted the other way around).  If a filter has a
    mirror wrt to the y axis, list the mirror's id here for use with
    copy and oppose. Otherwise, list NO_MIRROR.
*/
//@{
#define NO_MIRROR (0)
    int mirror_id;
    int next_id;
    int prev_id;
    char name[NAME_LENGTH];
#define NOT_POINTY (0)
#define NORTH (1)
#define SOUTH (2)
#define EAST (3)
#define WEST (4)
    int pointy_end;	/* Where is the pointy end? One of above. */
//@}

    JAW jaw_limits[2];	//!< Jaw limits for use of filter.
/** @name Transform matrix

    The filter is defined as a collection of prisms. The arbitrary
    decision was made that prisms are defined as polygons in x-y
    and translated along z. Hence for a wedge, z is perpendicular
    to the central axis. For a compensator, z is probably going to
    be parallel to the c. axis. That's what these two transforms
    are for - to move between the filter definition coordinate
    system and the beam (B) system. Don't forget that B is
    left-handed.  And don't forget to translate to the proper
    source distance.
*/
//@{
    float T_filter_to_beam[4][4];
    float T_beam_to_filter[4][4];
//@}

    float mu[MAX_CHUNK_COUNT];
    float mu_dx[MAX_CHUNK_COUNT];
    float mu_dr[MAX_CHUNK_COUNT];
    float mu_da[MAX_CHUNK_COUNT];
    float mu_dv[MAX_CHUNK_COUNT];
    float hvl_slope[MAX_CHUNK_COUNT];
    /*
      Store the linear attenuation coefficient (inverse cm, positive)
      for each prism in the prism's contours[foo].density.
      The min and max entries of the ANASTRUCT and CONTOURs are ignored
      except for the CONTOUR's min.z and max.z.
    */
/* @name Virtual wedge wedge-factor */
//@{
	float wfo_ao;
    float dwfo_da;
    float dwf_dfs_ao;
    float dwf_dfs_da;
//@}
    ANASTRUCT profiles; //!< The profile of the filter.
} FILTER;

/** @name Types of tray-hole shapes */
//@{
#define SQUARE_HOLE	(1)
#define ROUND_HOLE	(2)
#define SQUARE_SLOT	(3)
#define ROUND_SLOT	(4)
#define MLC_HOLE	(5)
//@}

/** TRAY_HOLE_struct is a data structure describing a hole in
    a tray that holds custom blocks.
*/
typedef struct TRAY_HOLE_struct
{
    int shape; //!< Shape of the hole: SQUARE_HOLE, ROUND_HOLE, etc.
    float x;   //!< x-position of the center of the hole.
    float y;   //!< y-position of the center of the hole.
    float height; //!< Height of the hole.
    float width;  //!< Width of the hole; Ignored for SQUARE and ROUND.
} TRAY_HOLE;

/** TRAY_struct is a data structure describing tray properties.
*/
typedef struct TRAY_struct
{
    int	version; //!< Version of this tray data structure.
    int unit_id; //!< Unit on which it can be used.
    int tray_id; //!< Unique tray id number.
    char name[NAME_LENGTH]; //!< Name of the tray.
    char code[NAME_LENGTH]; //!< Unique tray code.
    float tray_factor; //!< Dose attenuation factor of the tray.

/** @name tray sizes
	Dimmension of the tray, which can be asymmetric.
*/
//@{
    float xmin;
    float xmax;
    float ymin;
    float ymax;
//@}

    /**  Distance from the source of the beam to the side of tray
    where blocks are mounted on. */
    float tray_dist;

    /**  Distance from the source of the beam to the bottom of
    block for penumbra calc - different from tray_dist so block
    can hang from the tray. */
    float block_dist;

    int hole_count; //!< # of holes, which can be 0, but must be set.
    TRAY_HOLE *hole_list;	//!< Pointer to the list of holes. */
} TRAY;

/** ACCESSORY_struct is a data structure describing properites of
    an accessroy, such as a head holder or a breast board.
*/
typedef struct ACCESSORY_struct
{
    int	version; //!< Version of this accessroy data structure.
    int unit_id; //!< Unit on which it can be used.
    int accessory_id; //!< Unique accessory id number.
    int type; //!< Type of accessory.
    char name[NAME_LENGTH]; //!< Name of the tray.
    float accessory_factor; //!< Beam transmission factor.
} ACCESSORY;

/** SAR_TABLE_struct is a data structure describing the parameters of
    the tables defining the primary and scattering components of
    the dose distribution properties of a unit.
*/
typedef struct SAR_TABLE_struct
{
    int	version; //!< Version of this SAR-table data structure.
    int unit_id; //!< Unit on which it can be used.
    int type;    //!< Type of the table, one of TAR, TPR, TMR.

    float tran;	//!< Collimator transmission factor.
/** @name Table parameters and tables
*/
//@{
    int flatness_depth_count;
    int flatness_radius_count;
    float flatness_depth_vec[MAX_FLATNESS_DEPTHS];
    float flatness_radius_vec[MAX_FLATNESS_RADII];

    int depth_count;
    int radius_count;
    float depth_vec[MAX_SAR_DEPTHS];
    float radius_vec[MAX_SAR_RADII];
    float SAR[MAX_SAR_DEPTHS][MAX_SAR_RADII];

    int TAR0_depth_count;
    float TAR0_depth_vec[MAX_SAR_DEPTHS];
    float TAR0[MAX_SAR_DEPTHS];
//@}

/**
  One to MAX_FLATNESS_DEPTHS flatness profiles may be included here.
  The radius vector is at machine SAD and the flatness profiles are
  assumed to be depth scaled accordingly.
*/
    float flatness[MAX_FLATNESS_RADII][MAX_FLATNESS_DEPTHS];
/**
  This is a Phantom Scatter Factor (in principle the same as Khan's
  Sp) which is intended to correct for the change in scatter dose to
  the reference depth as a function of irradiated phantom area.  We do
  a Clarkson-type summation to determine the Sp for any given point.
  This table uses the same radii as the SAR table.  It is ignored for
  SAR/TAR calculations.
*/
    float Sp[MAX_SAR_RADII];
/** @name Sp0 and Sp_prime

  As delivered the Sp table is not suitable for
  the kind of Clarkson summation we intend to do with it,
  primarily because Sp(0) != 0.0  so what we do is to subtract
  the Sp(0) value from the table and thus build an integratable
  table stored here.  Then we just do a Clarkson-type integration
  and then add back the Sp0.
*/
//@{
    float Sp0;
    float Sp_prime[MAX_SAR_RADII];
//@}
} SAR_TABLE;

/** EXT_BEAM_struct is a data structure describing the properties of
    an external beam.
*/
typedef struct EXT_BEAM_struct
{
    unsigned long serial_number; /**< Unique instance id made with
                                      time stamp. */
    char name[NAME_LENGTH]; //!< Name of the beam.
    int version;  //!< Version of this beam data structure.
    int unit_id;  //!< Unit on which the beam is used.
    int x_jaw_count; //!< Number of jaws (leaf pairs) in x-direction.
    int y_jaw_count; //!< Number of jaws (leaf pairs) in y-direction.
    int segment_count; //!< # of segments to deliver the target dose.
    float *segment_weight; //!< Pointer to the segment weights.
    JAW_POSITION **jaw;	/**< Collimator setting: array of segment
                             pointers to JAW_POSITIONs  */
    PNT3D position;	/**<  Patient coordinates of central axis at
			      unit's definition distance. */
    PNT3D centroid;	//!< Center of rotation = position for ISOCENTRIC.
    float gantry_angle;	//!< "standard machine" reading for this beam.
    float coll_angle;	//!< "standard machine" reading for this beam.
    float table_angle;	//!< "standard machine" reading for this beam.

    float SSD; //!< Source-Surface-Distance given by virtual simulation.

    float T_pat_to_beam[4][4]; /**< Transformation matrix that takes
                               patient coordinates to beam coordinates.
                               */
    float T_beam_to_pat[4][4]; /**< The inverse matrix from beam
                               coordinates to patient coordinates. */

    int placement; //!< Patient placement: SSD_PLACEMENT or ISOCENTRIC.
    int extent;    //!< Beam movement type: STATIC or ARC.

/** Set LSB1 for dose inhomogeneity correction using CT numbers only.
*/
#define INH_CORR	(0x1)
/** Set LSB2 for dose inhomogeneity correction using CT numbers
    anywhere except the selected anastructs where given density
    numbers are used instead. */
#define CON_CORR	(0x2)
/** Bit OR of INH_CORR and CON_CORR to define dose inhomogneity
    correction requirement. */
    int corrections;
/** INH correction method: BATHO or EQ_TAR. */
    int inh_type;

#define NO_FILTER	(-1)
    int filter_id[MAX_FILTER_COUNT];
    int	filter_type[MAX_FILTER_COUNT];
    int vw_angle; //!< Virtual wedge angle.

#define NO_TRAY		(-1)
    int tray_id;	//!< NO_TRAY if id is absent.

#define NO_ACCESSORY	(-1)
    int accessory_id;	//!< NO_ACCESSORY if id is absent.

/** @name ARC

  The following three entries are of interest only for gantry-only
  arcs.  This section will get generalized later if the need arises to
  do compound-motion arcs.
*/
//@{
    float start_angle;
    float stop_angle;
    float angle_inc;
//@}

/**
  This is the beam outline at the SAD for the treatment unit. A
  vertex_count of 0 indicates no custom block. The z value of the
  CONTOUR must be set to the SAD for the treatment machine.
*/
    CONTOUR beam_outline;

/**
  This is the collimated version of the above outline.  It has a
  vertex_count > 2 and represents the intersection of the beam_outline
  with the collimator setting. If there is no beam_outline, then the
  collimated_beam_outline is determined by the collimator settings only.
*/
    CONTOUR collimated_beam_outline;

/**
  A custom filter (such as bolus or compensator) needs to have a set
  of profiles describing its shape, note, all general information
  about it is in the FILTER file.
*/

    ANASTRUCT	custom_filter[2];

} EXT_BEAM;

#endif

