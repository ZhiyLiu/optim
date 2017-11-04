

#ifndef PLAN_WINNT
#include <sys/time.h>
#endif

#ifdef PLAN_WINNT
#include <time.h>
#endif

#ifdef PLAN_AIX
#include <time.h>
#endif

#include <time.h>
/*#ifdef WIN32
#include <sys/time.h>
#else
#include <time.h>
#endif*/

#include <stdio.h>
#include "gen.h"
#include "libplan.h"



void
get_date(PDATE *pdate)
{
    time_t now;
    struct tm *l_time;

    time(&now);
    l_time = localtime(&now);
    pdate->day = l_time->tm_mday;
    pdate->month = l_time->tm_mon + 1;
    pdate->year = 2000 + l_time->tm_year % 100;
    pdate->dow = l_time->tm_wday;
    return;
}

void
get_time(PTIME *ptime)
{
    time_t now;
    struct tm *l_time;

    time(&now);
    l_time = localtime(&now);
    ptime->hour = l_time->tm_hour;
    ptime->minute = l_time->tm_min;
    ptime->second = l_time->tm_sec;
    return;
}


