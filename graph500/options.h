/* -*- mode: C; mode: folding; fill-column: 70; -*- */
/* Copyright 2010,  Georgia Institute of Technology, USA. */
/* See COPYING for license. */
#if !defined(OPTIONS_HEADER_)
#define OPTIONS_HEADER_

extern int VERBOSE;
extern int use_RMAT;
char *dumpname = "out.txt";
char *rootname = "root";

#define A_PARAM 0.57
#define B_PARAM 0.19
#define C_PARAM 0.19
/* Hence D = 0.05. */

double A = 0.57, B=0.19, C=0.19, D=0.05;

#define NBFS_max 64
int NBFS = 64;

#define default_SCALE ((int64_t)14)
#define default_edgefactor ((int64_t)16)

int64_t SCALE = default_SCALE;
int64_t edgefactor = default_edgefactor;

void get_options (int argc, char **argv);

#endif /* OPTIONS_HEADER_ */
