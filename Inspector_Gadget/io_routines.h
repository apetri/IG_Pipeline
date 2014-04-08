/*
 *  io_routines.h
 *  Inspector Gadget
 *
 *  Created by Jan Michael Kratochvil at Columbia University on 5/5/07.
 *  Copyright 2007. All rights reserved.
 *
 */

// Maximal allowable line length for read-in from input file:
#define MAXLINE 2000
// Maximal allowable path name length for paths from input file:
#define MAXPATHNAME 2000
// Maximal allowable name length for parameters and file names in input file:
#define MAXNAME 200
// The above values can be modified as needed for longer line and name lengths.
// WARNING: Make sure you also change the corresponding string array lengths in the struct analysis_parameters parameters in the file main.h!


int fgetline(FILE *stream, char s[], int lim);
void read_sampler(char filename[], int **readarray, int lines, int columns, int max_realizations, int subfield);
void read_analysis_parameter_file(FILE *stream, int file_type);

