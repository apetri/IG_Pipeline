/*
 *  fits.h
 *  FITS-Modifier
 *
 *  Created by Jan Michael Kratochvil on 1/26/09 at Columbia University.
 *  Copyright 2009 Jan Michael Kratochvil. All rights reserved.
 *
 */


void readFITSimage_f(char filename[], double *imagearray);
void writeFITSimage_f(char filename[], long naxis, long naxes[], double *imagearray);

void printerror(int status);

