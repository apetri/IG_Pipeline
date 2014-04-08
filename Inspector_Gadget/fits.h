/*
 *  fits.h
 *  Inspector Gadget
 *
 *  Created by Jan Kratochvil at Columbia University on 7/11/07.
 *  Copyright 2007. All rights reserved.
 *
 */

void readFITSheader (char filename[], int plane_number, struct fitsheader *FITSheader, struct plane_2D *Plane);
void readFITSpotential_singleplane(char filename[], double *potential_array);
void writePlaneFITSimage_f(char filename[], long naxis, long naxes[], double *imagearray, int plane_number);
void writeRayPlaneFITSimage_f(char filename[], long naxis, long naxes[], double *imagearray, int plane_number);
void writeWLmapFITSimage_f(char filename[], long naxis, long naxes[], double *imagearray,int plane_number);
void printerror(int status);

