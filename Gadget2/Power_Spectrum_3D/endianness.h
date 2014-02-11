/*
 *  endianness.h
 *  Gadget-2 Snapshot File Reader (low memory version with buffer - 11/03/2012).
 *
 *  Created by Jan Michael Kratochvil on 9/19/2009 at Columbia University.
 *  Copyright 2009 Jan Michael Kratochvil. All rights reserved.
 *
 */


int machineEndianness();
//int snapshotEndianness(char *fname, int file_number, int snapshot_number);
int snapshotEndianness(char *fname);

inline short ShortSwap( short s );
inline int IntSwap( int n );
inline float FloatSwap( float f );
inline double DoubleSwap( double f );

inline void singleint_fread(int *x, FILE *fd);
inline void int_fread(int *x, int size, int number, FILE *fd);
inline void singlefloat_fread(float *x, FILE *fd);
inline void float_fread(float *x, int size, int number, FILE *fd);
inline void singledouble_fread(double *x, FILE *fd);
inline void double_fread(double *x, int size, int number, FILE *fd);
inline void header_fread(struct io_header_1 *x, int size, int number, FILE *fd);

