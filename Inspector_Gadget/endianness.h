/*                                                                                                                                                                                       
 *  endianness.h                                                                                                                                                                         
 *  Inspector Gadget                                                                                                                                                                           
 *                                                                                                                                                                                       
 *  Created by Jan Kratochvil at Columbia University on 9/19/09.                                                                                                                                                
 *  Copyright 2009. All rights reserved.                                                                                                                               
 *                                                                                                                                                                                       
 */

int machineEndianness();

/*
 inline short ShortSwap( short s );
 inline int IntSwap( int n );
 inline float FloatSwap( float f );
 inline double DoubleSwap( double f );
 */

struct endian_parameters
{
	int byteswap;
};

extern struct endian_parameters endianness;


short ShortSwap( short s );
int IntSwap( int n );
long long LongSwap (long long i);
float FloatSwap( float f );
double DoubleSwap( double f );

void singleint_fread(int *x, FILE *fd);
void int_fread(int *x, size_t size, size_t number, FILE *fd);
void singlelong_fread(long long *x, FILE *fd);
void long_fread(long long *x, size_t size, size_t number, FILE *fd);
void singlefloat_fread(float *x, FILE *fd);
void float_fread(float *x, size_t size, size_t number, FILE *fd);
void singledouble_fread(double *x, FILE *fd);
void double_fread(double *x, size_t size, size_t number, FILE *fd);

void header_fread(struct io_header_1 *x, size_t size, size_t number, FILE *fd);
int snapshotEndianness(char *fname, int file_number);

/*
 inline void singleint_fread(int *x, FILE *fd);
 inline void int_fread(int *x, int size, int number, FILE *fd);
 inline void singlefloat_fread(float *x, FILE *fd);
 inline void float_fread(float *x, int size, int number, FILE *fd);
 inline void singledouble_fread(double *x, FILE *fd);
 inline void double_fread(double *x, int size, int number, FILE *fd);
 inline void header_fread(struct io_header_1 *x, int size, int number, FILE *fd);
 */


