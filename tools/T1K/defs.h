#ifndef _LSONG_DEFS_HEADER
#define _LSONG_DEFS_HEADER

#include <stdint.h>

//#define DEBUG

extern char nucToNum[26] ; 
extern char numToNuc[4] ;

#define MAX(x,y) (((x)>(y))?(x):(y))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define ABS(x) (((x)>(0))?(x):(-(x)))

typedef uint32_t index_t ; 

#define MAX_SEG_COUNT 127
struct _pair
{
	int a, b ;
} ;

struct _pair64
{
	int64_t a, b ;
} ;

struct _pairID
{
	int a ;
	double b ;
} ;

struct _pairIntDouble // integer, double
{
	int a ;
	double b ;
} ;

struct _triple
{
	int a, b, c ;
} ;

struct _tripleIID
{
	int a, b ;
	double c ;
} ;
/*struct _pair_b64
{
	int a ;
	uint64_t b ;
} ;*/

#endif
