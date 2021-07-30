#pragma once

#ifdef MOSAIC_EXPORTS
#define MOSAIC_API __declspec(dllexport)
#else
#define MOSAIC_API __declspec(dllimport)
#endif
#define PI	(3.1415926535897932384626)
#define	TO_ANG(x)	(180.0*x/PI)
#define	TO_RAD(x)	(PI*x/180.0)

