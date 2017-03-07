
/*******************************************************************************
 **
 *******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*******************************************************************************
 ** return atomic separation squared
 *******************************************************************************/
double atomicSeparation2(double ax, double ay, double az, double bx, double by, double bz, double xdim, double ydim, double zdim, int pbcx, int pbcy, int pbcz)
{
    double rx, ry, rz, r2;
    
    /* calculate separation */
    rx = ax - bx;
    ry = ay - by;
    rz = az - bz;
    
    /* handle PBCs here if required */
    if ( pbcx == 1 )
    {
        rx = rx - round( rx / xdim ) * xdim;
    }
    if ( pbcy == 1 )
    {
        ry = ry - round( ry / ydim ) * ydim;
    }
    if ( pbcz == 1 )
    {
        rz = rz - round( rz / zdim ) * zdim;
    }
    
    /* separation squared */
    r2 = rx * rx + ry * ry + rz * rz;
    
    return r2;
}



