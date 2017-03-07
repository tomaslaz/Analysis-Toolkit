/*******************************************************************************
 ** Copyright Tomas Lazauskas 2017
 *******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "utilities.h"
#include "boxeslib.h"

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

/*******************************************************************************
 * Search for defects and return the sub-system surrounding them
 *******************************************************************************/
int findDefects( int includeVacs, int includeInts, int includeAnts,
				 int* defectList, int* NDefectsByType,
				 int* vacancies, int* interstitials, int* antisites, int* onAntisites,
				 int inclSpecDim, int* inclSpec, int exclSpecInputDim, int* exclSpecInput, int exclSpecRefDim, int* exclSpecRef,
				 int NAtoms, char* specieList, int* specie, double* pos,
				 int refNAtoms, char* specieListRef, int* specieRef, double* refPos,
				 double *cellDims, int *PBC, double vacancyRadius, double inclusionRadius, double *minPos, double *maxPos,
                 int verboseLevel, int debugDefects)
{
    int i, exitLoop, k, j, index;
    double vacRad2, xpos, ypos, zpos;
    int boxNebList[27];
    char symtemp[3], symtemp2[3];
    int checkBox, refIndex, comp, boxIndex;
    double refxpos, refypos, refzpos;
    double sep2, incRad2, approxBoxWidth;
    int NDefects, NAntisites, NInterstitials, NVacancies;
    int *possibleVacancy, *possibleInterstitial;
    int *possibleAntisite, *possibleOnAntisite;
    int count, addToInt, skip;
    struct Boxes *boxes;

    if ( verboseLevel > 2)
    {
        printf("CLIB: finding defects\n");
        printf("  vac rad %f\n", vacancyRadius);
        printf("  inc rad %f\n", inclusionRadius);
    }

    /* approx width, must be at least vacRad */
    approxBoxWidth = 1.1 * vacancyRadius;

    /* box reference atoms */
    boxes = setupBoxes(approxBoxWidth, minPos, maxPos, PBC, cellDims);
    putAtomsInBoxes(refNAtoms, refPos, boxes);

    /* allocate local arrays for checking atoms */
    possibleVacancy = malloc( refNAtoms * sizeof(int) );
    if (possibleVacancy == NULL)
    {
        printf("ERROR: Boxes: could not allocate possibleVacancy\n");
        exit(1);
    }

    possibleInterstitial = malloc( NAtoms * sizeof(int) );
    if (possibleInterstitial == NULL)
    {
        printf("ERROR: Boxes: could not allocate possibleInterstitial\n");
        exit(1);
    }

    possibleAntisite = malloc( refNAtoms * sizeof(int) );
    if (possibleAntisite == NULL)
    {
        printf("ERROR: Boxes: could not allocate possibleAntisite\n");
        exit(1);
    }

    possibleOnAntisite = malloc( refNAtoms * sizeof(int) );
    if (possibleOnAntisite == NULL)
    {
        printf("ERROR: Boxes: could not allocate possibleOnAntisite\n");
        exit(1);
    }

    /* initialise arrays */
    for ( i=0; i<NAtoms; i++ )
    {
        possibleInterstitial[i] = 1;
    }
    for ( i=0; i<refNAtoms; i++ )
    {
        possibleVacancy[i] = 1;
        possibleAntisite[i] = 1;
    }

    /* build local specie list */

    vacRad2 = vacancyRadius * vacancyRadius;
    incRad2 = inclusionRadius * inclusionRadius;

    /* loop over all input atoms */
    for ( i=0; i<NAtoms; i++ )
    {
        int nboxes;

        xpos = pos[3*i];
        ypos = pos[3*i+1];
        zpos = pos[3*i+2];

        /* get box index of this atom */
        boxIndex = boxIndexOfAtom(xpos, ypos, zpos, boxes);

        /* find neighbouring boxes */
        nboxes = getBoxNeighbourhood(boxIndex, boxNebList, boxes);

        /* loop over neighbouring boxes */
        exitLoop = 0;
        for (j = 0; j < nboxes; j++)
        {
            if (exitLoop)
            {
                break;
            }

            checkBox = boxNebList[j];

            /* now loop over all reference atoms in the box */
            for ( k=0; k<boxes->boxNAtoms[checkBox]; k++ )
            {
                refIndex = boxes->boxAtoms[checkBox][k];

                /* if this vacancy has already been filled then skip to the next one */
                if ( possibleVacancy[refIndex] == 0 )
                {
                    continue;
                }

                refxpos = refPos[3*refIndex];
                refypos = refPos[3*refIndex+1];
                refzpos = refPos[3*refIndex+2];

                /* atomic separation of possible vacancy and possible interstitial */
                sep2 = atomicSeparation2(xpos, ypos, zpos, refxpos, refypos, refzpos,
                                         cellDims[0], cellDims[1], cellDims[2], PBC[0], PBC[1], PBC[2]);

                /* if within vacancy radius, is it an antisite or normal lattice point */
                if ( sep2 < vacRad2 )
                {
                    /* compare symbols */
                    symtemp[0] = specieList[2*specie[i]];
                    symtemp[1] = specieList[2*specie[i]+1];
                    symtemp[2] = '\0';

                    symtemp2[0] = specieListRef[2*specieRef[refIndex]];
                    symtemp2[1] = specieListRef[2*specieRef[refIndex]+1];
                    symtemp2[2] = '\0';

                    comp = strcmp( symtemp, symtemp2 );
                    if ( comp == 0 )
                    {
                        /* match, so not antisite */
                        possibleAntisite[refIndex] = 0;
                    }
                    else
                    {
                        possibleOnAntisite[refIndex] = i;
                    }

                    /* not an interstitial or vacancy */
                    possibleInterstitial[i] = 0;
                    possibleVacancy[refIndex] = 0;

                    /* no need to check further for this (no longer) possible interstitial */
                    exitLoop = 1;
                    break;
                }
            }
        }
    }

    /* free boxes structure */
    freeBoxes(boxes);

    /* now box input atoms, approx width must be at least incRad */
    approxBoxWidth = 1.1 * inclusionRadius;
    boxes = setupBoxes(approxBoxWidth, minPos, maxPos, PBC, cellDims);
    putAtomsInBoxes(NAtoms, pos, boxes);

    /* now classify defects */
    NVacancies = 0;
    NInterstitials = 0;
    NAntisites = 0;
    count = 0;
    for ( i=0; i<refNAtoms; i++ )
    {
        skip = 0;
        for (j=0; j<exclSpecRefDim; j++)
        {
            if (specieRef[i] == exclSpecRef[j])
            {
                skip = 1;
            }
        }

        if (skip == 1)
        {
            continue;
        }

        /* vacancies */
        if (possibleVacancy[i] == 1)
        {
            if (includeVacs == 1)
            {
                int nboxes;

                vacancies[NVacancies] = i;
                NVacancies++;

                /* find input atoms within inclusionRadius of this vacancy */
                refxpos = refPos[3*i];
                refypos = refPos[3*i+1];
                refzpos = refPos[3*i+2];
                boxIndex = boxIndexOfAtom(refxpos, refypos, refzpos, boxes);

                /* find neighbouring boxes */
                nboxes = getBoxNeighbourhood(boxIndex, boxNebList, boxes);

                for (j = 0; j < nboxes; j++ )
                {
                    checkBox = boxNebList[j];

                    /* loop over atoms in box */
                    for ( k=0; k<boxes->boxNAtoms[checkBox]; k++ )
                    {
                        index = boxes->boxAtoms[checkBox][k];

                        /* if already on defect list continue */
                        if (defectList[index] == 1)
                        {
                            continue;
                        }

                        /* if close to defect add to list */
                        xpos = pos[3*index];
                        ypos = pos[3*index+1];
                        zpos = pos[3*index+2];
                        sep2 = atomicSeparation2(xpos, ypos, zpos, refxpos, refypos, refzpos,
                                                 cellDims[0], cellDims[1], cellDims[2], PBC[0], PBC[1], PBC[2]);

                        if ( sep2 < incRad2 )
                        {
                            defectList[index] = 1;
                            count++;
                        }
                    }
                }
            }
        }

        /* antisites */
        else if ( (possibleAntisite[i] == 1) && (includeAnts == 1) )
        {
            int nboxes;

            antisites[NAntisites] = i;
            onAntisites[NAntisites] = possibleOnAntisite[i];
            NAntisites++;

            /* find input atoms within inclusionRadius of this antisite */
            refxpos = refPos[3*i];
            refypos = refPos[3*i+1];
            refzpos = refPos[3*i+2];
            boxIndex = boxIndexOfAtom(refxpos, refypos, refzpos, boxes);

            /* find neighbouring boxes */
            nboxes = getBoxNeighbourhood(boxIndex, boxNebList, boxes);

            for (j = 0; j < nboxes; j++)
            {
                checkBox = boxNebList[j];

                /* loop over atoms in box */
                for ( k=0; k<boxes->boxNAtoms[checkBox]; k++ )
                {
                    index = boxes->boxAtoms[checkBox][k];

                    /* if already on defect list continue */
                    if (defectList[index] == 1)
                    {
                        continue;
                    }

                    /* if close to defect add to list */
                    xpos = pos[3*index];
                    ypos = pos[3*index+1];
                    zpos = pos[3*index+2];
                    sep2 = atomicSeparation2(xpos, ypos, zpos, refxpos, refypos, refzpos,
                                             cellDims[0], cellDims[1], cellDims[2], PBC[0], PBC[1], PBC[2]);

                    if ( sep2 < incRad2 )
                    {
                        defectList[index] = 1;
                        count++;
                    }
                }
            }
        }
    }

    if (includeInts == 1)
    {
        for ( i=0; i<NAtoms; i++ )
        {
            skip = 0;
            for (j=0; j<exclSpecInputDim; j++)
            {
                if (specie[i] == exclSpecInput[j])
                {
                    skip = 1;
                    break;
                }
            }

            if (skip == 1)
            {
                continue;
            }

            addToInt = 0;
            for (j=0; j<inclSpecDim; j++)
            {
                if (specie[i] == inclSpec[j])
                {
                    addToInt = 1;
                    break;
                }
            }

            /* interstitials */
            if ( (possibleInterstitial[i] == 1) || (addToInt == 1) )
            {
                int nboxes;

                interstitials[NInterstitials] = i;
                NInterstitials++;
                defectList[i] = 1;
                count++;

                /* find input atoms within inclusionRadius of this interstitial */
                refxpos = pos[3*i];
                refypos = pos[3*i+1];
                refzpos = pos[3*i+2];
                boxIndex = boxIndexOfAtom(refxpos, refypos, refzpos, boxes);

                /* find neighbouring boxes */
                nboxes = getBoxNeighbourhood(boxIndex, boxNebList, boxes);

                for (j = 0; j < nboxes; j++)
                {
                    checkBox = boxNebList[j];

                    /* loop over atoms in box */
                    for ( k=0; k<boxes->boxNAtoms[checkBox]; k++ )
                    {
                        index = boxes->boxAtoms[checkBox][k];

                        /* if same atom continue */
                        if ( index == i )
                        {
                            continue;
                        }

                        /* if already on defect list continue */
                        if (defectList[index] == 1)
                        {
                            continue;
                        }

                        /* if close to defect add to list */
                        xpos = pos[3*index];
                        ypos = pos[3*index+1];
                        zpos = pos[3*index+2];
                        sep2 = atomicSeparation2(xpos, ypos, zpos, refxpos, refypos, refzpos,
                                                 cellDims[0], cellDims[1], cellDims[2], PBC[0], PBC[1], PBC[2]);

                        if ( sep2 < incRad2 )
                        {
                            defectList[index] = 1;
                            count++;
                        }
                    }
                }
            }
        }
    }

    /* shift indexes to zero */
    count = 0;
    for (i=0; i<NAtoms; i++)
    {
        if (defectList[i] == 1)
        {
            defectList[count] = i;
            count++;
        }
    }

    NDefects = NVacancies + NInterstitials + NAntisites;
    if ( verboseLevel > 2 )
    {
        printf("  found %d defects\n", NDefects);
        printf("    %d vacancies\n", NVacancies);
        printf("    %d interstitials\n", NInterstitials);
        printf("    %d antisites\n", NAntisites);
        printf("  including %d atoms in sub-system\n", count);
    }

    NDefectsByType[0] = NDefects;
    NDefectsByType[1] = NVacancies;
    NDefectsByType[2] = NInterstitials;
    NDefectsByType[3] = NAntisites;

    freeBoxes(boxes);
    free(possibleVacancy);
    free(possibleInterstitial);
    free(possibleAntisite);
    free(possibleOnAntisite);

    return count;
}
