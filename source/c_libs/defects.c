/*******************************************************************************
 ** Copyright Tomas Lazauskas 2017
 *******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include "defects.h"

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



/*******************************************************************************
 ** create and return pointer to Boxes structure
 ** #TODO: if PBCs are set min/max pos should be equal to cell dims
 *******************************************************************************/
struct Boxes * setupBoxes(double approxBoxWidth, double *minPos, double *maxPos, int *PBC, double *cellDims)
{
    int i;
    int maxBoxDim = cbrt((double) INT_MAX); // avoid overflow
    struct Boxes *boxes;


    /* allocate space for boxes struct */
    boxes = malloc(1 * sizeof(struct Boxes));
    if (boxes == NULL)
    {
        printf("ERROR: could not allocate boxes\n");
        exit(50);
    }

    /* setup boxes */
    for (i = 0; i < 3; i++)
    {
        double cellLength;

        /* store some parameters */
        boxes->PBC[i] = PBC[i];
        boxes->cellDims[i] = cellDims[i];

        /* if PBC box cell, otherwise box min-max pos */
        if (boxes->PBC[i] == 0)
        {
            boxes->minPos[i] = minPos[i];
            boxes->maxPos[i] = maxPos[i];
        }
        else
        {
            boxes->minPos[i] = 0.0;
            boxes->maxPos[i] = boxes->cellDims[i];
        }

        /* size of the region in this direction */
        cellLength = boxes->maxPos[i] - boxes->minPos[i];
        cellLength = (cellLength < 1.0) ? 1.0 : cellLength;

        /* number of boxes in this direction */
        if (approxBoxWidth * maxBoxDim <= cellLength)
            /* set to max num boxes */
            boxes->NBoxes[i] = maxBoxDim;
        else
        {
            boxes->NBoxes[i] = (int) (cellLength / approxBoxWidth);
            boxes->NBoxes[i] = (boxes->NBoxes[i] == 0) ? 1 : boxes->NBoxes[i];
        }

        /* length of box side */
        boxes->boxWidth[i] = cellLength / boxes->NBoxes[i];
    }
    /* total number of boxes */
    boxes->totNBoxes = boxes->NBoxes[0] * boxes->NBoxes[1] * boxes->NBoxes[2];

    /* amount of memory to allocate at a time */
    boxes->allocChunk = 16;

    /* allocate arrays for storing counters and atoms */
    boxes->boxNAtoms = calloc(boxes->totNBoxes, sizeof(int));
    if (boxes->boxNAtoms == NULL)
    {
        printf("ERROR: could not allocate boxNAtoms\n");
        exit(50);
    }
    boxes->boxAtoms = malloc(boxes->totNBoxes * sizeof(int *));
    if (boxes->boxAtoms == NULL)
    {
        printf("ERROR: could not allocate boxAtoms\n");
        exit(50);
    }

    return boxes;
}


/*******************************************************************************
 ** put atoms into boxes
 *******************************************************************************/
void putAtomsInBoxes(int NAtoms, double *pos, struct Boxes *boxes)
{
    int i;


    for (i = 0; i < NAtoms; i++)
    {
        int boxIndex = boxIndexOfAtom(pos[3*i], pos[3*i+1], pos[3*i+2], boxes);

        /* check if this box is empty or full */
        if (boxes->boxNAtoms[boxIndex] == 0)
        {
            /* allocate this box */
            boxes->boxAtoms[boxIndex] = malloc(boxes->allocChunk * sizeof(int));
            if (boxes->boxAtoms[boxIndex] == NULL)
            {
                printf("ERROR: could not allocate boxAtoms[%d]\n", boxIndex);
                exit(50);
            }
        }

        else if (boxes->boxNAtoms[boxIndex] % boxes->allocChunk == 0)
        {
            int newsize;

            /* realloc more space */
            newsize = boxes->boxNAtoms[boxIndex] + boxes->allocChunk;
            boxes->boxAtoms[boxIndex] = realloc(boxes->boxAtoms[boxIndex], newsize * sizeof(int));
            if (boxes->boxAtoms[boxIndex] == NULL)
            {
                printf("ERROR: could not reallocate boxAtoms[%d]: %d\n", boxIndex, newsize);
                exit(50);
            }
        }

        /* add atom to box */
        boxes->boxAtoms[boxIndex][boxes->boxNAtoms[boxIndex]] = i;
        boxes->boxNAtoms[boxIndex]++;
    }
}


/*******************************************************************************
 ** returns box index of given atom
 *******************************************************************************/
int boxIndexOfAtom(double xpos, double ypos, double zpos, struct Boxes *boxes)
{
    int posintx, posinty, posintz;
    int boxIndex;


    /* if atom is outside boxes min/max pos we wrap or translate it back depending on PBCs */
    if (xpos > boxes->maxPos[0] || xpos < boxes->minPos[0])
    {
        if (boxes->PBC[0] == 1)
        {
            /* wrap position */
            xpos = xpos - floor( xpos / boxes->cellDims[0] ) * boxes->cellDims[0];
        }
        else
        {
            if (xpos > boxes->maxPos[0])
            {
                /* put it in the end box */
                xpos = boxes->maxPos[0] - 0.5 * boxes->boxWidth[0];
            }
            else
            {
                /* put it in the end box */
                xpos = boxes->minPos[0] + 0.5 * boxes->boxWidth[0];
            }
        }
    }

    if (ypos > boxes->maxPos[1] || ypos < boxes->minPos[1])
    {
        if (boxes->PBC[1] == 1)
        {
            /* wrap position */
            ypos = ypos - floor( ypos / boxes->cellDims[1] ) * boxes->cellDims[1];
        }
        else
        {
            if (ypos > boxes->maxPos[1])
            {
                /* put it in the end box */
                ypos = boxes->maxPos[1] - 0.5 * boxes->boxWidth[1];
            }
            else
            {
                /* put it in the end box */
                ypos = boxes->minPos[1] + 0.5 * boxes->boxWidth[1];
            }
        }
    }

    if (zpos > boxes->maxPos[2] || zpos < boxes->minPos[2])
    {
        if (boxes->PBC[2] == 1)
        {
            /* wrap position */
            zpos = zpos - floor( zpos / boxes->cellDims[2] ) * boxes->cellDims[2];
        }
        else
        {
            if (zpos > boxes->maxPos[2])
            {
                /* put it in the end box */
                zpos = boxes->maxPos[2] - 0.5 * boxes->boxWidth[2];
            }
            else
            {
                /* put it in the end box */
                zpos = boxes->minPos[2] + 0.5 * boxes->boxWidth[2];
            }
        }
    }

    /* find box for atom */
    posintx = (int) ((xpos - boxes->minPos[0]) / boxes->boxWidth[0]);
    if (posintx >= boxes->NBoxes[0])
        posintx = boxes->NBoxes[0] - 1;

    posinty = (int) ((ypos - boxes->minPos[1]) / boxes->boxWidth[1]);
    if (posinty >= boxes->NBoxes[1])
        posinty = boxes->NBoxes[1] - 1;

    posintz = (int) ((zpos - boxes->minPos[2]) / boxes->boxWidth[2]);
    if (posintz >= boxes->NBoxes[2])
        posintz = boxes->NBoxes[2] - 1;

    /* this numbers by x then z then y */
    boxIndex = (int) (posintx + posintz * boxes->NBoxes[0] + posinty * boxes->NBoxes[0] * boxes->NBoxes[2]);

    if (boxIndex < 0)
    {
        printf("WARNING: boxIndexOfAtom (CLIB): boxIndex < 0: %d, %d %d %d\n", boxIndex, posintx, posinty, posintz);
        printf("         pos = %f %f %f, box widths = %f %f %f, NBoxes = %d %d %d\n", xpos, ypos, zpos, boxes->boxWidth[0], boxes->boxWidth[1],
                boxes->boxWidth[2], boxes->NBoxes[0], boxes->NBoxes[1], boxes->NBoxes[2]);
        printf("         min box pos = %f %f %f\n", boxes->minPos[0], boxes->minPos[1], boxes->minPos[2]);
    }

    return boxIndex;
}


/*******************************************************************************
 ** return i, j, k indices of box
 *******************************************************************************/
void boxIJKIndices(int dim1, int* ijkIndices, int boxIndex, struct Boxes *boxes)
{
    int xint, yint, zint, tmp;


    /* compute the indices */
    yint = (int) ( boxIndex / (boxes->NBoxes[0] * boxes->NBoxes[2]) );

    tmp = boxIndex - yint * boxes->NBoxes[0] * boxes->NBoxes[2];
    zint = (int) ( tmp / boxes->NBoxes[0] );

    xint = (int) (tmp - zint * boxes->NBoxes[0]);

    ijkIndices[0] = xint;
    ijkIndices[1] = yint;
    ijkIndices[2] = zint;
}


/*******************************************************************************
 ** returns the box index of box with given i,j,k indices
 *******************************************************************************/
int boxIndexFromIJK(int xindex, int yindex, int zindex, struct Boxes *boxes)
{
    int xint, yint, zint, box;


    xint = xindex;
    zint = boxes->NBoxes[0] * zindex;
    yint = boxes->NBoxes[0] * boxes->NBoxes[2] * yindex;

    box = (int) (xint + yint + zint);

    return box;
}


/*******************************************************************************
 ** returns neighbourhood of given box
 *******************************************************************************/
int getBoxNeighbourhood(int mainBox, int* boxNeighbourList, struct Boxes *boxes)
{
    int mainBoxIJK[3];
    int i, count;
    int lstart[3], lfinish[3];


    /* first get i,j,k indices of the main box */
    boxIJKIndices( 3, mainBoxIJK, mainBox, boxes );

    /* handle cases where we have less than three boxes in a direction */
    for (i = 0; i < 3; i++)
    {
        int num = boxes->NBoxes[i];

        if (num == 1)
        {
            lstart[i] = 0;
            lfinish[i] = 1;
        }
        else if (num == 2)
        {
            lstart[i] = 0;
            lfinish[i] = 2;
        }
        else
        {
            lstart[i] = -1;
            lfinish[i] = 2;
        }
    }

    /* loop in x direction */
    count = 0;
    for (i = lstart[0]; i < lfinish[0]; i++)
    {
        int j;
        int posintx = mainBoxIJK[0] + i;

        /* wrap as required */
        if (posintx < 0)
            posintx += boxes->NBoxes[0];
        else if (posintx >= boxes->NBoxes[0])
            posintx -= boxes->NBoxes[0];

        /* loop in y direction */
        for (j = lstart[1]; j < lfinish[1]; j++)
        {
            int k;
            int posinty = mainBoxIJK[1] + j;

            /* wrap as required */
            if (posinty < 0)
                posinty += boxes->NBoxes[1];
            else if (posinty >= boxes->NBoxes[1])
                posinty -= boxes->NBoxes[1];

            /* loop in z direction */
            for (k = lstart[2]; k < lfinish[2]; k++)
            {
                int index;
                int posintz = mainBoxIJK[2] + k;

                /* wrap */
                if (posintz < 0)
                    posintz += boxes->NBoxes[2];
                else if (posintz >= boxes->NBoxes[2])
                    posintz -= boxes->NBoxes[2];

                /* get index of box from this i,j,k */
                index = boxIndexFromIJK(posintx, posinty, posintz, boxes);
                boxNeighbourList[count] = index;
                count++;
            }
        }
    }

    return count;
}


/*******************************************************************************
 ** free boxes memory
 *******************************************************************************/
void freeBoxes(struct Boxes *boxes)
{
    int i;


    for (i = 0; i < boxes->totNBoxes; i++)
    {
        if (boxes->boxNAtoms[i]) free(boxes->boxAtoms[i]);
    }
    free(boxes->boxAtoms);
    free(boxes->boxNAtoms);
    free(boxes);
}
