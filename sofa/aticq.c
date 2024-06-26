#include "sofa.h"

void iauAticq(double ri, double di, iauASTROM *astrom,
              double *rc, double *dc)
/*
**  - - - - - - - - -
**   i a u A t i c q
**  - - - - - - - - -
**
**  Quick CIRS RA,Dec to ICRS astrometric place, given the star-
**  independent astrometry parameters.
**
**  Use of this function is appropriate when efficiency is important and
**  where many star positions are all to be transformed for one date.
**  The star-independent astrometry parameters can be obtained by
**  calling one of the functions iauApci[13], iauApcg[13], iauApco[13]
**  or iauApcs[13].
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     ri,di  double     CIRS RA,Dec (radians)
**     astrom iauASTROM* star-independent astrometry parameters:
**      pmt    double       PM time interval (SSB, Julian years)
**      eb     double[3]    SSB to observer (vector, au)
**      eh     double[3]    Sun to observer (unit vector)
**      em     double       distance from Sun to observer (au)
**      v      double[3]    barycentric observer velocity (vector, c)
**      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
**      bpn    double[3][3] bias-precession-nutation matrix
**      along  double       longitude + s' (radians)
**      xpl    double       polar motion xp wrt local meridian (radians)
**      ypl    double       polar motion yp wrt local meridian (radians)
**      sphi   double       sine of geodetic latitude
**      cphi   double       cosine of geodetic latitude
**      diurab double       magnitude of diurnal aberration vector
**      eral   double       "local" Earth rotation angle (radians)
**      refa   double       refraction constant A (radians)
**      refb   double       refraction constant B (radians)
**
**  Returned:
**     rc,dc  double     ICRS astrometric RA,Dec (radians)
**
**  Notes:
**
**  1) Only the Sun is taken into account in the light deflection
**     correction.
**
**  2) Iterative techniques are used for the aberration and light
**     deflection corrections so that the functions iauAtic13 (or
**     iauAticq) and iauAtci13 (or iauAtciq) are accurate inverses;
**     even at the edge of the Sun's disk the discrepancy is only about
**     1 nanoarcsecond.
**
**  Called:
**     iauS2c       spherical coordinates to unit vector
**     iauTrxp      product of transpose of r-matrix and p-vector
**     iauZp        zero p-vector
**     iauAb        stellar aberration
**     iauLdsun     light deflection by the Sun
**     iauC2s       p-vector to spherical
**     iauAnp       normalize angle into range +/- pi
**
**  This revision:   2013 October 9
**
**  SOFA release 2021-05-12
**
**  Copyright (C) 2021 IAU SOFA Board.  See notes at end.
*/
{
    int j, i;
    double pi[3], ppr[3], pnat[3], pco[3], w, d[3], before[3], r2, r,
            after[3];


/* CIRS RA,Dec to Cartesian. */
    iauS2c(ri, di, pi);

/* Bias-precession-nutation, giving GCRS proper direction. */
    iauTrxp(astrom->bpn, pi, ppr);

/* Aberration, giving GCRS natural direction. */
    iauZp(d);
    for (j = 0; j < 2; j++) {
        r2 = 0.0;
        for (i = 0; i < 3; i++) {
            w = ppr[i] - d[i];
            before[i] = w;
            r2 += w * w;
        }
        r = sqrt(r2);
        for (i = 0; i < 3; i++) {
            before[i] /= r;
        }
        iauAb(before, astrom->v, astrom->em, astrom->bm1, after);
        r2 = 0.0;
        for (i = 0; i < 3; i++) {
            d[i] = after[i] - before[i];
            w = ppr[i] - d[i];
            pnat[i] = w;
            r2 += w * w;
        }
        r = sqrt(r2);
        for (i = 0; i < 3; i++) {
            pnat[i] /= r;
        }
    }

/* Light deflection by the Sun, giving BCRS coordinate direction. */
    iauZp(d);
    for (j = 0; j < 5; j++) {
        r2 = 0.0;
        for (i = 0; i < 3; i++) {
            w = pnat[i] - d[i];
            before[i] = w;
            r2 += w * w;
        }
        r = sqrt(r2);
        for (i = 0; i < 3; i++) {
            before[i] /= r;
        }
        iauLdsun(before, astrom->eh, astrom->em, after);
        r2 = 0.0;
        for (i = 0; i < 3; i++) {
            d[i] = after[i] - before[i];
            w = pnat[i] - d[i];
            pco[i] = w;
            r2 += w * w;
        }
        r = sqrt(r2);
        for (i = 0; i < 3; i++) {
            pco[i] /= r;
        }
    }

/* ICRS astrometric RA,Dec. */
    iauC2s(pco, &w, dc);
    *rc = iauAnp(w);

/* Finished. */

/*----------------------------------------------------------------------
**
**  Copyright (C) 2021
**  Standards Of Fundamental Astronomy Board
**  of the International Astronomical Union.
**
**  =====================
**  SOFA Software License
**  =====================
**
**  NOTICE TO USER:
**
**  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING SIX TERMS AND
**  CONDITIONS WHICH APPLY TO ITS USE.
**
**  1. The Software is owned by the IAU SOFA Board ("SOFA").
**
**  2. Permission is granted to anyone to use the SOFA software for any
**     purpose, including commercial applications, free of charge and
**     without payment of royalties, subject to the conditions and
**     restrictions listed below.
**
**  3. You (the user) may copy and distribute SOFA source code to others,
**     and use and adapt its code and algorithms in your own software,
**     on a world-wide, royalty-free basis.  That portion of your
**     distribution that does not consist of intact and unchanged copies
**     of SOFA source code files is a "derived work" that must comply
**     with the following requirements:
**
**     a) Your work shall be marked or carry a statement that it
**        (i) uses routines and computations derived by you from
**        software provided by SOFA under license to you; and
**        (ii) does not itself constitute software provided by and/or
**        endorsed by SOFA.
**
**     b) The source code of your derived work must contain descriptions
**        of how the derived work is based upon, contains and/or differs
**        from the original SOFA software.
**
**     c) The names of all routines in your derived work shall not
**        include the prefix "iau" or "sofa" or trivial modifications
**        thereof such as changes of case.
**
**     d) The origin of the SOFA components of your derived work must
**        not be misrepresented;  you must not claim that you wrote the
**        original software, nor file a patent application for SOFA
**        software or algorithms embedded in the SOFA software.
**
**     e) These requirements must be reproduced intact in any source
**        distribution and shall apply to anyone to whom you have
**        granted a further right to modify the source code of your
**        derived work.
**
**     Note that, as originally distributed, the SOFA software is
**     intended to be a definitive implementation of the IAU standards,
**     and consequently third-party modifications are discouraged.  All
**     variations, no matter how minor, must be explicitly marked as
**     such, as explained above.
**
**  4. You shall not cause the SOFA software to be brought into
**     disrepute, either by misuse, or use for inappropriate tasks, or
**     by inappropriate modification.
**
**  5. The SOFA software is provided "as is" and SOFA makes no warranty
**     as to its use or performance.   SOFA does not and cannot warrant
**     the performance or results which the user may obtain by using the
**     SOFA software.  SOFA makes no warranties, express or implied, as
**     to non-infringement of third party rights, merchantability, or
**     fitness for any particular purpose.  In no event will SOFA be
**     liable to the user for any consequential, incidental, or special
**     damages, including any lost profits or lost savings, even if a
**     SOFA representative has been advised of such damages, or for any
**     claim by any third party.
**
**  6. The provision of any version of the SOFA software under the terms
**     and conditions specified herein does not imply that future
**     versions will also be made available under the same terms and
**     conditions.
*
**  In any published work or commercial product which uses the SOFA
**  software directly, acknowledgement (see www.iausofa.org) is
**  appreciated.
**
**  Correspondence concerning SOFA software should be addressed as
**  follows:
**
**      By email:  sofa@ukho.gov.uk
**      By post:   IAU SOFA Center
**                 HM Nautical Almanac Office
**                 UK Hydrographic Office
**                 Admiralty Way, Taunton
**                 Somerset, TA1 2DN
**                 United Kingdom
**
**--------------------------------------------------------------------*/
}
