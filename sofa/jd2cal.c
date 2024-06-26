#include "sofa.h"
#include "sofam.h"
#include <float.h>

int iauJd2cal(double dj1, double dj2,
              int *iy, int *im, int *id, double *fd)
/*
**  - - - - - - - - - -
**   i a u J d 2 c a l
**  - - - - - - - - - -
**
**  Julian Date to Gregorian year, month, day, and fraction of a day.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     dj1,dj2   double   Julian Date (Notes 1, 2)
**
**  Returned (arguments):
**     iy        int      year
**     im        int      month
**     id        int      day
**     fd        double   fraction of day
**
**  Returned (function value):
**               int      status:
**                           0 = OK
**                          -1 = unacceptable date (Note 1)
**
**  Notes:
**
**  1) The earliest valid date is -68569.5 (-4900 March 1).  The
**     largest value accepted is 1e9.
**
**  2) The Julian Date is apportioned in any convenient way between
**     the arguments dj1 and dj2.  For example, JD=2450123.7 could
**     be expressed in any of these ways, among others:
**
**            dj1             dj2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     Separating integer and fraction uses the "compensated summation"
**     algorithm of Kahan-Neumaier to preserve as much precision as
**     possible irrespective of the jd1+jd2 apportionment.
**
**  3) In early eras the conversion is from the "proleptic Gregorian
**     calendar";  no account is taken of the date(s) of adoption of
**     the Gregorian calendar, nor is the AD/BC numbering convention
**     observed.
**
**  References:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 12.92 (p604).
**
**     Klein, A., A Generalized Kahan-Babuska-Summation-Algorithm.
**     Computing, 76, 279-293 (2006), Section 3.
**
**  This revision:  2021 May 11
**
**  SOFA release 2021-05-12
**
**  Copyright (C) 2021 IAU SOFA Board.  See notes at end.
*/
{
/* Minimum and maximum allowed JD */
    const double DJMIN = -68569.5;
    const double DJMAX = 1e9;

    long jd, i, l, n, k;
    double dj, f1, f2, d, s, cs, v[2], x, t, f;


/* Verify date is acceptable. */
    dj = dj1 + dj2;
    if (dj < DJMIN || dj > DJMAX) return -1;

/* Separate day and fraction (where -0.5 <= fraction < 0.5). */
    d = dnint(dj1);
    f1 = dj1 - d;
    jd = (long) d;
    d = dnint(dj2);
    f2 = dj2 - d;
    jd += (long) d;

/* Compute f1+f2+0.5 using compensated summation (Klein 2006). */
    s = 0.5;
    cs = 0.0;
    v[0] = f1;
    v[1] = f2;
    for (i = 0; i < 2; i++) {
        x = v[i];
        t = s + x;
        cs += fabs(s) >= fabs(x) ? (s - t) + x : (x - t) + s;
        s = t;
        if (s >= 1.0) {
            jd++;
            s -= 1.0;
        }
    }
    f = s + cs;
    cs = f - s;

/* Deal with negative f. */
    if (f < 0.0) {

        /* Compensated summation: assume that |s| <= 1.0. */
        f = s + 1.0;
        cs += (1.0 - f) + s;
        s = f;
        f = s + cs;
        cs = f - s;
        jd--;
    }

/* Deal with f that is 1.0 or more (when rounded to double). */
    if ((f - 1.0) >= -DBL_EPSILON / 4.0) {

        /* Compensated summation: assume that |s| <= 1.0. */
        t = s - 1.0;
        cs += (s - t) - 1.0;
        s = t;
        f = s + cs;
        if (-DBL_EPSILON / 2.0 < f) {
            jd++;
            f = gmax(f, 0.0);
        }
    }

/* Express day in Gregorian calendar. */
    l = jd + 68569L;
    n = (4L * l) / 146097L;
    l -= (146097L * n + 3L) / 4L;
    i = (4000L * (l + 1L)) / 1461001L;
    l -= (1461L * i) / 4L - 31L;
    k = (80L * l) / 2447L;
    *id = (int) (l - (2447L * k) / 80L);
    l = k / 11L;
    *im = (int) (k + 2L - 12L * l);
    *iy = (int) (100L * (n - 49L) + i + l);
    *fd = f;

/* Success. */
    return 0;

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
