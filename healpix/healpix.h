//#include "healpix_base.h"
#pragma once

#ifndef HEALPIX_H
#define HEALPIX_H

#include "vec3.h"
#include "arr.h"

#ifndef _MYDEFS_H
#define _MYDEFS_H

#include <math.h>
#include <stdint.h>

// these are no longer in math.h for the c99 standard, so
// define these if not already defined

#define wlog(...) fprintf(stderr, __VA_ARGS__)

#ifndef M_PI
# define M_E		2.7182818284590452354	/* e */
# define M_LOG2E	1.4426950408889634074	/* log_2 e */
# define M_LOG10E	0.43429448190325182765	/* log_10 e */
# define M_LN2		0.69314718055994530942	/* log_e 2 */
# define M_LN10		2.30258509299404568402	/* log_e 10 */
# define M_PI		3.14159265358979323846	/* pi */
# define M_PI_2		1.57079632679489661923	/* pi/2 */
# define M_PI_4		0.78539816339744830962	/* pi/4 */
# define M_1_PI		0.31830988618379067154	/* 1/pi */
# define M_2_PI		0.63661977236758134308	/* 2/pi */
# define M_2_SQRTPI	1.12837916709551257390	/* 2/sqrt(pi) */
# define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
# define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */

#endif
#define I int


# define M_TWO_PI   6.28318530717958647693 /* 2*pi */
# define M_TWOTHIRD 0.66666666666666666666
static const double twothird = 2.0 / 3.0;
static const double pi = 3.141592653589793238462643383279502884197;
static const double twopi = 6.283185307179586476925286766559005768394;
static const double halfpi = 1.570796326794896619231321691639751442099;
static const double inv_halfpi = 0.6366197723675813430755350534900574;
const double inv_twopi = 1.0 / twopi;
#define D2R  0.017453292519943295
#define R2D  57.295779513082323
static const short utab[] = {
#define Z(a) 0x##a##0, 0x##a##1, 0x##a##4, 0x##a##5
#define Y(a) Z(a##0), Z(a##1), Z(a##4), Z(a##5)
#define X(a) Y(a##0), Y(a##1), Y(a##4), Y(a##5)
        X(0), X(1), X(4), X(5)
#undef X
#undef Y
#undef Z
};
typedef int64_t int64;

#define SYSTEM_EQ 0

#endif
using namespace std;
I nside_, npface_, ncap_, npix_;
I order_;
double fact1_, fact2_;

int hpix_radec_degrees_to_thetaphi_radians(double ra, double dec, double *theta, double *phi) {

    int status = 0;

    if (ra < 0.0 || ra > 360.) {
        goto _hpix_conv_bail;
    }
    if (dec < -90. || dec > 90.) {
        goto _hpix_conv_bail;
    }

    *phi = ra * D2R;
    *theta = -dec * D2R + M_PI_2;

    status = 1;

    _hpix_conv_bail:

    return status;
}

/*
int64 hpix_eq2pix(double ra, double dec,int *status) {
    int64 nside=1;
    int64 ipix=0;
    double theta=0, phi=0;
    *status = hpix_radec_degrees_to_thetaphi_radians(ra, dec, &theta, &phi);
    if (! (*status) ) {
        return ipix;
    }
    double z = cos(theta);
    double za = fabs(z);
    // in [0,4)
    double tt = fmod(phi, M_TWO_PI)/M_PI_2;
    if (za <= M_TWOTHIRD) {
        double temp1 = nside*(.5 + tt);
        double temp2 = nside*.75*z;
        int64 jp = (int64)(temp1-temp2); // index of  ascending edge content
        int64 jm = (int64)(temp1+temp2); // index of descending edge content
        int64 ir = nside + 1 + jp - jm;  // in {1,2n+1} (ring number counted from z=2/3)
        int64 kshift = 1 - (ir % 2);      // kshift=1 if ir even, 0 otherwise
        int64 nl4 = 4*nside;
        int64 ip = (int64)( ( jp+jm - nside + kshift + 1 ) / 2); // in {0,4n-1}
        ip = ip % nl4;
        ipix = 2*nside*(nside-1) + nl4*(ir-1) + ip;
    } else {
        // North & South polar caps
        double tp = tt - (int64)(tt);   // MODULO(tt,1.0_dp)
        double tmp = nside * sqrt( 3.0*(1.0 - za) );
        int64 jp = (int64)(tp*tmp);              // increasing edge content index
        int64 jm = (int64)((1.0 - tp) * tmp); // decreasing edge content index
        int64 ir = jp + jm + 1;        // ring number counted from the closest pole
        int64 ip = (int64)( tt * ir);     // in {0,4*ir-1}
        if (ip >= 4*ir) {
            ip = ip - 4*ir;
        }
        if (z>0.) {
            ipix = 2*ir*(ir-1) + ip;
        } else {
            ipix = 12*nside*nside; - 2*ir*(ir+1) + ip;
        }
    }
    return ipix;
}
 */
static int xyf2nest(int nside, int ix, int iy, int face_num) {
    return (face_num * nside * nside) +
           (utab[ix & 0xff] | (utab[ix >> 8] << 16)
            | (utab[iy & 0xff] << 1) | (utab[iy >> 8] << 17));
}

//static double fmodulo (double v1, double v2)
//{
//    if (v1>=0)
//        return (v1<v2) ? v1 : fmod(v1,v2);
//    double tmp=fmod(v1,v2)+v2;
//    return (tmp==v2) ? 0. : tmp;
//}
static int ang2pix_nest_z_phi(long nside_, double z, double phi) {
    double za = fabs(z);
    double tt = fmodulo(phi, twopi) * inv_halfpi; /* in [0,4) */
    int face_num, ix, iy;

    if (za <= twothird) /* Equatorial region */
    {
        double temp1 = nside_ * (0.5 + tt);
        double temp2 = nside_ * (z * 0.75);
        int jp = (int) (temp1 - temp2); /* index of  ascending edge content */
        int jm = (int) (temp1 + temp2); /* index of descending edge content */
        int ifp = jp / nside_;  /* in {0,4} */
        int ifm = jm / nside_;
        face_num = (ifp == ifm) ? (ifp | 4) : ((ifp < ifm) ? ifp : (ifm + 8));

        ix = jm & (nside_ - 1);
        iy = nside_ - (jp & (nside_ - 1)) - 1;
    } else /* polar region, za > 2/3 */
    {
        int ntt = (int) tt, jp, jm;
        double tp, tmp;
        if (ntt >= 4) ntt = 3;
        tp = tt - ntt;
        tmp = nside_ * sqrt(3 * (1 - za));

        jp = (int) (tp * tmp); /* increasing edge content index */
        jm = (int) ((1.0 - tp) * tmp); /* decreasing edge content index */
        if (jp >= nside_) jp = nside_ - 1; /* for points too close to the boundary */
        if (jm >= nside_) jm = nside_ - 1;
        if (z >= 0) {
            face_num = ntt;  /* in {0,3} */
            ix = nside_ - jm - 1;
            iy = nside_ - jp - 1;
        } else {
            face_num = ntt + 8; /* in {8,11} */
            ix = jp;
            iy = jm;
        }
    }

    return xyf2nest(nside_, ix, iy, face_num);
}

void ang2pix_nest(long nside, double theta, double phi, long *ipix) {
    //   UTIL_ASSERT((theta>=0)&&(theta<=pi),"theta completetree of range");
    *ipix = ang2pix_nest_z_phi(nside, cos(theta), phi);
}

long eq2pix_nest(long nside, double ra, double dec) {
    long t;
    double theta = 0, phi = 0;

    hpix_radec_degrees_to_thetaphi_radians(ra, dec, &theta, &phi);
    ang2pix_nest(nside, theta, phi, &t);
    return t;
}

static int isqrt(int v) { return (int) (sqrt(v + 0.5)); }

static void pix2ang_ring_z_phi(int nside_, int pix, double *z, double *phi) {
    long ncap_ = nside_ * (nside_ - 1) * 2;
    long npix_ = 12 * nside_ * nside_;
    double fact2_ = 4. / npix_;
    if (pix < ncap_) /* North Polar cap */
    {
        int iring = (1 + isqrt(1 + 2 * pix)) >> 1; /* counted from North pole */
        int iphi = (pix + 1) - 2 * iring * (iring - 1);

        *z = 1.0 - (iring * iring) * fact2_;
        *phi = (iphi - 0.5) * halfpi / iring;
    } else if (pix < (npix_ - ncap_)) /* Equatorial region */
    {
        double fact1_ = (nside_ << 1) * fact2_;
        int ip = pix - ncap_;
        int iring = ip / (4 * nside_) + nside_; /* counted from North pole */
        int iphi = ip % (4 * nside_) + 1;
        /* 1 if iring+nside is odd, 1/2 otherwise */
        double fodd = ((iring + nside_) & 1) ? 1 : 0.5;

        int nl2 = 2 * nside_;
        *z = (nl2 - iring) * fact1_;
        *phi = (iphi - fodd) * pi / nl2;
    } else /* South Polar cap */
    {
        int ip = npix_ - pix;
        int iring = (1 + isqrt(2 * ip - 1)) >> 1; /* counted from South pole */
        int iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));

        *z = -1.0 + (iring * iring) * fact2_;
        *phi = (iphi - 0.5) * halfpi / iring;
    }
}

static double haversine(double ra1, double dec1,
                        double ra2, double dec2) {
    // cout<<"dec"<<lat1<<"ra"<<lon1<<"dec"<<lat2<<"ra"<<lon2<<endl;
    // distance between latitudes
    // and longitudes
    double lat1=dec1;
    double lon1=ra1;
    double lat2=dec2;
    double lon2=ra2;
    double dLat = (lat2 - lat1) *
                  M_PI / 180.0;
    double dLon = (lon2 - lon1) *
                  M_PI / 180.0;

    // convert to radians
    lat1 = (lat1) * M_PI / 180.0;
    lat2 = (lat2) * M_PI / 180.0;

    // apply formulae
    double a = pow(sin(dLat / 2), 2) +
               pow(sin(dLon / 2), 2) *
               cos(lat1) * cos(lat2);

    double r = 2 * asin(sqrt(a));//between 0 and 2*pi
    //convert to degree
    double rad = std::fmod(r, 2.0 * pi);
    if (rad < 0)
    {
        rad *= -1;
    }
    double deg = rad * (180. / pi);

    return deg;
}

static const short ctab[] = {
#define Z(a) a,a+1,a+256,a+257
#define Y(a) Z(a),Z(a+2),Z(a+512),Z(a+514)
#define X(a) Y(a),Y(a+4),Y(a+1024),Y(a+1028)
        X(0), X(8), X(2048), X(2056)
#undef X
#undef Y
#undef Z
};

static void nest2xyf(int nside, int pix, int *ix, int *iy, int *face_num) {
    int npface_ = nside * nside, raw;
    *face_num = pix / npface_;
    pix &= (npface_ - 1);
    raw = (pix & 0x5555) | ((pix & 0x55550000) >> 15);
    *ix = ctab[raw & 0xff] | (ctab[raw >> 8] << 4);
    pix >>= 1;
    raw = (pix & 0x5555) | ((pix & 0x55550000) >> 15);
    *iy = ctab[raw & 0xff] | (ctab[raw >> 8] << 4);
}

static const int jrll[] = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};
static const int jpll[] = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7};

static void pix2ang_nest_z_phi(int nside_, int pix, double *z, double *phi) {
    int nl4 = nside_ * 4;
    int npix_ = 12 * nside_ * nside_;
    double fact2_ = 4. / npix_;
    int face_num, ix, iy, nr, kshift;

    nest2xyf(nside_, pix, &ix, &iy, &face_num);
    int jr = (jrll[face_num] * nside_) - ix - iy - 1;

    if (jr < nside_) {
        nr = jr;
        *z = 1 - nr * nr * fact2_;
        kshift = 0;
    } else if (jr > 3 * nside_) {
        nr = nl4 - jr;
        *z = nr * nr * fact2_ - 1;
        kshift = 0;
    } else {
        double fact1_ = (nside_ << 1) * fact2_;
        nr = nside_;
        *z = (2 * nside_ - jr) * fact1_;
        kshift = (jr - nside_) & 1;
    }

    int jp = (jpll[face_num] * nr + ix - iy + 1 + kshift) / 2;
    if (jp > nl4) jp -= nl4;
    if (jp < 1) jp += nl4;

    *phi = (jp - (kshift + 1) * 0.5) * (halfpi / nr);
}

long nside2npix(const long nside) { return 12 * nside * nside; }

void pix2ang_nest(long nside, long ipix, double *theta, double *phi) {
    double z;
    pix2ang_nest_z_phi(nside, ipix, &z, phi);
    *theta = acos(z);
}

void find(int nside, double ra, double dec, double r) {
    int tot = nside2npix(nside);
    int tot2 = 0;
    for (int i = 0; i < tot; i++) {
        double theta, phi;
        pix2ang_nest(nside, i, &theta, &phi);
        double dec2 = (theta - M_PI_2) / (-D2R);
        double ra2 = phi / D2R;
        long t;
        ang2pix_nest(nside, theta, phi, &t);
        if (haversine(dec, ra, dec2, ra2) < (r * 2 * pi) / 360) {
            tot2++;
        }
    }
}

double calcu_dis(int nside, int pix, double ra, double dec) {
    //  cout << nside << "  " << pix << "  " << ra << "  " << dec << "  " << endl;

    double theta, phi;
    pix2ang_nest(nside, pix, &theta, &phi);
    double dec2 = (theta - M_PI_2) / (-D2R);
    double ra2 = phi / D2R;

    double tot = haversine(dec, ra, dec2, ra2) * 360 / (2 * pi);
    //   cout << totInMemory << endl;
    return tot;
}

void find2(int nside, double ra, double dec, double r) {
    int tot = nside2npix(nside);
    int tot2 = 0;
    for (int i = 0; i < tot; i++) {
        double theta, phi;
        pix2ang_nest(nside, i, &theta, &phi);
        double dec2 = (theta - M_PI_2) / (-D2R);
        double ra2 = phi / D2R;
        long t;
        ang2pix_nest(nside, theta, phi, &t);
        // cout<<i<<" "<<ra2<<" "<<dec2<<" "<<haversine(dec,ra,dec2,ra2)<<" "<<eq2pix(2,ra2,dec2)<<" "<<t<<endl;
        if (calcu_dis(nside, i, ra, dec) < r) {
            tot2++;
        }
    }
}

inline double v_angle(const vec3 &v1, const vec3 &v2) {
    using namespace std;
    return atan2(crossprod(v1, v2).Length(), dotprod(v1, v2));
}

double max_pixrad(int nside_) {
    vec3 va, vb;
    va.set_z_phi(2. / 3., pi / (4 * nside_));
    double t1 = 1. - 1. / nside_;
    t1 *= t1;
    vb.set_z_phi(1 - t1 / 3, 0);
    return v_angle(va, vb);
}

vector<int> healpix_vector;

void append(int start, int end, int id = 0) {
//    cout<<start<<" "<<end<<" "<<id<<endl;
    for (int i = start; i < end; i++) {
        healpix_vector.push_back(i);
    }
}

I loc2pix(double z, double phi, double sth, bool have_sth) {
    double za = abs(z);
    double tt = fmodulo(phi * inv_halfpi, 4.0); // in [0,4)

    if (1)//scheme_==RING
    {
        if (za <= twothird) // Equatorial region
        {
            I nl4 = 4 * nside_;
            double temp1 = nside_ * (0.5 + tt);
            double temp2 = nside_ * z * 0.75;
            I jp = I(temp1 - temp2); // index of  ascending edge content
            I jm = I(temp1 + temp2); // index of descending edge content

            // ring number counted from z=2/3
            I ir = nside_ + 1 + jp - jm; // in {1,2n+1}
            I kshift = 1 - (ir & 1); // kshift=1 if ir even, 0 otherwise

            I t1 = jp + jm - nside_ + kshift + 1 + nl4 + nl4;
            I ip = (order_ > 0) ?
                   (t1 >> 1) & (nl4 - 1) : ((t1 >> 1) % nl4); // in {0,4n-1}

            return ncap_ + (ir - 1) * nl4 + ip;
        } else  // North & South polar caps
        {
            double tp = tt - I(tt);
            double tmp = ((za < 0.99) || (!have_sth)) ?
                         nside_ * sqrt(3 * (1 - za)) :
                         nside_ * sth / sqrt((1. + za) / 3.);

            I jp = I(tp * tmp); // increasing edge content index
            I jm = I((1.0 - tp) * tmp); // decreasing edge content index

            I ir = jp + jm + 1; // ring number counted from the closest pole
            I ip = I(tt * ir); // in {0,4*ir-1}
//            planck_assert((ip >= 0) && (ip < 4 * ir), "must not happen");
            //ip %= 4*ir;

            return (z > 0) ? 2 * ir * (ir - 1) + ip : npix_ - 2 * ir * (ir + 1) + ip;
        }
    } else // scheme_ == NEST
    {
        if (za <= twothird) // Equatorial region
        {
            double temp1 = nside_ * (0.5 + tt);
            double temp2 = nside_ * (z * 0.75);
            I jp = I(temp1 - temp2); // index of  ascending edge content
            I jm = I(temp1 + temp2); // index of descending edge content
            I ifp = jp >> order_;  // in {0,4}
            I ifm = jm >> order_;
            int face_num = (ifp == ifm) ? (ifp | 4) : ((ifp < ifm) ? ifp : (ifm + 8));

            int ix = jm & (nside_ - 1),
                    iy = nside_ - (jp & (nside_ - 1)) - 1;
            return xyf2nest(nside_, ix, iy, face_num);
        } else // polar region, za > 2/3
        {
            int ntt = min(3, int(tt));
            double tp = tt - ntt;
            double tmp = ((za < 0.99) || (!have_sth)) ?
                         nside_ * sqrt(3 * (1 - za)) :
                         nside_ * sth / sqrt((1. + za) / 3.);

            I jp = I(tp * tmp); // increasing edge content index
            I jm = I((1.0 - tp) * tmp); // decreasing edge content index
            jp = min(jp, nside_ - 1); // for points too close to the boundary
            jm = min(jm, nside_ - 1);
            return (z >= 0) ?
                   xyf2nest(nside_, nside_ - jm - 1, nside_ - jp - 1, ntt) : xyf2nest(nside_, jp, jm, ntt + 8);
        }
    }
}

I zphi2pix(double z, double phi) { return loc2pix(z, phi, 0., false); }

I ring_above(double z) {
    double az = abs(z);
    if (az <= twothird) // equatorial region
        return I(nside_ * (2 - 1.5 * z));
    I iring = I(nside_ * sqrt(3 * (1 - az)));
    return (z > 0) ? iring : 4 * nside_ - iring - 1;
}

void get_ring_info_small(I ring, I &startpix, I &ringpix, bool &shifted) {
    if (ring < nside_) {
        shifted = true;
        ringpix = 4 * ring;
        startpix = 2 * ring * (ring - 1);
    } else if (ring < 3 * nside_) {
        shifted = ((ring - nside_) & 1) == 0;
        ringpix = 4 * nside_;
        startpix = ncap_ + (ring - nside_) * ringpix;
    } else {
        shifted = true;
        I nr = 4 * nside_ - ring;
        ringpix = 4 * nr;
        startpix = npix_ - 2 * nr * (nr + 1);
    }
}

double ring2z(I ring) {
    if (ring < nside_)
        return 1 - ring * ring * fact2_;
    if (ring <= 3 * nside_)
        return (2 * nside_ - ring) * fact1_;
    ring = 4 * nside_ - ring;
    return ring * ring * fact2_ - 1;
}

static inline int special_div(int a, int b) {
    int t = (a >= (b << 1));
    a -= t * (b << 1);
    return (t << 1) + (a >= b);
}

static void ring2xyf(int nside_, int pix, int *ix, int *iy, int *face_num) {
    int iring, iphi, kshift, nr;
    int ncap_ = 2 * nside_ * (nside_ - 1);
    int npix_ = 12 * nside_ * nside_;
    int nl2 = 2 * nside_;

    if (pix < ncap_) /* North Polar cap */
    {
        iring = (1 + isqrt(1 + 2 * pix)) >> 1; /* counted from North pole */
        iphi = (pix + 1) - 2 * iring * (iring - 1);
        kshift = 0;
        nr = iring;
        *face_num = special_div(iphi - 1, nr);
    } else if (pix < (npix_ - ncap_)) /* Equatorial region */
    {
        int ip = pix - ncap_;
        iring = (ip / (4 * nside_)) + nside_; /* counted from North pole */
        iphi = (ip % (4 * nside_)) + 1;
        kshift = (iring + nside_) & 1;
        nr = nside_;
        unsigned int ire = iring - nside_ + 1;
        unsigned int irm = nl2 + 2 - ire;
        int ifm = (iphi - ire / 2 + nside_ - 1) / nside_;
        int ifp = (iphi - irm / 2 + nside_ - 1) / nside_;
        *face_num = (ifp == ifm) ? (ifp | 4) : ((ifp < ifm) ? ifp : (ifm + 8));
    } else /* South Polar cap */
    {
        int ip = npix_ - pix;
        iring = (1 + isqrt(2 * ip - 1)) >> 1; /* counted from South pole */
        iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
        kshift = 0;
        nr = iring;
        iring = 2 * nl2 - iring;
        *face_num = 8 + special_div(iphi - 1, nr);
    }

    int irt = iring - (jrll[*face_num] * nside_) + 1;
    int ipt = 2 * iphi - jpll[*face_num] * nr - kshift - 1;
    if (ipt >= nl2) ipt -= 8 * nside_;

    *ix = (ipt - irt) >> 1;
    *iy = (-(ipt + irt)) >> 1;
}

static int xyf2ring(int nside_, int ix, int iy, int face_num) {
    int nl4 = 4 * nside_;
    int jr = (jrll[face_num] * nside_) - ix - iy - 1, jp;

    int nr, kshift, n_before;
    if (jr < nside_) {
        nr = jr;
        n_before = 2 * nr * (nr - 1);
        kshift = 0;
    } else if (jr > 3 * nside_) {
        nr = nl4 - jr;
        n_before = 12 * nside_ * nside_ - 2 * (nr + 1) * nr;
        kshift = 0;
    } else {
        int ncap_ = 2 * nside_ * (nside_ - 1);
        nr = nside_;
        n_before = ncap_ + (jr - nside_) * nl4;
        kshift = (jr - nside_) & 1;
    }

    jp = (jpll[face_num] * nr + ix - iy + 1 + kshift) / 2;
    if (jp > nl4)
        jp -= nl4;
    else if (jp < 1) jp += nl4;

    return n_before + jp - 1;
}

//static void pix2ang_ring_z_phi(int nside_, int pix, double *z, double *phi) {
//    long ncap_ = nside_ * (nside_ - 1) * 2;
//    long npix_ = 12 * nside_ * nside_;
//    double fact2_ = 4. / npix_;
//    if (pix < ncap_) /* North Polar cap */
//    {
//        int iring = (1 + isqrt(1 + 2 * pix)) >> 1; /* counted from North pole */
//        int iphi = (pix + 1) - 2 * iring * (iring - 1);
//
//        *z = 1.0 - (iring * iring) * fact2_;
//        *phi = (iphi - 0.5) * halfpi / iring;
//    } else if (pix < (npix_ - ncap_)) /* Equatorial region */
//    {
//        double fact1_ = (nside_ << 1) * fact2_;
//        int ip = pix - ncap_;
//        int iring = ip / (4 * nside_) + nside_; /* counted from North pole */
//        int iphi = ip % (4 * nside_) + 1;
//        /* 1 if iring+nside is odd, 1/2 otherwise */
//        double fodd = ((iring + nside_) & 1) ? 1 : 0.5;
//
//        int nl2 = 2 * nside_;
//        *z = (nl2 - iring) * fact1_;
//        *phi = (iphi - fodd) * pi / nl2;
//    } else /* South Polar cap */
//    {
//        int ip = npix_ - pix;
//        int iring = (1 + isqrt(2 * ip - 1)) >> 1; /* counted from South pole */
//        int iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
//
//        *z = -1.0 + (iring * iring) * fact2_;
//        *phi = (iphi - 0.5) * halfpi / iring;
//    }
//}

inline double cosdist_zphi(double z1, double phi1, double z2, double phi2) {
    using namespace std;
    return z1 * z2 + cos(phi1 - phi2) * sqrt((1. - z1 * z1) * (1. - z2 * z2));
}

bool check_pixel_ring(I pix, I nr, I ipix1, int fct, double cz, double cphi, double cosrp2, I cpix) {
    if (pix >= nr) pix -= nr;
    if (pix < 0) pix += nr;
    pix += ipix1;
    if (pix == cpix) return false; // disk center in pixel => overlap
    int px, py, pf;

    ring2xyf(nside_, pix, &px, &py, &pf);
    for (int i = 0; i < fct - 1; ++i) // go along the 4 edges
    {
        I ox = fct * px, oy = fct * py;
        double pz, pphi;
        pix2ang_ring_z_phi(fct * nside_, xyf2ring(fct * nside_, ox + i, oy, pf), &pz, &pphi);
        if (cosdist_zphi(pz, pphi, cz, cphi) > cosrp2) // overlap
            return false;
        pix2ang_ring_z_phi(fct * nside_, xyf2ring(fct * nside_, ox + fct - 1, oy + i, pf), &pz, &pphi);
        if (cosdist_zphi(pz, pphi, cz, cphi) > cosrp2) // overlap
            return false;
        pix2ang_ring_z_phi(fct * nside_, xyf2ring(fct * nside_, ox + fct - 1 - i, oy + fct - 1, pf), &pz, &pphi);
        if (cosdist_zphi(pz, pphi, cz, cphi) > cosrp2) // overlap
            return false;
        pix2ang_ring_z_phi(fct * nside_, xyf2ring(fct * nside_, ox, oy + fct - 1 - i, pf), &pz, &pphi);
        if (cosdist_zphi(pz, pphi, cz, cphi) > cosrp2) // overlap
            return false;
    }
    return true;
}

inline int ilog2(I arg) {
#ifdef __GNUC__
    if (arg == 0) return 0;
    if (sizeof(I) == sizeof(int))
        return 8 * sizeof(int) - 1 - __builtin_clz(arg);
    if (sizeof(I) == sizeof(long))
        return 8 * sizeof(long) - 1 - __builtin_clzl(arg);
    if (sizeof(I) == sizeof(long long))
        return 8 * sizeof(long long) - 1 - __builtin_clzll(arg);
#endif
    int res = 0;
    while (arg > 0xFFFF) {
        res += 16;
        arg >>= 16;
    }
    if (arg > 0x00FF) {
        res |= 8;
        arg >>= 8;
    }
    if (arg > 0x000F) {
        res |= 4;
        arg >>= 4;
    }
    if (arg > 0x0003) {
        res |= 2;
        arg >>= 2;
    }
    if (arg > 0x0001) { res |= 1; }
    return res;
}

inline I ifloor(double arg) {
    using namespace std;
    return I(floor(arg));
}

int compress_bits(int v) {
    int raw = (v & 0x5555) | ((v & 0x55550000) >> 15);
    return ctab[raw & 0xff] | (ctab[raw >> 8] << 4);
}

int compress_bits(int64 v) {
    int64 raw = v & 0x5555555555555555ull;
    raw |= raw >> 15;
    return ctab[raw & 0xff] | (ctab[(raw >> 8) & 0xff] << 4)
           | (ctab[(raw >> 32) & 0xff] << 16) | (ctab[(raw >> 40) & 0xff] << 20);
}

void nest2xyf(I pix, int &ix, int &iy, int &face_num, int order_, int npface_) {
    face_num = pix >> (2 * order_);
    pix &= (npface_ - 1);
    ix = compress_bits(pix);
    iy = compress_bits(pix >> 1);
}

void pix2loc(I pix, double &z, double &phi, double &sth, bool &have_sth, I nside_) {
    I npface_, ncap_, npix_;
    I order_;
    double fact1_, fact2_;
    order_ = ((nside_) & (nside_ - 1)) ? -1 : ilog2(nside_);
    npface_ = nside_ * nside_;
    ncap_ = (npface_ - nside_) << 1;
    npix_ = 12 * npface_;
    fact2_ = 4. / npix_;
    fact1_ = (nside_ << 1) * fact2_;

    have_sth = false;
    int face_num, ix, iy;
    nest2xyf(pix, ix, iy, face_num, order_, npface_);

    I jr = (I(jrll[face_num]) << order_) - ix - iy - 1;

    I nr;
    if (jr < nside_) {
        nr = jr;
        double tmp = (nr * nr) * fact2_;
        z = 1 - tmp;
        if (z > 0.99) {
            sth = sqrt(tmp * (2. - tmp));
            have_sth = true;
        }
    } else if (jr > 3 * nside_) {
        nr = nside_ * 4 - jr;
        double tmp = (nr * nr) * fact2_;
        z = tmp - 1;
        if (z < -0.99) {
            sth = sqrt(tmp * (2. - tmp));
            have_sth = true;
        }
    } else {
        nr = nside_;
        z = (2 * nside_ - jr) * fact1_;
    }

    I tmp = I(jpll[face_num]) * nr + ix - iy;
//        planck_assert(tmp<8*nr,"must not happen");
    if (tmp < 0) tmp += 8 * nr;
    phi = (nr == nside_) ? 0.75 * halfpi * tmp * fact1_ :
          (0.5 * halfpi * tmp) / nr;

}

void pix2zphi_NESTED(I pix, double &z, double &phi, I nside_) {
    bool dum_b;
    double dum_d;
    pix2loc(pix, z, phi, dum_d, dum_b, nside_);
}

void thetaphi2xyz(double theta, double phi, double &x, double &y, double &z) {
    double st = sin(theta);
    x = st * cos(phi), y = st * sin(phi), z = cos(theta);
}

void append(I x) {
    append(x, x + 1);
}

/* Short note on the "zone":
   zone = 0: pixel lies completely outside the queried shape
          1: pixel may overlap with the shape, pixel center is outside
          2: pixel center is inside the shape, but maybe not the complete pixel
          3: pixel lies completely inside the shape */
void
check_pixel(int o, int order_, int omax, int zone, I pix, vector<pair<I, int> > &stk, bool inclusive, int &stacktop) {
    // cout<<"o "<<o<<"order_ "<<order_<<"omax "<<omax<<"zone "<<zone<<"pix "<<pix<<endl;
    if (zone == 0) return;

    if (o < order_) {
        if (zone >= 3) {
            int sdist = 2 * (order_ - o); // the "bit-shift distance" between map orders
            append(pix << sdist, (pix + 1) << sdist); // output all subpixels
        } else // (1<=zone<=2)
            for (int i = 0; i < 4; ++i)
                stk.push_back(make_pair(4 * pix + 3 - i, o + 1)); // add children
    } else if (o > order_) // this implies that inclusive==true
    {
        if (zone >= 2) // pixel center in shape
        {
            append(pix >> (2 * (o - order_))); // output the parent pixel at order_
            stk.resize(stacktop); // unwind the stack
        } else // (zone==1): pixel center in safety range
        {
            if (o < omax) // check sublevels
                for (int i = 0; i < 4; ++i) // add children in reverse order
                    stk.push_back(make_pair(4 * pix + 3 - i, o + 1));
            else // at resolution limit
            {
                append(pix >> (2 * (o - order_))); // output the parent pixel at order_
                stk.resize(stacktop); // unwind the stack
            }
        }
    } else // o==order_
    {
        if (zone >= 2)
            append(pix);
        else if (inclusive) // and (zone>=1)
        {
            if (order_ < omax) // check sublevels
            {
                stacktop = stk.size(); // remember current stack position
                for (int i = 0; i < 4; ++i) // add children in reverse order
                    stk.push_back(make_pair(4 * pix + 3 - i, o + 1));
            } else // at resolution limit
                append(pix); // output the pixel
        }
    }
}

void find3(double ra, double dec, double radius, int nside, int fact, bool ring) {//RING
    healpix_vector.clear();


    nside_ = nside;
    npface_ = nside_ * nside_;
    ncap_ = (npface_ - nside_) << 1;
    npix_ = 12 * npface_;
    fact2_ = 4. / npix_;
    fact1_ = (nside_ << 1) * fact2_;
    int order_max = 29;
    order_ = ((nside) & (nside - 1)) ? -1 : ilog2(nside);
    bool inclusive = (fact != 0);
    if (ring) {
        int fct = 1;
        int order_max = 29;
        order_ = ((nside) & (nside - 1)) ? -1 : ilog2(nside);
        if (inclusive) {
            // planck_assert (((I(1) << order_max) / nside_) >= fact,
            //             "invalid oversampling factor");
            fct = fact;
        }
        double rsmall, rbig;
        if (fct > 1) {

            rsmall = radius + max_pixrad(fct * nside_);
            rbig = radius + max_pixrad(nside_);
        } else
            rsmall = rbig = inclusive ? radius + max_pixrad(nside_) : radius;
        if (rsmall >= pi) {
            append(0, npix_, 1);
            return;
        }

        rbig = min(pi, rbig);

        double cosrsmall = cos(rsmall);
        double cosrbig = cos(rbig);
        double theta = 0, phi = 0;

        hpix_radec_degrees_to_thetaphi_radians(ra, dec, &theta, &phi);

        double z0 = cos(theta);
        double xa = 1. / sqrt((1 - z0) * (1 + z0));

        I cpix = zphi2pix(z0, phi);
        double rlat1 = theta - rsmall;
        double zmax = cos(rlat1);

        I irmin = ring_above(zmax) + 1;

        if ((rlat1 <= 0) && (irmin > 1)) // north pole in the disk
        {
            I sp, rp;
            bool dummy;
            get_ring_info_small(irmin - 1, sp, rp, dummy);
            append(0, sp + rp, 2);
        }

        if ((fct > 1) && (rlat1 > 0)) irmin = max(I(1), irmin - 1);

        double rlat2 = theta + rsmall;
        double zmin = cos(rlat2);
        I irmax = ring_above(zmin);

        if ((fct > 1) && (rlat2 < pi)) irmax = min(4 * nside_ - 1, irmax + 1);

        for (I iz = irmin; iz <= irmax; ++iz) {
            double z = ring2z(iz);
            double x = (cosrbig - z * z0) * xa;
            double ysq = 1 - z * z - x * x;
            double dphi = -1;
            if (ysq <= 0) // no intersection, ring completely inside or outside
                dphi = (fct == 1) ? 0 : pi - 1e-15;
            else
                dphi = atan2(sqrt(ysq), x);
            if (dphi > 0) {
                I nr, ipix1;
                bool shifted = false;
                get_ring_info_small(iz, ipix1, nr, shifted);
                double shift = shifted ? 0.5 : 0.;

                I ipix2 = ipix1 + nr - 1; // highest pixel number in the ring

                I ip_lo = ifloor(nr * inv_twopi * (phi - dphi) - shift) + 1;
                I ip_hi = ifloor(nr * inv_twopi * (phi + dphi) - shift);

                if (fct > 1) {
                    while ((ip_lo <= ip_hi) && check_pixel_ring
                            (ip_lo, nr, ipix1, fct, z0, phi, cosrsmall, cpix))
                        ++ip_lo;
                    while ((ip_hi > ip_lo) && check_pixel_ring
                            (ip_hi, nr, ipix1, fct, z0, phi, cosrsmall, cpix))
                        --ip_hi;
                }

                if (ip_lo <= ip_hi) {
                    if (ip_hi >= nr) {
                        ip_lo -= nr;
                        ip_hi -= nr;
                    }
                    if (ip_lo < 0) {
                        append(ipix1, ipix1 + ip_hi + 1, 3);
                        append(ipix1 + ip_lo + nr, ipix2 + 1, 4);
                    } else
                        append(ipix1 + ip_lo, ipix1 + ip_hi + 1, 5);

                }
            }
        }
        if ((rlat2 >= pi) && (irmax + 1 < 4 * nside_)) // south pole in the disk
        {
            I sp, rp;
            bool dummy;
            get_ring_info_small(irmax + 1, sp, rp, dummy);
            append(sp, npix_, 6);
        }
    } else {

        if (radius >= pi) // disk covers the whole sphere
        {
            append(0, npix_);
            return;
        }

        int oplus = 0;
        if (inclusive) {
            //  planck_assert ((I(1)<<(order_max-order_))>=fact, "invalid oversampling factor");
            //planck_assert ((fact&(fact-1))==0,"oversampling factor must be a power of 2");
            oplus = ilog2(fact);
        }
        int omax = order_ + oplus; // the order up to which we test

        //vec3 vptg(ptg);
        //arr<I> base(omax+1);
        arr<double> crpdr(omax + 1), crmdr(omax + 1);
        double cosrad = cos(radius);
        for (int o = 0; o <= omax; ++o) // prepare data at the required orders
        {
            // base[o].Set(o,NEST);
            double dr = max_pixrad(I(1) << o); // safety distance
            crpdr[o] = (radius + dr > pi) ? -1. : cos(radius + dr);
            crmdr[o] = (radius - dr < 0.) ? 1. : cos(radius - dr);
        }
        vector<pair<I, int> > stk; // stack for pixel numbers and their orders
        stk.reserve(12 + 3 * omax); // reserve maximum size to avoid reallocation
        for (int i = 0; i < 12; ++i) // insert the 12 base pixels in reverse order
            stk.push_back(make_pair(I(11 - i), 0));

        int stacktop = 0; // a place to save a stack position

        while (!stk.empty()) // as long as there are pixels on the stack
        {
            // cout<<stk.size()<<endl;
            // pop current pixel number and order from the stack
            I pix = stk.back().first;
            int o = stk.back().second;
            stk.pop_back();
            //edit by zqq
            double theta_r = 0, phi_r = 0;
            hpix_radec_degrees_to_thetaphi_radians(ra, dec, &theta_r, &phi_r);
            double x_r, y_r, z_r;
            thetaphi2xyz(theta_r, phi_r, x_r, y_r, z_r);
            //end zqq

            double z, phi;
            pix2zphi_NESTED(pix, z, phi, I(1) << o);
            // cosine of angular distance between pixel center and disk center
            //  cout<<"pix "<<pix<<" o "<<o<<" z "<<z<<" phi "<<phi<<endl;

            double cangdist = cosdist_zphi(z_r, phi_r, z, phi);

            if (cangdist > crpdr[o]) {
                int zone = (cangdist < cosrad) ? 1 : ((cangdist <= crmdr[o]) ? 2 : 3);

                check_pixel(o, order_, omax, zone, pix, stk, inclusive,
                            stacktop);
            }
        }
    }
}

static int ang2pix_ring_z_phi(int nside_, double z, double phi) {
    double za = fabs(z);
    double tt = fmodulo(phi, twopi) * inv_halfpi; /* in [0,4) */

    if (za <= twothird) /* Equatorial region */
    {
        double temp1 = nside_ * (0.5 + tt);
        double temp2 = nside_ * z * 0.75;
        int jp = (int) (temp1 - temp2); /* index of  ascending edge content */
        int jm = (int) (temp1 + temp2); /* index of descending edge content */

        /* ring number counted from z=2/3 */
        int ir = nside_ + 1 + jp - jm; /* in {1,2n+1} */
        int kshift = 1 - (ir & 1); /* kshift=1 if ir even, 0 otherwise */

        int ip = (jp + jm - nside_ + kshift + 1) / 2; /* in {0,4n-1} */
        ip = imodulo(ip, 4 * nside_);

        return nside_ * (nside_ - 1) * 2 + (ir - 1) * 4 * nside_ + ip;
    } else  /* North & South polar caps */
    {
        double tp = tt - (int) (tt);
        double tmp = nside_ * sqrt(3 * (1 - za));

        int jp = (int) (tp * tmp); /* increasing edge content index */
        int jm = (int) ((1.0 - tp) * tmp); /* decreasing edge content index */

        int ir = jp + jm + 1; /* ring number counted from the closest pole */
        int ip = (int) (tt * ir); /* in {0,4*ir-1} */
        ip = imodulo(ip, 4 * ir);

        if (z > 0)
            return 2 * ir * (ir - 1) + ip;
        else
            return 12 * nside_ * nside_ - 2 * ir * (ir + 1) + ip;
    }
}

long eq2pix_ring(int nside, double ra, double dec) {
    long t;
    double theta = 0, phi = 0;
    hpix_radec_degrees_to_thetaphi_radians(ra, dec, &theta, &phi);
    long ipix = ang2pix_ring_z_phi(nside, cos(theta), phi);
    return ipix;
}

vector<int> query_disc(int Nside, double ra, double dec, double range, bool r) {
    find3(ra, dec, range, Nside, 4, r);

    vector<int> v;
    std::copy(healpix_vector.begin(), healpix_vector.end(), std::back_inserter(v));

    return v;
}
/*
int main() {
//    double theta,phi;
//    pix2ang_ring_z_phi(2,0,&theta,&phi);
//    cout<<(theta-M_PI_2)/(-D2R)<<" "<<phi/D2R<<endl;
//    int totInMemory = nside2npix(2048);
    //  cout << totInMemory << endl;
    //  find2(64, 0, 0, 10);
//	for(int ra=0;ra<359;ra+=3){
//		for(int dec=90;dec>-90;dec-=3){
//			int s=0;
//
//	long t;
//    t =eq2pix(2,ra,dec);
//
////            double theta=0, phi=0;
////
////             hpix_radec_degrees_to_thetaphi_radians(ra, dec, &theta, &phi);
////
////            ang2pix_nest(1,theta,phi,&t);
//	if(t<10)
//	cout<<0;
//	cout<<t;
//		}
//		cout<<endl;
//	}
//cin>>phi;
    // void T_Healpix_Base::T_Healpix_Base(100, 1);
    int nside;
    double ra, dec, r;
    //  cin>>nside>>ra>>dec>>r;
    nside = 32;
    ra = 108.4;
    dec = -8.7;
    r = 10.0;
    find3(ra, dec, r * (pi / 180.0), nside, 4);
    sort(healpix_vector.begin(), healpix_vector.end());
    for (auto i = 0; i < healpix_vector.size(); i++) {
        cout << healpix_vector[i] << endl;
    }
    return 0;
}
*/
//int main(){
//    double theta, phi;
//    pix2ang_nest(2048, 45859308, &theta, &phi);
//    double dec2 = (theta - M_PI_2) / (-D2R);
//    double ra2 = phi / D2R;
//    cout<<ra2<<" "<<dec2<<endl;
//
//}
#undef I
#endif