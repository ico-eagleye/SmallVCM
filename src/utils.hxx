/*
 * This is published under Apache 2.0
 */

#ifndef __UTILS_HXX__
#define __UTILS_HXX__

#include <vector>
#include <cmath>
#include "math.hxx"

// sRGB luminance
float Luminance(const Vec3f& aRGB)
{
    return 0.212671f * aRGB.x +
        0.715160f * aRGB.y +
        0.072169f * aRGB.z;
}

float FresnelDielectric(float aCosInc, float mIOR)
{
    if(mIOR < 0) return 1.f;

    float etaIncOverEtaTrans;

    if(aCosInc < 0.f)
    {
        aCosInc = -aCosInc;
        etaIncOverEtaTrans = mIOR;
    }
    else
    {
        etaIncOverEtaTrans = 1.f / mIOR;
    }

    const float sinTrans2 = Sqr(etaIncOverEtaTrans) * (1.f - Sqr(aCosInc));
    const float cosTrans = std::sqrt(std::max(0.f, 1.f - sinTrans2));

    const float term1 = etaIncOverEtaTrans * cosTrans;
    const float rParallel =
        (aCosInc - term1) / (aCosInc + term1);

    const float term2 = etaIncOverEtaTrans * aCosInc;
    const float rPerpendicular =
        (term2 - cosTrans) / (term2 + cosTrans);

    return 0.5f * (Sqr(rParallel) + Sqr(rPerpendicular));
}

// reflect vector through (0,0,1)
Vec3f reflect001(const Vec3f& aVector)
{
    return Vec3f(aVector.x, aVector.y, -aVector.z);
}

//////////////////////////////////////////////////////////////////////////
// Cosine lobe hemisphere sampling

Vec3f SamplePowerCosHemisphereW(
    const Vec2f  &aSamples,
    const float  aPower,
    float        *oPdfW)
{
    const float term1 = 2.f * PI_F * aSamples.x;
    const float term2 = std::pow(aSamples.y, 1.f / (aPower + 1.f));
    const float term3 = std::sqrt(1.f - term2 * term2);

    if(oPdfW)
    {
        *oPdfW = (aPower + 1.f) * std::pow(term2, aPower) * (0.5f * INV_PI_F);
    }

    return Vec3f(
        std::cos(term1) * term3,
        std::sin(term1) * term3,
        term2);
}

float EvalPowerCosHemispherePdfW(
    const Vec3f  &aNormal,
    const Vec3f  &aDirection,
    const float  aPower)
{
    const float cosTheta = std::max(0.f, Dot(aNormal, aDirection));

    return (aPower + 1.f) * std::pow(cosTheta, aPower) * (INV_PI_F * 0.5f);
}

//////////////////////////////////////////////////////////////////////////
/// Sample direction in the upper hemisphere with cosine-proportional pdf
/** The returned PDF is with respect to solid angle measure */
Vec3f SampleCosHemisphereW(
    const Vec2f  &aSamples,
    float        *oPdfW)
{
    const float term1 = 2.f * PI_F * aSamples.x;
    const float term2 = std::sqrt(1.f - aSamples.y);

    const Vec3f ret(
        std::cos(term1) * term2,
        std::sin(term1) * term2,
        std::sqrt(aSamples.y));

    if(oPdfW)
    {
        *oPdfW = ret.z * INV_PI_F;
    }

    return ret;
}

float EvalCosHemispherePdfW(
    const Vec3f  &aNormal,
    const Vec3f  &aDirection)
{
    return std::max(0.f, Dot(aNormal, aDirection)) * INV_PI_F;
}

// Sample Triangle
// returns barycentric coordinates
Vec2f SampleUniformTriangle(
    const Vec2f &aSamples)
{
    const float term = std::sqrt(aSamples.x);

    return Vec2f(1.f - term, aSamples.y * term);
}

//////////////////////////////////////////////////////////////////////////
// Pdf utils using the new naming scheme
// W - solind angle
// A - area
float FactorWtoA_dist2(
    const float aDist2,
    const float aCosThere)
{
    return std::abs(aCosThere) / aDist2;
}

float FactorWtoA(
    const float aDist,
    const float aCosThere)
{
    return FactorWtoA_dist2(Sqr(aDist), aCosThere);
}


float PdfWtoA_dist2(
    const float aPdfW,
    const float aDist2,
    const float aCosThere)
{
    return aPdfW * FactorWtoA_dist2(aDist2, aCosThere);
}

float PdfWtoA(
    const float aPdfW,
    const float aDist,
    const float aCosThere)
{
    return PdfWtoA_dist2(aPdfW, Sqr(aDist), aCosThere);
}

float FactorAtoW_dist2(
    const float aDist2,
    const float aCosThere)
{
    return aDist2 / std::abs(aCosThere);
}

float FactorAtoW(
    const float aDist,
    const float aCosThere)
{
    return FactorAtoW_dist2(Sqr(aDist), aCosThere);
}


float PdfAtoW_dist2(
    const float aPdfA,
    const float aDist2,
    const float aCosThere)
{
    return aPdfA * FactorAtoW_dist2(aDist2, aCosThere);
}

float PdfAtoW(
    const float aPdfA,
    const float aDist,
    const float aCosThere)
{
    return PdfAtoW_dist2(aPdfA, Sqr(aDist), aCosThere);
}

// Mis power (1 for balance heuristic)
float Mis(float aPdf) { return aPdf; }

// Mis weight for 2 pdfs
float Mis2(float aSamplePdf, float aOtherPdf)
{
    return Mis(aSamplePdf) / (Mis(aSamplePdf) + Mis(aOtherPdf));
}

#endif //__UTILS_HXX__