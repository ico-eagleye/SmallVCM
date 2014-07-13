/*
 * Copyright (C) 2012, Tomas Davidovic (http://www.davidovic.cz)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * (The above is MIT License: http://en.wikipedia.org/wiki/MIT_License)
 */

#ifndef __VERTEXCM_HXX__
#define __VERTEXCM_HXX__
#include <vector>
#include <cmath>
#include <cassert>
#include "renderer.hxx"
#include "bsdf.hxx"
#include "rng.hxx"
#include "hashgrid.hxx"
#include "debug.h"


#define CONNECT_VERTICES 1
#define CONNECT_CAMERA 1
#define CONNECT_LIGHTS0 1 // getRadiance
#define CONNECT_LIGHTS1 1 // illuminate

////////////////////////////////////////////////////////////////////////////////
// A NOTE ON PATH MIS WEIGHT EVALUATION
////////////////////////////////////////////////////////////////////////////////
//
// We compute path MIS weights iteratively as we trace the light and eye
// sub-paths. We cache three floating points quantities at each sub-path vertex:
//
//   dVCM  dVC  dVM
//
// These quantities represent partial weights associated with the sub-path. When
// we connect or merge one vertex to another, we use these quantities to quickly
// evaluate the MIS weight for the full path we have constructed. This scheme is
// presented in the technical report
//
//   "Implementing Vertex Connection and Merging"
//   http://www.iliyan.com/publications/ImplementingVCM
//
// The MIS code in the VertexCM class references the corresponding equations in
// the report in the form
//
//   [tech. rep. (##)]
//
// where ## is the equation number. 
//

class VertexCM : public AbstractRenderer
{
    // The sole point of this structure is to make carrying around the ray baggage easier.
    struct SubPathState
    {
        Vec3f mOrigin;             // Path origin
        Vec3f mDirection;          // Where to go next
        Vec3f mThroughput;         // Path throughput
        uint  mPathLength    : 30; // Number of path segments, including this
        uint  mIsFiniteLight :  1; // Just generate by finite light
        uint  mSpecularPath  :  1; // All scattering events so far were specular

        float dVCM; // MIS quantity used for vertex connection and merging
        float dVC;  // MIS quantity used for vertex connection
        float dVM;  // MIS quantity used for vertex merging
    };

    // Path vertex, used for merging and connection
    template<bool tFromLight>
    struct PathVertex
    {
        Vec3f mHitpoint;   // Position of the vertex
        Vec3f mThroughput; // Path throughput (including emission)
        uint  mPathLength; // Number of segments between source and vertex

        // Stores all required local information, including incoming direction.
        BSDF<tFromLight> mBsdf;

        float dVCM; // MIS quantity used for vertex connection and merging
        float dVC;  // MIS quantity used for vertex connection
        float dVM;  // MIS quantity used for vertex merging

        // Used by HashGrid
        const Vec3f& GetPosition() const
        {
            return mHitpoint;
        }
    };

    typedef PathVertex<false> CameraVertex;
    typedef PathVertex<true>  LightVertex;

    typedef BSDF<false>       CameraBSDF;
    typedef BSDF<true>        LightBSDF;

    // Range query used for PPM, BPT, and VCM. When HashGrid finds a vertex
    // within range -- Process() is called and vertex
    // merging is performed. BSDF of the camera vertex is used.
    class RangeQuery
    {
    public:

        RangeQuery(
            const VertexCM     &aVertexCM,
            const Vec3f        &aCameraPosition,
            const CameraBSDF   &aCameraBsdf,
            const SubPathState &aCameraState
        ) : 
            mVertexCM(aVertexCM),
            mCameraPosition(aCameraPosition),
            mCameraBsdf(aCameraBsdf),
            mCameraState(aCameraState),
            mContrib(0)
        {}

        const Vec3f& GetPosition() const { return mCameraPosition; }

        const Vec3f& GetContrib() const { return mContrib; }

        void Process(const LightVertex& aLightVertex)
        {
            // Reject if full path length below/above min/max path length
            if((aLightVertex.mPathLength + mCameraState.mPathLength > mVertexCM.mMaxPathLength) ||
               (aLightVertex.mPathLength + mCameraState.mPathLength < mVertexCM.mMinPathLength))
                 return;

            // Retrieve light incoming direction in world coordinates
            const Vec3f lightDirection = aLightVertex.mBsdf.WorldDirFix();

            float cosCamera, cameraBsdfDirPdfW, cameraBsdfRevPdfW;
            const Vec3f cameraBsdfFactor = mCameraBsdf.Evaluate(
                mVertexCM.mScene, lightDirection, cosCamera, &cameraBsdfDirPdfW,
                &cameraBsdfRevPdfW);

            if(cameraBsdfFactor.IsZero())
                return;

            cameraBsdfDirPdfW *= mCameraBsdf.ContinuationProb();

            // Even though this is pdf from camera BSDF, the continuation probability
            // must come from light BSDF, because that would govern it if light path
            // actually continued
            cameraBsdfRevPdfW *= aLightVertex.mBsdf.ContinuationProb();		// vmarz!: review

            // Partial light sub-path MIS weight [tech. rep. (38)]
            const float wLight = aLightVertex.dVCM * mVertexCM.mMisVcWeightFactor +
                aLightVertex.dVM * mVertexCM.Mis(cameraBsdfDirPdfW);
            // vmarz?: why cameraDirPdf when formula has s-2, e.g. pdf of the vertex before merged vertex?
            // because pdf is reverse in relation to light paht Y, e.g. it's pdf of (s-2)<-(s-1)<-(s)

            // Partial eye sub-path MIS weight [tech. rep. (39)]
            const float wCamera = mCameraState.dVCM * mVertexCM.mMisVcWeightFactor +
                mCameraState.dVM * mVertexCM.Mis(cameraBsdfRevPdfW);
            // vmarz: check reasoning above why cameraRevPdf used

            // Full path MIS weight [tech. rep. (37)]. No MIS for PPM
            const float misWeight = mVertexCM.mPpm ?
                1.f :
                1.f / (wLight + 1.f + wCamera);

            mContrib += misWeight * cameraBsdfFactor * aLightVertex.mThroughput;
        }

    private:

        const VertexCM     &mVertexCM;
        const Vec3f        &mCameraPosition;
        const CameraBSDF   &mCameraBsdf;
        const SubPathState &mCameraState;
        Vec3f              mContrib;
    };

public:

    enum AlgorithmType
    {
        // light vertices contribute to camera,
        // No MIS weights (dVCM, dVM, dVC all ignored)
        kLightTrace = 0,

        // Camera and light vertices merged on first non-specular surface from camera.
        // Cannot handle mixed specular + non-specular materials.
        // No MIS weights (dVCM, dVM, dVC all ignored)
        kPpm,

        // Camera and light vertices merged on along full path.
        // dVCM and dVM used for MIS
        kBpm,

        // Standard bidirectional path tracing
        // dVCM and dVC used for MIS
        kBpt,

        // Vertex connection and mering
        // dVCM, dVM, and dVC used for MIS
        kVcm
    };

public:

    VertexCM(
        const Scene&  aScene,
        AlgorithmType aAlgorithm,
        const float   aRadiusFactor,
        const float   aRadiusAlpha,
        int           aSeed = 1234
    ) :
        AbstractRenderer(aScene),
        mRng(aSeed),
        mLightTraceOnly(false),
        mUseVC(false),
        mUseVM(false),
        mPpm(false)
    {
        switch(aAlgorithm)
        {
        case kLightTrace:
            mLightTraceOnly = true;
            break;
        case kPpm:
            mPpm   = true;
            mUseVM = true;
            break;
        case kBpm:
            mUseVM = true;
            break;
        case kBpt:
            mUseVC = true;
            break;
        case kVcm:
            mUseVC = true;
            mUseVM = true;
            break;
        default:
            printf("Unknown algorithm requested\n");
            break;
        }

        if(mPpm)
        {
            // We will check the scene to make sure it does not contain mixed
            // specular and non-specular materials
            for(int i = 0; i < mScene.GetMaterialCount(); ++i)
            {
                const Material &mat = mScene.GetMaterial(i);

                const bool hasNonSpecular =
                    (mat.mDiffuseReflectance.Max() > 0) ||
                    (mat.mPhongReflectance.Max() > 0);

                const bool hasSpecular =
                    (mat.mMirrorReflectance.Max() > 0) ||
                    (mat.mIOR > 0);

                if(hasNonSpecular && hasSpecular)
                {
                    printf(
                        "*WARNING* Our PPM implementation cannot handle materials mixing\n"
                        "Specular and NonSpecular BSDFs. The extension would be\n"
                        "fairly straightforward. In SampleScattering for camera sub-paths\n"
                        "limit the considered events to Specular only.\n"
                        "Merging will use non-specular components, scattering will be specular.\n"
                        "If there is no specular component, the ray will terminate.\n\n");

                    printf("We are now switching from *PPM* to *BPM*, which can handle the scene\n\n");

                    mPpm = false;
                    break;
                }
            }
        }

        mBaseRadius  = aRadiusFactor * mScene.mSceneSphere.mSceneRadius;
        mRadiusAlpha = aRadiusAlpha;
    }

    virtual void RunIteration(int aIteration)
    {
        if (DEBUG_THREAD_ID == omp_get_thread_num())
            dbgPrintf("\n\n\nITERATION: %d ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n", aIteration);
        // While we have the same number of pixels (camera paths)
        // and light paths, we do keep them separate for clarity reasons
        const int resX = int(mScene.mCamera.mResolution.x);
        const int resY = int(mScene.mCamera.mResolution.y);
        const int pathCount = resX * resY;
        mScreenPixelCount = float(resX * resY);
        mLightSubPathCount = float(resX * resY);

        // Setup our radius, 1st iteration has aIteration == 0, thus offset
        float radius = mBaseRadius;
        radius /= std::pow(float(aIteration + 1), 0.5f * (1 - mRadiusAlpha));
        // Purely for numeric stability
        radius = std::max(radius, 1e-7f);
        const float radiusSqr = Sqr(radius);

        // Factor used to normalize vertex merging contribution.
        // We divide the summed up energy by disk radius and number of light paths
        mVmNormalization = 1.f / (radiusSqr * PI_F * mLightSubPathCount);
        // vmarz: 1/(PI*r*r) in mVmNormalization coming from P_vm [tech. rep. (10)]
        // vmarz?: why mLightSubPathCount? because of N_vm in [tech. rep. (11)]

        // MIS weight constant [tech. rep. (20)], with n_VC = 1 and n_VM = mLightPathCount
        const float etaVCM = (PI_F * radiusSqr) * mLightSubPathCount; // / n_VC =1 ;
        mMisVmWeightFactor = mUseVM ? Mis(etaVCM)       : 0.f;
        mMisVcWeightFactor = mUseVC ? Mis(1.f / etaVCM) : 0.f;
        if (DEBUG_THREAD_ID == omp_get_thread_num())
            dbgPrintf("         etaVCM % 14f misVcWeightFac % 14f         radius % 14f\n", etaVCM, mMisVcWeightFactor, radius);

        // Clear path ends, nothing ends anywhere
        mPathEnds.resize(pathCount);
        memset(&mPathEnds[0], 0, mPathEnds.size() * sizeof(int));

        // Remove all light vertices and reserve space for some
        mLightVertices.reserve(pathCount);
        mLightVertices.clear();

        //////////////////////////////////////////////////////////////////////////
        // Generate light paths
        //////////////////////////////////////////////////////////////////////////
        for(int pathIdx = 0; pathIdx < pathCount; pathIdx++)
        {
            DBG_PRINTFI(&pathIdx, "\n\nLIGHT PASS: %d ----------------------------------------------------------------------------------\n", aIteration);
            SubPathState lightState;
            GenerateLightSample(lightState, &pathIdx);

            //////////////////////////////////////////////////////////////////////////
            // Trace light path
            for(;; ++lightState.mPathLength)
            {
                // Offset ray origin instead of setting tmin due to numeric
                // issues in ray-sphere intersection. The isect.dist has to be
                // extended by this EPS_RAY after hit point is determined
                Ray ray(lightState.mOrigin + lightState.mDirection * EPS_RAY,
                    lightState.mDirection, 0);
                Isect isect(1e36f);

                if(!mScene.Intersect(ray, isect))
                {
                    DBG_PRINTFI(&pathIdx, "MISS \n", lightState.mPathLength);
                    break;
                }

                const Vec3f hitPoint = ray.org + ray.dir * isect.dist;
                isect.dist += EPS_RAY;
                DBG_PRINTFI(&pathIdx, "\nHIT-L %d - %d\n", aIteration, lightState.mPathLength);
                DBG_PRINTFI(&pathIdx, "       hitPoint % 14f % 14f % 14f \n", hitPoint.x, hitPoint.y, hitPoint.z);

                LightBSDF bsdf(ray, isect, mScene);
                if(!bsdf.IsValid())
                {
                    break;
                    DBG_PRINTFI(&pathIdx, " Hit BSDF INVALID \n");
                }

                DBG_PRINTFI(&pathIdx, "   directionFix % 14f % 14f % 14f \n", bsdf.mLocalDirFix.x, bsdf.mLocalDirFix.y, bsdf.mLocalDirFix.z);
                DBG_PRINTFI(&pathIdx, "Light Hit MIS update \n");
                DBG_PRINTFI(&pathIdx, "            dVC % 14f            dVM % 14f           dVCM % 14f \n", 
                    lightState.dVC, lightState.dVM, lightState.dVCM);

                // Update the MIS quantities before storing them at the vertex.
                // These updates follow the initialization in GenerateLightSample() or
                // SampleScattering(), and together implement equations [tech. rep. (31)-(33)]
                // or [tech. rep. (34)-(36)], respectively.
                {
                    // Infinite lights use MIS handled via solid angle integration,
                    // so do not divide by the distance for such lights [tech. rep. Section 5.1]
                    if(lightState.mPathLength > 1 || lightState.mIsFiniteLight == 1)
                        lightState.dVCM *= Mis(Sqr(isect.dist)); 
                    // vmarz: from g in 1/p1 (or 1/pi)
                    //        for dVC and dVM sqr(dist) terms cancel out, see explanation in SampleScattering()
                    DBG_PRINTFI(&pathIdx, "         U dVCM % 14f    *= sqrDist  % 14f \n", lightState.dVCM, Mis(Sqr(isect.dist)) );

                    // vmarz: from g in p1 (or 1/pi)
                    lightState.dVCM /= Mis(std::abs(bsdf.CosThetaFix())); // vmarz?: why abs here?
                    lightState.dVC  /= Mis(std::abs(bsdf.CosThetaFix())); // not really needed since bsdf initialization
                    lightState.dVM  /= Mis(std::abs(bsdf.CosThetaFix())); // rejects rays when abs(mLocalDirFix.z) < EPS_COSINE by not setting materialID and
                }                                                         // causing bsdf.IsValid() to false which is checked few lines above for continuation
                DBG_PRINTFI(&pathIdx, "   /cosThetaFix % 14f \n", bsdf.CosThetaFix() );
                DBG_PRINTFI(&pathIdx, "          U dVC % 14f          U dVM % 14f         U dVCM % 14f \n", 
                    lightState.dVC, lightState.dVM, lightState.dVCM);

                // Store vertex, unless BSDF is purely specular, which prevents
                // vertex connections and merging
                if(!bsdf.IsDelta() && (mUseVC || mUseVM))
                {
                    LightVertex lightVertex;
                    lightVertex.mHitpoint   = hitPoint;
                    lightVertex.mThroughput = lightState.mThroughput;
                    lightVertex.mPathLength = lightState.mPathLength;
                    lightVertex.mBsdf       = bsdf;

                    lightVertex.dVCM = lightState.dVCM;
                    lightVertex.dVC  = lightState.dVC;
                    lightVertex.dVM  = lightState.dVM;

                    mLightVertices.push_back(lightVertex);
                }

                // vmarz: mandatory?
                // Connect to camera, unless BSDF is purely specular
#if CONNECT_CAMERA
                if(!bsdf.IsDelta() && (mUseVC || mLightTraceOnly))
                {
                    if(lightState.mPathLength + 1 >= mMinPathLength)
                        ConnectToCamera(lightState, hitPoint, bsdf, &pathIdx);
                }
#endif
                // Terminate if the path would become too long after scattering
                if(lightState.mPathLength + 2 > mMaxPathLength) // vmarz?: why +2 ?
                {
                    DBG_PRINTFI(&pathIdx, "MAX PATH LENGTH %d \n", mMaxPathLength);
                    break;
                }

                // Continue random walk
                if(!SampleScattering(bsdf, hitPoint, lightState, &pathIdx))
                    break;
            }

            mPathEnds[pathIdx] = (int)mLightVertices.size();
        }

        //////////////////////////////////////////////////////////////////////////
        // Build hash grid
        //////////////////////////////////////////////////////////////////////////

        // Only build grid when merging (VCM, BPM, and PPM)
        if(mUseVM)
        {
            // The number of cells is somewhat arbitrary, but seems to work ok
            mHashGrid.Reserve(pathCount);
            mHashGrid.Build(mLightVertices, radius);
        }

        //////////////////////////////////////////////////////////////////////////
        // Generate camera paths
        //////////////////////////////////////////////////////////////////////////

        // Unless rendering with traditional light tracing
        for(int pathIdx = 0; (pathIdx < pathCount) && (!mLightTraceOnly); ++pathIdx)
        {
            DBG_PRINTFI(&pathIdx, "\n\nCAMERA PASS: %d ----------------------------------------------------------------------------------\n");
            SubPathState cameraState;
            const Vec2f screenSample = GenerateCameraSample(pathIdx, cameraState);
            Vec3f color(0);

            //////////////////////////////////////////////////////////////////////
            // Trace camera path
            for(;; ++cameraState.mPathLength)
            {
                // Offset ray origin instead of setting tmin due to numeric
                // issues in ray-sphere intersection. The isect.dist has to be
                // extended by this EPS_RAY after hit point is determined
                Ray ray(cameraState.mOrigin + cameraState.mDirection * EPS_RAY,
                    cameraState.mDirection, 0);

                Isect isect(1e36f);

                // Get radiance from environment
                if(!mScene.Intersect(ray, isect))
                {
                    if(mScene.GetBackground() != NULL)
                    {
                        if(cameraState.mPathLength >= mMinPathLength)
                        {
                            Vec3f contrib = cameraState.mThroughput * GetLightRadiance(mScene.GetBackground(), cameraState,
                                Vec3f(0), ray.dir);
                            color += contrib;
                        }
                    }

                    break;
                }

                const Vec3f hitPoint = ray.org + ray.dir * isect.dist;
                isect.dist += EPS_RAY;
                DBG_PRINTFI(&pathIdx, "\nHIT-C %d - %d \n", aIteration, cameraState.mPathLength);
                DBG_PRINTFI(&pathIdx, "       hitPoint % 14f % 14f % 14f \n", hitPoint.x, hitPoint.y, hitPoint.z);

                CameraBSDF bsdf(ray, isect, mScene);
                if(!bsdf.IsValid())
                {
                    DBG_PRINTFI(&pathIdx, " Hit BSDF INVALID \n");
                    break;
                }

                DBG_PRINTFI(&pathIdx, "   directionFix % 14f % 14f % 14f \n", bsdf.mLocalDirFix.x, bsdf.mLocalDirFix.y, bsdf.mLocalDirFix.z);
                DBG_PRINTFI(&pathIdx, "Camera Hit MIS update \n");
                DBG_PRINTFI(&pathIdx, "            dVC % 14f            dVM % 14f           dVCM % 14f \n", 
                    cameraState.dVC, cameraState.dVM, cameraState.dVCM);
                // Update the MIS quantities, following the initialization in
                // GenerateLightSample() or SampleScattering(). Implement equations
                // [tech. rep. (31)-(33)] or [tech. rep. (34)-(36)], respectively.
                {
                    cameraState.dVCM *= Mis(Sqr(isect.dist));
                    DBG_PRINTFI(&pathIdx, "         U dVCM % 14f    *= sqrDist  % 14f \n", cameraState.dVCM, Mis(Sqr(isect.dist)) );
                    cameraState.dVCM /= Mis(std::abs(bsdf.CosThetaFix()));
                    cameraState.dVC  /= Mis(std::abs(bsdf.CosThetaFix()));
                    cameraState.dVM  /= Mis(std::abs(bsdf.CosThetaFix()));
                    DBG_PRINTFI(&pathIdx, "   /cosThetaFix % 14f \n", bsdf.CosThetaFix() );
                    DBG_PRINTFI(&pathIdx, "          U dVC % 14f          U dVM % 14f         U dVCM % 14f \n", 
                        cameraState.dVC, cameraState.dVM, cameraState.dVCM);
                }

#if CONNECT_LIGHTS0
                // Light source has been hit; terminate afterwards, since
                // our light sources do not have reflective properties
                if(isect.lightID >= 0)
                {
                    const AbstractLight *light = mScene.GetLightPtr(isect.lightID);
                
                    if(cameraState.mPathLength >= mMinPathLength)
                    {
                        color += cameraState.mThroughput *
                            GetLightRadiance(light, cameraState, hitPoint, ray.dir);
                    }
                    
                    break;
                }
#endif

                // Terminate if eye sub-path is too long for connections or merging
                if(cameraState.mPathLength >= mMaxPathLength)
                {
                    DBG_PRINTFI(&pathIdx, "MAX PATH LENGTH %d \n", mMaxPathLength);
                    break;
                }

#if CONNECT_LIGHTS1
                ////////////////////////////////////////////////////////////////
                // Vertex connection: Connect to a light source
                if(!bsdf.IsDelta() && mUseVC)
                {
                    if(cameraState.mPathLength + 1>= mMinPathLength)
                    {
                        Vec3f contrib = cameraState.mThroughput *
                            DirectIllumination(cameraState, hitPoint, bsdf, &pathIdx);
                        DBG_PRINTFI(&pathIdx, "camera.Throughpt % 14f % 14f % 14f \n", cameraState.mThroughput.x, cameraState.mThroughput.y, cameraState.mThroughput.z);
                        DBG_PRINTFI(&pathIdx, "contri*througput % 14f % 14f % 14f \n", contrib.x, contrib.y, contrib.z);
                        color += contrib;
                    }
                }
#endif

#if CONNECT_VERTICES
                ////////////////////////////////////////////////////////////////
                // Vertex connection: Connect to light vertices
                if(!bsdf.IsDelta() && mUseVC)
                {
                    // For VC, each light sub-path is assigned to a particular eye
                    // sub-path, as in traditional BPT. It is also possible to
                    // connect to vertices from any light path, but MIS should
                    // be revisited.
                    
                    // vmarz?: doesn't it imply need to revisit MIS also if using Light Vertex Cache?
                    // I guess just means need to be computed correctly, e.g. cases of delayed computation of some factors
                    const Vec2i range(
                        (pathIdx == 0) ? 0 : mPathEnds[pathIdx-1],
                        mPathEnds[pathIdx]);

                    for(int i = range.x; i < range.y; i++)
                    {
                        const LightVertex &lightVertex = mLightVertices[i];

                        if(lightVertex.mPathLength + 1 +
                           cameraState.mPathLength < mMinPathLength)
                            continue;

                        // Light vertices are stored in increasing path length
                        // order; once we go above the max path length, we can
                        // skip the rest
                        if(lightVertex.mPathLength + 1 +
                           cameraState.mPathLength > mMaxPathLength)
                            break;

                        Vec3f connectContrib = cameraState.mThroughput * lightVertex.mThroughput *
                            ConnectVertices(lightVertex, bsdf, hitPoint, cameraState, &pathIdx);
                        DBG_PRINTFI(&pathIdx, "contri*througput % 14f % 14f % 14f \n", connectContrib.x, connectContrib.y, connectContrib.z);
                        DBG_PRINTFI(&pathIdx, "camera.Throughpt % 14f % 14f % 14f \n", cameraState.mThroughput.x, cameraState.mThroughput.y, cameraState.mThroughput.z);
                        DBG_PRINTFI(&pathIdx, "vertex.Throughpt % 14f % 14f % 14f \n", lightVertex.mThroughput.x, lightVertex.mThroughput.y, lightVertex.mThroughput.z);
                        color += connectContrib;
                    }
                }
#endif

                ////////////////////////////////////////////////////////////////
                // Vertex merging: Merge with light vertices
                if(!bsdf.IsDelta() && mUseVM)
                {
                    RangeQuery query(*this, hitPoint, bsdf, cameraState);
                    mHashGrid.Process(mLightVertices, query);
                    color += cameraState.mThroughput * mVmNormalization * query.GetContrib();
                    // vmarz: path pdfs already divided into camera Throughput (and light Throughput in RangeQuery.Process)
                    // vmarz: 1/(PI*r*r) in mVmNormalization coming from P_vm [tech. rep. (10)]

                    // PPM merges only at the first non-specular surface from camera
                    if(mPpm) break;
                }

                if(!SampleScattering(bsdf, hitPoint, cameraState))
                    break;

                //if (IS_DEBUG_IDX(&pathIdx))
                //{
                //    color = Vec3f(10000.f);
                //}
            }

            mFramebuffer.AddColor(screenSample, color);

            int a = 1; // vmarz dummy to stay in scope
        }

        mIterations++;
    }

private:

    // Mis power, we use balance heuristic 
    float Mis(float aPdf) const     // vmarz: No Mis2() as in pathtracer.cxx because whole MIS weight partitioned.
    {                               // Here only chosen power heuristic factor is potentially applied
        //return std::pow(aPdf, /*power*/);
        return aPdf;
    }

    //////////////////////////////////////////////////////////////////////////
    // Camera tracing methods
    //////////////////////////////////////////////////////////////////////////

    // Generates new camera sample given a pixel index
    Vec2f GenerateCameraSample(
        const int    aPixelIndex,
        SubPathState &oCameraState,
        const int    *idx = NULL)
    {
        DBG_PRINTFI(idx, "GenerateCameraSample(): \n");
        const Camera &camera = mScene.mCamera;
        const int resX = int(camera.mResolution.x);
        const int resY = int(camera.mResolution.y);

        // Determine pixel (x, y)
        const int x = aPixelIndex % resX;
        const int y = aPixelIndex / resX;

        // Jitter pixel position
        const Vec2f sample = Vec2f(float(x), float(y)) + mRng.GetVec2f();

        // Generate ray
        const Ray primaryRay = camera.GenerateRay(sample);

        // Compute pdf conversion factor from area on image plane to solid angle on ray
        const float cosAtCamera = Dot(camera.mForward, primaryRay.dir);
        const float imagePointToCameraDist = camera.mImagePlaneDist / cosAtCamera;
        const float imageToSolidAngleFactor = Sqr(imagePointToCameraDist) / cosAtCamera;

        // We put the virtual image plane at such a distance from the camera origin
        // that the pixel area is one and thus the image plane sampling pdf is 1.
        // The solid angle ray pdf is then equal to the conversion factor from
        // image plane area density to ray solid angle density
        const float cameraPdfW = imageToSolidAngleFactor /* * areaPdf */; // vmarz: areaPdf = 1/area = 1/1 = 1
                                                                          //        cameraPdf = areaPdf * p1 = p1_ro * g1 [g1 added after tracing]
        oCameraState.mOrigin       = primaryRay.org;                      //        p0_connect = areaPdf = 1 
        oCameraState.mDirection    = primaryRay.dir;
        oCameraState.mThroughput   = Vec3f(1);

        oCameraState.mPathLength   = 1;
        oCameraState.mSpecularPath = 1;

        // Eye sub-path MIS quantities. Implements [tech. rep. (31)-(33)] partially.
        // The evaluation is completed after tracing the camera ray in the eye sub-path loop.
        oCameraState.dVCM = Mis( /* p0_connect * */ mLightSubPathCount / cameraPdfW); // vmarz: dVCM = (p0connect/p0trace)*(nLightSamples/p1)
        oCameraState.dVC  = 0;
        oCameraState.dVM  = 0;

        DBG_PRINTFI(idx, "         origin % 14f % 14f % 14f \n", oCameraState.mOrigin.x, oCameraState.mOrigin.y, oCameraState.mOrigin.z);
        DBG_PRINTFI(idx, "      direction % 14f % 14f % 14f \n", oCameraState.mDirection.x, oCameraState.mDirection.y, oCameraState.mDirection.z);
        DBG_PRINTFI(idx, "dVCM=subPathCnt % 14d / cameraPdfW % 14f \n", mLightSubPathCount, cameraPdfW);
        DBG_PRINTFI(idx, "            dVC % 14f            dVM % 14f           dVCM % 14f \n", oCameraState.dVC, oCameraState.dVM, oCameraState.dVCM);

        return sample;
    }

    // Returns the radiance of a light source when hit by a random ray,
    // multiplied by MIS weight. Can be used for both Background and Area lights.
    //
    // For Background lights:
    //    Has to be called BEFORE updating the MIS quantities.
    //    Value of aHitpoint is irrelevant (passing Vec3f(0))
    //
    // For Area lights:
    //    Has to be called AFTER updating the MIS quantities.
    Vec3f GetLightRadiance(
        const AbstractLight *aLight,
        const SubPathState  &aCameraState,
        const Vec3f         &aHitpoint,
        const Vec3f         &aRayDirection) const
    {
        // We sample lights uniformly
        const int   lightCount    = mScene.GetLightCount();
        const float lightPickProb = 1.f / lightCount;

        float directPdfA, emissionPdfW;
        const Vec3f radiance = aLight->GetRadiance(mScene.mSceneSphere,
            aRayDirection, aHitpoint, &directPdfA, &emissionPdfW);

        if(radiance.IsZero())
            return Vec3f(0);

        // If we see light source directly from camera, no weighting is required
        if(aCameraState.mPathLength == 1)
            return radiance;

        // When using only vertex merging, we want purely specular paths
        // to give radiance (cannot get it otherwise). Rest is handled
        // by merging and we should return 0.
        if(mUseVM && !mUseVC)
            return aCameraState.mSpecularPath ? radiance : Vec3f(0);

        directPdfA   *= lightPickProb; // vmarz: p0connect in tech. rep
        emissionPdfW *= lightPickProb; // vmarz: p0trace in tech. rep

        // Partial eye sub-path MIS weight [tech. rep. (43)].
        // If the last hit was specular, then dVCM == 0.
        const float wCamera = Mis(directPdfA) * aCameraState.dVCM +
            Mis(emissionPdfW) * aCameraState.dVC;

        // Partial light sub-path weight is 0 [tech. rep. (42)].

        // Full path MIS weight [tech. rep. (37)].
        const float misWeight = 1.f / (1.f + wCamera);
        
        return misWeight * radiance;
    }

    // Connects camera vertex to randomly chosen light point.
    // Returns emitted radiance multiplied by path MIS weight.
    // Has to be called AFTER updating the MIS quantities.
    Vec3f DirectIllumination(
        const SubPathState &aCameraState,
        const Vec3f        &aHitpoint,
        const CameraBSDF   &aBsdf,
        const int          *idx = NULL)
    {
        DBG_PRINTFI(idx, "\DirectIllumination(): \n");
        // We sample lights uniformly
        const int   lightCount    = mScene.GetLightCount();
        const float lightPickProb = 1.f / lightCount;

        const int   lightID       = int(mRng.GetFloat() * lightCount);
        const Vec2f rndPosSamples = mRng.GetVec2f();

        const AbstractLight *light = mScene.GetLightPtr(lightID);

        Vec3f directionToLight;
        float distance;
        float directPdfW, emissionPdfW, cosAtLight;
        const Vec3f radiance = light->Illuminate(mScene.mSceneSphere, aHitpoint,
            rndPosSamples, directionToLight, distance, directPdfW,
            &emissionPdfW, &cosAtLight);

        if (IS_DEBUG_IDX(idx))
        {
            int a = 1;
        }
        
        DBG_PRINTFI(idx, "       radiance % 14f % 14f % 14f \n", radiance.x, radiance.y, radiance.z);
        DBG_PRINTFI(idx, "directionToLigt % 14f % 14f % 14f      distance % 14f \n", directionToLight.x, directionToLight.y, directionToLight.z, distance);

        // If radiance == 0, other values are undefined, so have to early exit
        if(radiance.IsZero())
        {
            DBG_PRINTFI(idx, "radiance ZERO \n");
            return Vec3f(0);
        }

        DBG_PRINTFI(idx, "     cosAtLight % 14f   emissionPdfW % 14f     directPdfW % 14f \n", cosAtLight, emissionPdfW, directPdfW);

        float bsdfDirPdfW, bsdfRevPdfW, cosToLight;
        const Vec3f bsdfFactor = aBsdf.Evaluate(mScene,
            directionToLight, cosToLight, &bsdfDirPdfW, &bsdfRevPdfW);

        DBG_PRINTFI(idx, "     cosToLight % 14f \n", cosToLight);
        DBG_PRINTFI(idx, "     bsdfFactor % 14f % 14f % 14f \n", bsdfFactor.x, bsdfFactor.y, bsdfFactor.z);

        if(bsdfFactor.IsZero())
        {
            DBG_PRINTFI(idx, "bsdfFactor ZERO \n");
            return Vec3f(0);
        }

        const float continuationProbability = aBsdf.ContinuationProb();
        
        // If the light is delta light, we can never hit it
        // by BSDF sampling, so the probability of this path is 0
        bsdfDirPdfW *= light->IsDelta() ? 0.f : continuationProbability;

        bsdfRevPdfW *= continuationProbability;

        // Partial light sub-path MIS weight [tech. rep. (44)].
        // Note that wLight is a ratio of area pdfs. But since both are on the
        // light source, their distance^2 and cosine terms cancel out.
        // Therefore we can write wLight as a ratio of solid angle pdfs,
        // both expressed w.r.t. the same shading point.
        const float wLight = Mis(bsdfDirPdfW / (lightPickProb * directPdfW));

        DBG_PRINTFI(idx, "         wLight =    bsdfDirPdfW / ( lightPickProb *     directPdfW )\n");
        DBG_PRINTFI(idx, " % 14f = % 14f / (% 14f * % 14f ) \n", wLight, bsdfDirPdfW, lightPickProb, directPdfW);

        // Partial eye sub-path MIS weight [tech. rep. (45)].
        //
        // In front of the sum in the parenthesis we have Mis(ratio), where
        //    ratio = emissionPdfA / directPdfA,
        // with emissionPdfA being the product of the pdfs for choosing the
        // point on the light source and sampling the outgoing direction.
        // What we are given by the light source instead are emissionPdfW
        // and directPdfW. Converting to area pdfs and plugging into ratio:
        //    emissionPdfA = emissionPdfW * cosToLight / dist^2
        //    directPdfA   = directPdfW * cosAtLight / dist^2
        //    ratio = (emissionPdfW * cosToLight / dist^2) / (directPdfW * cosAtLight / dist^2)
        //    ratio = (emissionPdfW * cosToLight) / (directPdfW * cosAtLight)
        //
        // Also note that both emissionPdfW and directPdfW should be
        // multiplied by lightPickProb, so it cancels out.
        const float wCamera = Mis(emissionPdfW * cosToLight / (directPdfW * cosAtLight)) * (
            mMisVmWeightFactor + aCameraState.dVCM + aCameraState.dVC * Mis(bsdfRevPdfW));

        DBG_PRINTFI(idx, "        wCamera = (  emissionPdfW *     cosToLight / (    directPdfW *     cosAtLight )) * ( vmWeightFactor +    camera.dVCM +     camera.dVC *    bsdfRevPdfW ) \n");
        DBG_PRINTFI(idx, " % 14f = (% 14f * % 14f / (% 14f * % 14f)) * (% 14f + % 14f + % 14f + % 14f) \n",
            wCamera, emissionPdfW, cosToLight, directPdfW, cosAtLight, mMisVmWeightFactor, aCameraState.dVCM, aCameraState.dVC, bsdfRevPdfW);

        // Full path MIS weight [tech. rep. (37)]
        const float misWeight = 1.f / (wLight + 1.f + wCamera);
        DBG_PRINTFI(idx, "      misWeight % 14f \n", misWeight);


        // vmarz: radiance not scaled by cosAtLight, also not in Illuminate function
        Vec3f contrib = (cosToLight / (lightPickProb * directPdfW)) * (radiance * bsdfFactor);

        DBG_PRINTFI(idx, "unweigh contrb % 14f % 14f % 14f \n", contrib.x, contrib.y, contrib.z);
        DBG_PRINTFI(idx, "unweigh contrb = (     cosToLight / ( lightPickProb *     directPdfW)) * (      radiance *     bsdfFactor ) \n");
        DBG_PRINTFI(idx, "unweigh contrb = ( % 14f / (% 14f * % 14f )) * (% 14f * % 14f ) \n",
            cosToLight, lightPickProb, directPdfW, radiance, bsdfFactor);
        contrib *= misWeight;
        DBG_PRINTFI(idx, " weight contrib % 14f % 14f % 14f \n", contrib.x, contrib.y, contrib.z);

        if(contrib.IsZero() || mScene.Occluded(aHitpoint, directionToLight, distance))
            return Vec3f(0);

        return contrib;
    }

    // Connects an eye and a light vertex. Result multiplied by MIS weight, but
    // not multiplied by vertex throughputs. Has to be called AFTER updating MIS
    // constants. 'direction' is FROM eye TO light vertex.
    Vec3f ConnectVertices(
        const LightVertex  &aLightVertex,
        const CameraBSDF   &aCameraBsdf,
        const Vec3f        &aCameraHitpoint,
        const SubPathState &aCameraState,
        const int          *idx = NULL) const
    {
        DBG_PRINTFI(idx, "\nConnectVertices(): \n");
        // Get the connection
        Vec3f direction   = aLightVertex.mHitpoint - aCameraHitpoint;
        const float dist2 = direction.LenSqr();
        float  distance   = std::sqrt(dist2);
        direction        /= distance;

        // Evaluate BSDF at camera vertex
        float cosCamera, cameraBsdfDirPdfW, cameraBsdfRevPdfW;
        const Vec3f cameraBsdfFactor = aCameraBsdf.Evaluate(
            mScene, direction, cosCamera, &cameraBsdfDirPdfW,
            &cameraBsdfRevPdfW);

        DBG_PRINTFI(idx, "  connect point % 14f % 14f % 14f        pathLen %d\n",
            aLightVertex.mHitpoint.x, aLightVertex.mHitpoint.y, aLightVertex.mHitpoint.z, aLightVertex.mPathLength);
        DBG_PRINTFI(idx, "      direction % 14f % 14f % 14f      distance % 14f \n", direction.x, direction.y, direction.z, distance);
        DBG_PRINTFI(idx, "cameraBsdfFactr % 14f % 14f % 14f \n", cameraBsdfFactor.x, cameraBsdfFactor.y, cameraBsdfFactor.z);
        DBG_PRINTFI(idx, "      cosCamera % 14f camBsdfDirPdfW % 14f camBsdfRevPdfW % 14f \n", cosCamera, cameraBsdfDirPdfW, cameraBsdfRevPdfW);

        if(cameraBsdfFactor.IsZero())
        {
            DBG_PRINTFI(idx, "cameraBsdfFactor ZERO \n");
            return Vec3f(0);
        }

        // Camera continuation probability (for Russian roulette)
        const float cameraCont = aCameraBsdf.ContinuationProb();
        cameraBsdfDirPdfW *= cameraCont;
        cameraBsdfRevPdfW *= cameraCont;
        DBG_PRINTFI(idx, "                     RR scaled camBsdfDirPdfW % 14f camBsdfRevPdfW % 14f \n", cameraBsdfDirPdfW, cameraBsdfRevPdfW);

        // Evaluate BSDF at light vertex
        float cosLight, lightBsdfDirPdfW, lightBsdfRevPdfW;
        const Vec3f lightBsdfFactor = aLightVertex.mBsdf.Evaluate(
            mScene, -direction, cosLight, &lightBsdfDirPdfW,
            &lightBsdfRevPdfW);

        DBG_PRINTFI(idx, "lightBsdfFactor % 14f % 14f % 14f \n", lightBsdfFactor.x, lightBsdfFactor.y, lightBsdfFactor.z);
        DBG_PRINTFI(idx, "       cosLight % 14f lgtBsdfDirPdfW % 14f lgtBsdfRevPdfW % 14f \n", cosLight, lightBsdfDirPdfW, lightBsdfRevPdfW);

        if(lightBsdfFactor.IsZero())
        {
            DBG_PRINTFI(idx, "lightBsdfFactor ZERO \n");
            return Vec3f(0);
        }

        // Light continuation probability (for Russian roulette)
        const float lightCont = aLightVertex.mBsdf.ContinuationProb();
        lightBsdfDirPdfW *= lightCont;
        lightBsdfRevPdfW *= lightCont;
        DBG_PRINTFI(idx, "                     RR scaled lgtBsdfDirPdfW % 14f lgtBsdfRevPdfW % 14f \n", lightBsdfDirPdfW, lightBsdfRevPdfW);

        // Compute geometry term
        const float geometryTerm = cosLight * cosCamera / dist2;
        DBG_PRINTFI(idx, "  geometryTerm % 14f \n", geometryTerm);
        if(geometryTerm < 0)
        {
            DBG_PRINTFI(idx, "geometryTerm < ZERO \n");
            return Vec3f(0);
        }

        // Convert pdfs to area pdf
        const float cameraBsdfDirPdfA = PdfWtoA(cameraBsdfDirPdfW, distance, cosLight);
        const float lightBsdfDirPdfA  = PdfWtoA(lightBsdfDirPdfW,  distance, cosCamera);
        DBG_PRINTFI(idx, " camBsdfDirPdfA = (camBsdfDirPdfW *       cosLight) / sqr (      distance) \n");
        DBG_PRINTFI(idx, " % 14f = (% 14f * % 14f) / sqr (% 14f) \n", cameraBsdfDirPdfA, cameraBsdfDirPdfW, cosLight, distance);
        DBG_PRINTFI(idx, " lgtBsdfDirPdfA = (lgtBsdfDirPdfW *      cosCamera) / sqr (      distance) \n");
        DBG_PRINTFI(idx, " % 14f = (% 14f * % 14f) / sqr (% 14f) \n", lightBsdfDirPdfA, lightBsdfDirPdfW, cosCamera, distance);

        // Partial light sub-path MIS weight [tech. rep. (40)]
        const float wLight = Mis(cameraBsdfDirPdfA) * (
            mMisVmWeightFactor + aLightVertex.dVCM + aLightVertex.dVC * Mis(lightBsdfRevPdfW));
        // vmarz: lightBsdfRevPdfW is Reverse with respect to light path, e.g. in eye path progression 
        // direction (note same arrow dirs in formula)
        // note (40) and (41) uses light subpath Y and camera subpath z
        DBG_PRINTFI(idx, "         wLight = camBsdfDirPdfA * (VmWeightFactor +     light.dVCM +      light.dVC * lgtBsdfRevPdfW) \n");
        DBG_PRINTFI(idx, " % 14f = % 14f * (% 14f + % 14e + % 14e * % 14f) \n", 
            wLight, cameraBsdfDirPdfA, mMisVmWeightFactor, aLightVertex.dVCM, aLightVertex.dVC, lightBsdfRevPdfW);

        // Partial eye sub-path MIS weight [tech. rep. (41)]
        const float wCamera = Mis(lightBsdfDirPdfA) * (
            mMisVmWeightFactor + aCameraState.dVCM + aCameraState.dVC * Mis(cameraBsdfRevPdfW));
        DBG_PRINTFI(idx, "        wCamera = lgtBsdfDirPdfA * (VmWeightFactor +    camera.dVCM +     camera.dVC * camBsdfRevPdfW) \n");
        DBG_PRINTFI(idx, " % 14f = % 14f * (% 14f + % 14e + % 14e * % 14f) \n", 
            wLight, lightBsdfDirPdfA, mMisVmWeightFactor, aCameraState.dVCM, aCameraState.dVC, cameraBsdfRevPdfW);

        // Full path MIS weight [tech. rep. (37)]
        const float misWeight = 1.f / (wLight + 1.f + wCamera);

        Vec3f contrib = misWeight * geometryTerm * cameraBsdfFactor * lightBsdfFactor; 
        // vmarz: 1) Where is divide by path pdf? A: it is predivided into throughput at every scattering
        //        2) Should divide contrib also by by lightVertexPickPdf (numConnect/numVertices) in case of use of LightVertexCache ?
        //           In this case numVertices = numLightSubpathVertices and numConnect=numLightSubpathVertices, therefore cancel out ?
        DBG_PRINTFI(idx, "      misWeight % 14f \n", misWeight);
        DBG_PRINTFI(idx, "unweigh contrib % 14f % 14f % 14f \n", contrib.x, contrib.y, contrib.z);
        DBG_PRINTFI(idx, "unweigh contrib = geometryTerm * cameraBsdfFactor * lightBsdfFactor\n");
        contrib *= misWeight;
        DBG_PRINTFI(idx, " weight contrib % 14f % 14f % 14f \n", contrib.x, contrib.y, contrib.z);

        if(contrib.IsZero() || mScene.Occluded(aCameraHitpoint, direction, distance))
        {
            DBG_PRINTFI(idx, "OCCLUDED \n");
            return Vec3f(0);
        }

        return contrib;
    }

    //////////////////////////////////////////////////////////////////////////
    // Light tracing methods
    //////////////////////////////////////////////////////////////////////////

    // Samples light emission
    void GenerateLightSample(SubPathState &oLightState, int *idx = NULL)
    {
        DBG_PRINTFI(idx, "GenerateLightSample(): \n");
        // We sample lights uniformly
        const int   lightCount    = mScene.GetLightCount();
        const float lightPickProb = 1.f / lightCount;

        const int   lightID       = int(mRng.GetFloat() * lightCount);
        const Vec2f rndDirSamples = mRng.GetVec2f();
        const Vec2f rndPosSamples = mRng.GetVec2f();

        const AbstractLight *light = mScene.GetLightPtr(lightID);

        float emissionPdfW, directPdfW, cosLight;
        oLightState.mThroughput = light->Emit(mScene.mSceneSphere, rndDirSamples, rndPosSamples,
            oLightState.mOrigin, oLightState.mDirection,
            emissionPdfW, &directPdfW, &cosLight);
        // vmarz?: AreaLight->Emit sets directPdfW to oDirectPdfA = mInvArea; not really pdf w.r.t solid angle?
        DBG_PRINTFI(idx, "         origin % 14f % 14f % 14f \n", oLightState.mOrigin.x, oLightState.mOrigin.y, oLightState.mOrigin.z);
        DBG_PRINTFI(idx, "      direction % 14f % 14f % 14f \n", oLightState.mDirection.x, oLightState.mDirection.y, oLightState.mDirection.z);
        DBG_PRINTFI(idx, "       emission % 14f % 14f % 14f \n", oLightState.mThroughput.x, oLightState.mThroughput.y, oLightState.mThroughput.z);
        DBG_PRINTFI(idx, "       cosLight % 14f   emissionPdfW % 14f     directPdfW % 14f \n", cosLight, emissionPdfW, directPdfW);

        emissionPdfW *= lightPickProb;
        directPdfW   *= lightPickProb;
        DBG_PRINTFI(idx, "       scaled w light pic prob - emissionPdfW % 14f   U directPdfW % 14f \n", emissionPdfW, directPdfW);

        oLightState.mThroughput    /= emissionPdfW;
        oLightState.mPathLength    = 1;
        oLightState.mIsFiniteLight = light->IsFinite() ? 1 : 0;

        // Light sub-path MIS quantities. Implements [tech. rep. (31)-(33)] partially.
        // The evaluation is completed after tracing the emission ray in the light sub-path loop.
        // Delta lights are handled as well [tech. rep. (48)-(50)].
        {
            oLightState.dVCM = Mis(directPdfW / emissionPdfW); // vmarz: dVCM_1 = p0_connect / ( p0_trace * p1 )   [connect/trace potentially different points sampling techniques]
                                                               //        directPdfW = p0_connect = areaSamplePdf * lightPickPdf
                                                               //        emissionPdfW = p0_trace * p1 
            if(!light->IsDelta())                              //           p0_trace = areaSamplePdf * lightPickPdf
            {                                                  //           p1 = cosineSamplePdf * g1 = (cos / Pi) * g1 [g1 added after tracing]
                const float usedCosLight = light->IsFinite() ? cosLight : 1.f;
                oLightState.dVC = Mis(usedCosLight / emissionPdfW);  // vmarz: dVC_1 = _g0 / ( p0_trace * p1 )   [connect/trace potentially different points sampling techniques]
            }                                                        //        usedCosLight is part of _g0   [ _g - reverse pdf conversion factor!, uses outgoing cos not incident at next vertex]
            else                                                     //           [sqr(dist) from _g0 added after tracing]
            {                                                        //        emissionPdfW = p0_trace * p1 
                oLightState.dVC = 0.f;
            }

            oLightState.dVM = oLightState.dVC * mMisVcWeightFactor;
            // vmarz: dVM_1 = dVC_1 / etaVCM
            //        [sqr(dist) from _g0 added after tracing]
        }
        
        DBG_PRINTFI(idx, "     throughput % 14f % 14f % 14f \n", oLightState.mThroughput.x, oLightState.mThroughput.y, oLightState.mThroughput.z);
        DBG_PRINTFI(idx, "            dVC % 14f            dVM % 14f           dVCM % 14f \n", oLightState.dVC, oLightState.dVM, oLightState.dVCM);
    }

    // Computes contribution of light sample to camera by splatting is onto the
    // framebuffer. Multiplies by throughput (obviously, as nothing is returned).
    void ConnectToCamera(
        const SubPathState &aLightState,
        const Vec3f        &aHitpoint,
        const LightBSDF    &aBsdf,
        int *idx = NULL )
    {
        const Camera &camera    = mScene.mCamera;
        Vec3f directionToCamera = camera.mPosition - aHitpoint;

        // Check point is in front of camera
        if(Dot(camera.mForward, -directionToCamera) <= 0.f)
            return;

        // Check it projects to the screen (and where)
        const Vec2f imagePos = camera.WorldToRaster(aHitpoint);
        if(!camera.CheckRaster(imagePos))
            return;

        // Compute distance and normalize direction to camera
        const float distEye2 = directionToCamera.LenSqr();
        const float distance = std::sqrt(distEye2);
        directionToCamera   /= distance;

        // Get the BSDF
        float cosToCamera, bsdfDirPdfW, bsdfRevPdfW;
        const Vec3f bsdfFactor = aBsdf.Evaluate(mScene,
            directionToCamera, cosToCamera, &bsdfDirPdfW, &bsdfRevPdfW);

        int dbgPixel  = ( DEBUG_PIX && IS_DEBUG_PIX(imagePos));
        int dbgLaunch = (!DEBUG_PIX && IS_DEBUG_IDX(idx));
        int dbgCond = dbgPixel || dbgLaunch;

        DBG_PRINTFC(dbgCond, "ConnectToCamera():    debugPixel: %d  pixelPos %d %d   launchIdx x %u y %u \n",
            dbgPixel, int(imagePos.x), int(imagePos.y), IDX_X(*idx), IDX_Y(*idx) );
        DBG_PRINTFC(dbgCond, "    dirToCamera % 14f % 14f % 14f      distance % 14f \n", directionToCamera.x, directionToCamera.y, directionToCamera.z, distance);
        DBG_PRINTFC(dbgCond, "     bsdfFactor % 14f % 14f % 14f \n", bsdfFactor.x, bsdfFactor.y, bsdfFactor.z);
        DBG_PRINTFC(dbgCond, "    cosToCamera % 14f    bsdfDirPdfW % 14f    bsdfRevPdfW % 14f \n", cosToCamera, bsdfDirPdfW, bsdfRevPdfW);

        if(bsdfFactor.IsZero())
        {
            DBG_PRINTFC(dbgCond, "bsdfFactor ZERO \n ");
            return;
        }

        bsdfRevPdfW *= aBsdf.ContinuationProb();

        // Compute pdf conversion factor from image plane area to surface area
        const float cosAtCamera = Dot(camera.mForward, -directionToCamera);
        const float imagePointToCameraDist = camera.mImagePlaneDist / cosAtCamera;
        const float imageToSolidAngleFactor = Sqr(imagePointToCameraDist) / cosAtCamera;
        const float imageToSurfaceFactor = imageToSolidAngleFactor * std::abs(cosToCamera) / Sqr(distance);

        DBG_PRINTFC(dbgCond, " imgToPtCamDist =  imagePlaneDist /    cosAtCamera) \n");
        DBG_PRINTFC(dbgCond, " % 14f = % 14f  * % 14f \n", imagePointToCameraDist, camera.mImagePlaneDist, cosAtCamera);
        DBG_PRINTFC(dbgCond, " imgToSolAngFac =  sqr(imgToPtCamDist) /    cosAtCamera) \n");
        DBG_PRINTFC(dbgCond, " % 14f = sqr(% 14f) / % 14f \n", imageToSolidAngleFactor, imagePointToCameraDist, cosAtCamera);
        DBG_PRINTFC(dbgCond, "  imgToSurfFact = imgSolAngleFac * abs(   cosToCamera)) / sqr(      distance) \n");
        DBG_PRINTFC(dbgCond, " % 14f = % 14f * abs(% 14f) / sqr(% 14f) \n", imageToSurfaceFactor, imageToSolidAngleFactor, cosToCamera, distance);

        // We put the virtual image plane at such a distance from the camera origin
        // that the pixel area is one and thus the image plane sampling pdf is 1.
        // The area pdf of aHitpoint as sampled from the camera is then equal to
        // the conversion factor from image plane area density to surface area density
        const float cameraPdfA = imageToSurfaceFactor; // * 1.f 
        DBG_PRINTFC(dbgCond, "    cameraPdfA % 14f \n", cameraPdfA);

        // Partial light sub-path weight [tech. rep. (46)]. Note the division by
        // mLightPathCount, which is the number of samples this technique uses.
        // This division also appears a few lines below in the framebuffer accumulation.
        const float wLight = Mis(cameraPdfA / mLightSubPathCount) * (
            mMisVmWeightFactor + aLightState.dVCM + aLightState.dVC * Mis(bsdfRevPdfW));
        DBG_PRINTFC(dbgCond, "         wLight = (    cameraPdfA / lightPathCount) * (vmWeightFactor +     light.dVCM +      light.dVC *    bsdfRevPdfW) \n");
        DBG_PRINTFC(dbgCond, " % 14f = (% 14f / % 14f) * (% 14f + % 14e + % 14e * % 14f) \n", 
            wLight, cameraPdfA, mLightSubPathCount, mMisVmWeightFactor, aLightState.dVCM, aLightState.dVC, bsdfRevPdfW);

        // Partial eye sub-path weight is 0 [tech. rep. (47)]

        // Full path MIS weight [tech. rep. (37)]. No MIS for traditional light tracing.
        const float misWeight = mLightTraceOnly ? 1.f : (1.f / (wLight + 1.f));

        const float surfaceToImageFactor = 1.f / imageToSurfaceFactor;

        DBG_PRINTFC(dbgCond, "       misWeight % 14f         wLight % 14f\n", misWeight, wLight);
        DBG_PRINTFC(dbgCond, "  srfToImgFactor % 14f 1/imgSurfFactr % 14f lghtPathCount % 14f \n",  surfaceToImageFactor, imageToSurfaceFactor, mLightSubPathCount);

        // We divide the contribution by surfaceToImageFactor to convert the (already
        // divided) pdf from surface area to image plane area, w.r.t. which the
        // pixel integral is actually defined. We also divide by the number of samples
        // this technique makes, which is equal to the number of light sub-paths
        // vmarz: "to convert the (already divided) pdf from surface area.."
        //         already divided where? throughput?
        Vec3f contrib =  aLightState.mThroughput * bsdfFactor /
            (mLightSubPathCount * surfaceToImageFactor ) ;
        DBG_PRINTFC(dbgCond, " light.throughpt % 14f % 14f % 14f           depth % 14u \n", aLightState.mThroughput.x, aLightState.mThroughput.y, aLightState.mThroughput.z, aLightState.mPathLength)
        DBG_PRINTFC(dbgCond, " unweigh contrib % 14f % 14f % 14f \n", contrib.x, contrib.y, contrib.z);
        DBG_PRINTFC(dbgCond, " unweigh contrib = light.throughpt *     bsdfFactor / (lightPathCount * srfToImgFactor \n");
        contrib *= misWeight;
        DBG_PRINTFC(dbgCond, " weight contrib % 14f % 14f % 14f \n\n", contrib.x, contrib.y, contrib.z)

        if(!contrib.IsZero())
        {
            if(mScene.Occluded(aHitpoint, directionToCamera, distance))
                return;

            mFramebuffer.AddColor(imagePos, contrib);
        }
    }

    // Samples a scattering direction camera/light sample according to BSDF.
    // Returns false for termination
    template<bool tLightSample>
    bool SampleScattering(
        const BSDF<tLightSample> &aBsdf,
        const Vec3f              &aHitPoint,
        SubPathState             &aoState,
        int                      *idx = NULL)
    {
        DBG_PRINTFI(idx, "SampleScattering(): \n");

        // x,y for direction, z for component. No rescaling happens
        Vec3f rndTriplet  = mRng.GetVec3f();
        float bsdfDirPdfW, cosThetaOut;
        uint  sampledEvent;

        Vec3f bsdfFactor = aBsdf.Sample(mScene, rndTriplet, aoState.mDirection,
            bsdfDirPdfW, cosThetaOut, &sampledEvent);

        DBG_PRINTFI(idx, "      direction % 14f % 14f % 14f \n", aoState.mDirection.x, aoState.mDirection.y, aoState.mDirection.z);
        DBG_PRINTFI(idx, "     bsdfFactor % 14f % 14f % 14f \n", bsdfFactor.x, bsdfFactor.y, bsdfFactor.z);
        DBG_PRINTFI(idx, "    cosThetaOut % 14f \n", cosThetaOut);

        if(bsdfFactor.IsZero())
        {
            return false;
            DBG_PRINTFI(idx, "bsdfFactor Zero \n");
        }

        // If we sampled specular event, then the reverse probability
        // cannot be evaluated, but we know it is exactly the same as
        // forward probability, so just set it. If non-specular event happened,
        // we evaluate the pdf
        float bsdfRevPdfW = bsdfDirPdfW;
        if((sampledEvent & LightBSDF::kSpecular) == 0)
            bsdfRevPdfW = aBsdf.Pdf(mScene, aoState.mDirection, true);

        // Russian roulette
        const float contProb = aBsdf.ContinuationProb();
        const float rrSample = mRng.GetFloat();
        DBG_PRINTFI(idx, "       contProb % 14f             RR % 14f \n", contProb, rrSample);
        if(rrSample > contProb)
        {
            DBG_PRINTFI(idx, "RR STOP \n");
            return false;
        }

        DBG_PRINTFI(idx, "    bsdfDirPdfW % 14f    bsdfRevPdfW % 14f \n", bsdfDirPdfW, bsdfRevPdfW);
        bsdfDirPdfW *= contProb;
        bsdfRevPdfW *= contProb;

        DBG_PRINTFI(idx, "  U bsdfDirPdfW % 14f  U bsdfRevPdfW % 14f\n", bsdfDirPdfW, bsdfRevPdfW);
        DBG_PRINTFI(idx, "            dVC % 14f            dVM % 14f           dVCM % 14f \n",
            aoState.dVC, aoState.dVM, aoState.dVCM);

        const float dVC = aoState.dVC;
        const float dVM = aoState.dVM;
        const float dVCM = aoState.dVCM;

        // Sub-path MIS quantities for the next vertex. Only partial - the
        // evaluation is completed when the actual hit point is known,
        // i.e. after tracing the ray, in the sub-path loop.
        if(sampledEvent & LightBSDF::kSpecular)
        {
            // Specular scattering case [tech. rep. (53)-(55)] (partially, as noted above)
            aoState.dVCM = 0.f;
            //aoState.dVC *= Mis(cosThetaOut / bsdfDirPdfW) * Mis(bsdfRevPdfW);
            //aoState.dVM *= Mis(cosThetaOut / bsdfDirPdfW) * Mis(bsdfRevPdfW);
            assert(bsdfDirPdfW == bsdfRevPdfW);
            aoState.dVC *= Mis(cosThetaOut);
            aoState.dVM *= Mis(cosThetaOut);

            aoState.mSpecularPath &= 1;
        }
        else
        {
            // Implements [tech. rep. (34)-(36)] (partially, as noted above)
            aoState.dVC = Mis(cosThetaOut / bsdfDirPdfW) * ( // vmarz: dVC = (g_i-1 / pi) * (etaVCM + dVCM_i-1 + _p_ro_i-2 * dVC_i-1)
                aoState.dVC * Mis(bsdfRevPdfW) +             //        cosThetaOut part of g_i-1  [ _g reverse pdf conversion!, uses outgoing cosTheta]
                aoState.dVCM + mMisVmWeightFactor);          //          !! sqr(dist) terms for _g_i-1 and gi of pi are the same and cancel out, hence NOT scaled after tracing]
                                                             //        pi = bsdfDirPdfW * g1
            aoState.dVM = Mis(cosThetaOut / bsdfDirPdfW) * ( //        bsdfDirPdfW = _p_ro_i    [part of pi]
                aoState.dVM * Mis(bsdfRevPdfW) +             //        bsdfRevPdfW = _p_ro_i-2
                aoState.dVCM * mMisVcWeightFactor + 1.f);    // 
                                                             //        dVM = (g_i-1 / pi) * (1 + dVCM_i-1/etaVCM + _p_ro_i-2 * dVM_i-1)
            aoState.dVCM = Mis(1.f / bsdfDirPdfW);           //        cosThetaOut part of g_i-1 [_g reverse pdf conversion!, uses outgoing cosTheta]
                                                             //          !! sqr(dist) terms for _g_i-1 and gi of pi are the same and cancel out, hence NOT scaled after tracing]
            aoState.mSpecularPath &= 0;                      //
        }                                                    //        dVCM = 1 / pi
                                                             //        pi = bsdfDirPdfW * g1 = _p_ro_i * g1 [only for dVCM sqe(dist) terms do not cancel out and are added after tracing]
        aoState.mOrigin  = aHitPoint;
        aoState.mThroughput *= bsdfFactor * (cosThetaOut / bsdfDirPdfW);
        // vmarz: GEOMETRY TERM NOTE
        //        cosAtHitPoint / sqr(dist) not known at this point until tracking, but still important to note they won't be actually used
        //        since they cancel out due to division of sampled dir solid angle pdf by solid angle to area pdf conversion factor.
        //        See PBR 765

        DBG_PRINTFI(idx, "     throughput % 14f = bsdfFactor * (cosThetaOut / bsdfDirPdfW) \n", aoState.mThroughput.x, aoState.mThroughput.y, aoState.mThroughput.z);
        DBG_PRINTFI(idx, "          U dVC = (   cosThetaOut /    bsdfDirPdfW) * (           dVC *    bsdfRevPdfW +           dVCM + VmWeightFactor) \n");
        DBG_PRINTFI(idx, " % 14f = (% 14f / % 14f) * (% 14e * % 14f + % 14e + % 14f) \n", 
            aoState.dVC, cosThetaOut, bsdfDirPdfW, dVC, bsdfRevPdfW, aoState.dVCM, mMisVmWeightFactor);

        DBG_PRINTFI(idx, "          U dVM = (   cosThetaOut /    bsdfDirPdfW) * (           dVM *    bsdfRevPdfW +           dVCM + VcWeightFactor + 1) \n");
        DBG_PRINTFI(idx, " % 14f = (% 14f / % 14f) * (% 14e * % 14f + % 14e + % 14f + 1) \n", 
            aoState.dVM, cosThetaOut, bsdfDirPdfW, dVM, bsdfRevPdfW, aoState.dVCM, mMisVcWeightFactor);
        DBG_PRINTFI(idx, "         U dVCM = (1 /    bsdfDirPdfW) \n");
        DBG_PRINTFI(idx, " % 14f = (1 / %14f) \n",  dVCM, bsdfDirPdfW);
        DBG_PRINTFI(idx, "          U dVC % 14f          U dVM % 14f         U dVCM % 14f \n", aoState.dVC, aoState.dVM, aoState.dVCM);
        return true;
    }

private:

    bool  mUseVM;             // Vertex merging (of some form) is used
    bool  mUseVC;             // Vertex connection (BPT) is used
    bool  mLightTraceOnly;    // Do only light tracing
    bool  mPpm;               // Do PPM, same terminates camera after first merge

    float mRadiusAlpha;       // Radius reduction rate parameter
    float mBaseRadius;        // Initial merging radius
    float mMisVmWeightFactor; // Weight of vertex merging (used in VC)
    float mMisVcWeightFactor; // Weight of vertex connection (used in VM)
    float mScreenPixelCount;  // Number of pixels
    float mLightSubPathCount; // Number of light sub-paths
    float mVmNormalization;   // 1 / (Pi * radius^2 * light_path_count)

    std::vector<LightVertex> mLightVertices; //!< Stored light vertices

    // For light path belonging to pixel index [x] it stores
    // where it's light vertices end (begin is at [x-1])
    std::vector<int> mPathEnds;
    HashGrid         mHashGrid;

    Rng              mRng;
};

#endif //__VERTEXCM_HXX__
