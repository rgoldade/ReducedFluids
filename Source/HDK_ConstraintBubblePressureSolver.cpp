#include "HDK_ConstraintBubblePressureSolver.h"

#include "tbb/tbb.h"

#include <GA/GA_PageHandle.h>
#include <GA/GA_Iterator.h>
#include <GA/GA_SplittableRange.h>

#include <GU/GU_Detail.h>

#include <PRM/PRM_Include.h>

#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_FieldUtils.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_VolumetricConnectedComponentBuilder.h>

#include <UT/UT_DSOVersion.h>
#include <UT/UT_ParallelUtil.h>
#include <UT/UT_PerfMonAutoEvent.h>

#include <UT/UT_StopWatch.h>

void
initializeSIM(void *)
{
   IMPLEMENT_DATAFACTORY(HDK_ConstraintBubblePressureSolver);
}

// Standard constructor, note that BaseClass was crated by the
// DECLARE_DATAFACTORY and provides an easy way to chain through
// the class hierarchy.
HDK_ConstraintBubblePressureSolver::HDK_ConstraintBubblePressureSolver(const SIM_DataFactory *factory)
    : BaseClass(factory)
{
}

HDK_ConstraintBubblePressureSolver::~HDK_ConstraintBubblePressureSolver()
{
}

const SIM_DopDescription*
HDK_ConstraintBubblePressureSolver::getDopDescription()
{
    static PRM_Name	theSurfaceName(GAS_NAME_SURFACE, "Surface Field");
    static PRM_Default 	theSurfaceDefault(0,"surface");

    static PRM_Name	theSurfacePressureName("surfacepressure","Surface Pressure Field");
    static PRM_Default	theSurfacePressureDefault(0, "surfacepressure");

    static PRM_Name	theVelocityName(GAS_NAME_VELOCITY, "Velocity Field");
    static PRM_Default	theVelocityDefault(0, "vel");

    static PRM_Name    theSolidSurfaceName(GAS_NAME_COLLISION, "Solid Field");
    static PRM_Default theSolidSurfaceDefault(0, "collision");

    static PRM_Name	theSolidVelocityName(GAS_NAME_COLLISIONVELOCITY, "Solid Velocity Field");
    static PRM_Default  theSolidVelocityDefault(0, "collisionvel");

    static PRM_Name 	theCutCellWeightsName("cutCellWeights", "Cut-cell Weights Field");
    static PRM_Default  theCutCellWeightsDefault(0, "collisionweights");

    static PRM_Name 	thePressureName(GAS_NAME_PRESSURE, "Pressure");
    static PRM_Default 	thePressureNameDefault(0, "pressure");

    static PRM_Name    theUseOldPressureName("useOldPressure", "Use old pressure as an initial guess");

    static PRM_Name 	theDensityName(GAS_NAME_DENSITY, "Liquid Density Field");
    static PRM_Default 	theDensityNameDefault(0, "massdensity");

    static PRM_Name 	theValidFacesName("validFaces", "Valid Faces Field");

    static PRM_Name 	theToleranceName(SIM_NAME_TOLERANCE, "Error Tolerance");
    static PRM_Default 	theToleranceDefault(1e-5);

    static PRM_Name 	theMaxIterations("maxIterations", "Max Solver Iterations");
    static PRM_Default 	theMaxIterationsDefault(2500);

    static PRM_Name	theUseMGPreconditionerName("useMGPreconditioner", "Use Multigrid Preconditioner");

    static PRM_Name	theUseOpenBoundariesName("useOpenBoundaries", "Use Open Boundaries");

    static PRM_Name	theTrackBubbleIDsName("trackBubbleIDs", "Track Bubble IDs");

    static PRM_Name	theTrackBubbleVolumesName("trackBubbleVolumes", "Track Bubble Volumes");
    static PRM_Conditional theTrackBubbleVolumesDisable("{ trackBubbleIDs == 0 }");

    static PRM_Name	theCorrectVolumeDriftName("correctVolumeDrift", "Correct Volume Drift");
    static PRM_Conditional theCorrectVolumeDriftDisable("{ trackBubbleIDs == 0 }{ trackBubbleVolumes == 0 }");

    static PRM_Name	theGeometryName(GAS_NAME_GEOMETRY, "Bubble ID Tracking Geometry");
    static PRM_Default  theGeometryNameDefault(0, "Geometry");

    static PRM_Name theBubbleVolumeGeometryName("bubbleVolumeGeometry", "Bubble Volume Tracking Geometry");
    static PRM_Default  theBubbleVolumeGeometryNameDefault(0, "BubbleVolume");

    static PRM_Template	theTemplates[] =
    {
        PRM_Template(PRM_STRING, 1, &theSurfaceName, &theSurfaceDefault),

    	PRM_Template(PRM_STRING, 1, &theSurfacePressureName, &theSurfacePressureDefault),

    	PRM_Template(PRM_STRING, 1, &theVelocityName, &theVelocityDefault),

    	PRM_Template(PRM_STRING, 1, &theSolidSurfaceName, &theSolidSurfaceDefault),

    	PRM_Template(PRM_STRING, 1, &theSolidVelocityName, &theSolidVelocityDefault),

    	PRM_Template(PRM_STRING, 1, &theCutCellWeightsName, &theCutCellWeightsDefault),

        PRM_Template(PRM_STRING, 1, &thePressureName, &thePressureNameDefault),

	PRM_Template(PRM_TOGGLE, 1, &theUseOldPressureName, PRMoneDefaults),

        PRM_Template(PRM_STRING, 1, &theDensityName, &theDensityNameDefault),

        PRM_Template(PRM_STRING, 1, &theValidFacesName),

    	PRM_Template(PRM_FLT, 1, &theToleranceName, &theToleranceDefault),

    	PRM_Template(PRM_FLT, 1, &theMaxIterations, &theMaxIterationsDefault),

	PRM_Template(PRM_TOGGLE, 1, &theUseMGPreconditionerName, PRMoneDefaults),

	PRM_Template(PRM_TOGGLE, 1, &theUseOpenBoundariesName, PRMzeroDefaults),

        PRM_Template(PRM_TOGGLE, 1, &theTrackBubbleIDsName, PRMzeroDefaults),

    	PRM_Template(PRM_TOGGLE, 1, &theTrackBubbleVolumesName, PRMzeroDefaults,
                        0, 0, 0, 0, 1, 0, &theTrackBubbleVolumesDisable),

    	PRM_Template(PRM_TOGGLE, 1, &theCorrectVolumeDriftName, PRMzeroDefaults,
                        0, 0, 0, 0, 1, 0, &theCorrectVolumeDriftDisable),

        PRM_Template(PRM_STRING, 1, &theGeometryName, &theGeometryNameDefault),

    	PRM_Template(PRM_STRING, 1, &theBubbleVolumeGeometryName, &theBubbleVolumeGeometryNameDefault),

    	PRM_Template()
    };

    static SIM_DopDescription theDopDescription(true,
						"HDK_ConstraintBubblePressureSolver",
						"HDK Constraint Bubble Pressure Solver",
						"$OS",
						classname(),
						theTemplates);

    setGasDescription(theDopDescription);

    return &theDopDescription;
}

bool HDK_ConstraintBubblePressureSolver::solveGasSubclass(SIM_Engine &engine,
							SIM_Object *obj,
							SIM_Time time,
							SIM_Time timestep)
{
    const SolveReal dt = timestep;

    const SIM_VectorField *solidVelocity = getConstVectorField(obj, GAS_NAME_COLLISIONVELOCITY);

    // Load liquid velocity

    SIM_VectorField *velocity = getVectorField(obj, GAS_NAME_VELOCITY);

    if (velocity == nullptr)
    {
        addError(obj, SIM_MESSAGE, "Velocity field missing", UT_ERROR_WARNING);
        return false;
    }
    else if (!velocity->isFaceSampled())
    {
        addError(obj, SIM_MESSAGE, "Velocity field must be a staggered grid", UT_ERROR_WARNING);
        return false;
    }

    // Load cut-cell weights

    const SIM_VectorField *cutCellWeightsField = getConstVectorField(obj, "cutCellWeights");
    
    if (cutCellWeightsField == nullptr)
    {
        addError(obj, SIM_MESSAGE, "Cut-cell weights field missing", UT_ERROR_WARNING);
        return false;
    }
    else if (!cutCellWeightsField->isAligned(velocity))
    {
        addError(obj, SIM_MESSAGE, "Cut-cell weights must align with velocity samples", UT_ERROR_WARNING);
        return false;
    }
    
    std::array<const SIM_RawField*, 3> cutCellWeights;
    for (int axis : {0, 1, 2})
	cutCellWeights[axis] = cutCellWeightsField->getField(axis);
    
    // Load valid fluid faces

    SIM_VectorField *validFaces = getVectorField(obj, "validFaces");

    if (validFaces == nullptr)
    {
        addError(obj, SIM_MESSAGE, "No 'valid' field found", UT_ERROR_ABORT);
        return false;
    }
    else if (!validFaces->isAligned(velocity))
    {
        addError(obj, SIM_MESSAGE, "Valid field sampling needs to match velocity field", UT_ERROR_ABORT);
        return false;
    }

    // Load liquid SDF

    const SIM_ScalarField *liquidSurfaceField = getConstScalarField(obj, GAS_NAME_SURFACE);

    if (liquidSurfaceField == nullptr)
    {
        addError(obj, SIM_MESSAGE, "Surface field is missing. There is nothing to represent the liquid", UT_ERROR_WARNING);
        return false;
    }

    const SIM_RawField &liquidSurface = *liquidSurfaceField->getField();

    // Load liquid pressure

    SIM_ScalarField *pressureField = getScalarField(obj, GAS_NAME_PRESSURE, true);
    SIM_RawField localPressure, *pressure;

    const UT_Vector3I voxelRes = velocity->getTotalVoxelRes();

    if (pressureField == nullptr)
    {
	localPressure.init(SIM_SAMPLE_CENTER,
			    velocity->getOrig(),
			    velocity->getSize(),
			    voxelRes[0], voxelRes[1], voxelRes[2]);
	localPressure.makeConstant(0);
	pressure = &localPressure;
    }
    else
    {
	pressureField->matchField(liquidSurfaceField);
	pressure = pressureField->getField();
    }

    // Load solid SDF

    const SolveReal dx = velocity->getVoxelSize().maxComponent();

    const SIM_ScalarField *solidSurfaceField = getConstScalarField(obj, GAS_NAME_COLLISION);
    const SIM_RawField *solidSurface;
    SIM_RawField localSolidSurface;

    if (solidSurfaceField == nullptr)
    {
        // Treat as all fluid.
        localSolidSurface.init(SIM_SAMPLE_CENTER,
				velocity->getOrig(),
				velocity->getSize(),
				voxelRes[0], voxelRes[1], voxelRes[2]);

        // Solid level set convention in Houdini is negative outside and positive inside.
        localSolidSurface.makeConstant(-10. * dx);
        solidSurface = &localSolidSurface;
    }
    else
	solidSurface = solidSurfaceField->getField();

    const SIM_ScalarField *surfacePressureField = getConstScalarField(obj, "surfacepressure");

    if (surfacePressureField != nullptr && !surfacePressureField->getField()->isAligned(&liquidSurface))
    {
	addError(obj, SIM_MESSAGE, "Surface pressure must align with the liquid surface field", UT_ERROR_WARNING);
	return false;
    }

    // Load liquid density

    const SIM_ScalarField *liquidDensityField = getConstScalarField(obj, GAS_NAME_DENSITY);

    if (liquidDensityField == nullptr)
    {
        addError(obj, SIM_MESSAGE, "There is no liquid density to simulate with", UT_ERROR_WARNING);
        return false;
    }

    if (!liquidDensityField->getField()->isAligned(&liquidSurface))
    {
        addError(obj, SIM_MESSAGE, "Density must align with the surface volume", UT_ERROR_WARNING);
        return false;
    }

    fpreal32 constantLiquidDensity = 0.;
    if (!liquidDensityField->getField()->field()->isConstant(&constantLiquidDensity))
    {
	addError(obj, SIM_MESSAGE, "Variable density is not currently supported", UT_ERROR_WARNING);
        return false;
    }

    const SolveReal liquidDensity = constantLiquidDensity;

    std::cout << "//\n//\n// Starting constraint bubble solver\n//\n//" << std::endl;
    
    ////////////////////////////////////////////
    //
    // Build material labels
    //
    ////////////////////////////////////////////

    SIM_RawIndexField materialCellLabels;

    {
	UT_PerfMonAutoSolveEvent event(this, "Build cell material labels");
	std::cout << "\n// Build cell material labels" << std::endl;

	HDK::Utilities::buildMaterialCellLabels(materialCellLabels,
						liquidSurface,
						*solidSurface,
						cutCellWeights);

    	// Set labels specific to this solver
	setMaterialCellLabels(materialCellLabels);
    }

    ////////////////////////////////////////////
    //
    // Build a mapping of solvable pressure cells
    // to rows in the matrix.
    //
    ////////////////////////////////////////////

    exint liquidCellCount = 0;
    SIM_RawIndexField liquidCellIndices;

    {
    	UT_PerfMonAutoSolveEvent event(this, "Build liquid cell labels");

    	std::cout << "\n// Build liquid cell labels" << std::endl;

	liquidCellCount = buildLiquidCellIndices(liquidCellIndices, materialCellLabels);

    	UT_WorkBuffer extrainfo;
	extrainfo.sprintf("liquid DOFs=%d", int(liquidCellCount));
	event.setExtraInfo(extrainfo.buffer());
    }

    if (liquidCellCount == 0)
    {
        addError(obj, SIM_MESSAGE, "No liquid cells found", UT_ERROR_WARNING);
        return false;
    }

    std::cout << "    Liquid cell count: " << liquidCellCount << std::endl;

    ////////////////////////////////////////////
    //
    // Build valid faces to indicate active velocity elements
    //
    ////////////////////////////////////////////

    {
    	std::cout << "\n// Build valid flags" << std::endl;
    	UT_PerfMonAutoSolveEvent event(this, "Build valid flags");

    	buildValidFaces(*validFaces,
			materialCellLabels,
			cutCellWeights);
    }

    ////////////////////////////////////////////
    //
    // Integrate surface pressure into the velocity field
    //
    ////////////////////////////////////////////    

    if (surfacePressureField != nullptr)
    {
	UT_PerfMonAutoSolveEvent event(this, "Integrate Surface Pressures in Velocity Field");

	for (int axis : {0,1,2})
	{
	    integrateSurfacePressure(*velocity->getField(axis),
					*validFaces->getField(axis),
					*surfacePressureField->getField(),
					liquidSurface,
					liquidCellIndices,
					constantLiquidDensity,
					dt, 
					axis);
	}
    }

    ////////////////////////////////////////////
    //
    // Flood fill the air regions and give each bubble
    // a degree-of-freedom to be treated as a monolithic
    // solvable cell.
    //
    ////////////////////////////////////////////

    exint bubbleRegionCount;
    SIM_RawIndexField bubbleRegionIndices;
    UT_Array<exint> bubbleRegionSizes;

    {
    	std::cout << "\n// Build bubble region labels" << std::endl;
    	UT_PerfMonAutoSolveEvent event(this, "Build bubble region labels");

    	bubbleRegionCount = buildBubbleRegionIndices(bubbleRegionIndices,
							bubbleRegionSizes,
							materialCellLabels,
							cutCellWeights);

    	UT_WorkBuffer extrainfo;
    	extrainfo.sprintf("bubble DOFs=%d", int(bubbleRegionCount));
    	event.setExtraInfo(extrainfo.buffer());

    	std::cout << "  Number of bubbles: " << bubbleRegionCount << std::endl;

    	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    std::cout << "    Bubble " << bubbleRegion << ". Size: " << bubbleRegionSizes[bubbleRegion] << std::endl;
    }

    ////////////////////////////////////////////
    //
    // Build connected components for non-solid regions.
    // This joins air and liquid regions together and is
    // needed to project out the null space when simulating
    // the entire volume. This is also useful for removing
    // air regions that have no air-liquid boundary from
    // the linear solver.
    //
    ////////////////////////////////////////////

    SIM_RawIndexField combinedRegionIndices;
    exint combinedRegionCount;

    {
    	std::cout << "\n// Building combined connected components" << std::endl;
    	UT_PerfMonAutoSolveEvent event(this, "Build combined connected components");

	combinedRegionCount = buildCombinedRegionIndices(combinedRegionIndices, materialCellLabels, cutCellWeights);

	std::cout << "  Combined region count: " << combinedRegionCount << std::endl;
    }

    ////////////////////////////////////////////
    //
    // Build a list of bubbles for each combined region and
    // check if the region contains any liquid cells.
    // By building a list of bubbles for each region,
    // we can easily check if the region contains
    // a removed bubble. We don't want to enforce a bubble
    // constraint if there is no liquid
    // so we can remove stray bubbles on the list.
    //
    ////////////////////////////////////////////

    UT_Array<UT_Array<bool>> isBubbleInCombinedRegion;
    UT_Array<bool> isLiquidInCombinedRegion;

    {
    	std::cout << "\n// Build list of bubbles per combined region." << std::endl;

    	UT_PerfMonAutoSolveEvent event(this, "Build list of bubbles per combined region");

    	buildRegionBooleans(isBubbleInCombinedRegion,
			    isLiquidInCombinedRegion,
			    materialCellLabels,
			    liquidCellIndices,
			    bubbleRegionIndices,
			    combinedRegionIndices,
			    combinedRegionCount,
			    bubbleRegionCount);
    }

    ////////////////////////////////////////////
    //
    // If bubble tracking is active, let's build
    // a table of old bubbles linking to new bubbles
    // and update the bubble ID attribute on the particles.
    //
    ////////////////////////////////////////////

    bool doTrackBubbleVolumes = getTrackBubbleVolumes();

    // If we're tracking volumes we have to track IDs.
    bool doTrackBubbleIDs = doTrackBubbleVolumes || getTrackBubbleIDs();

    UT_Array<bool> isBubbleTracked;
    isBubbleTracked.setSize(bubbleRegionCount);
    isBubbleTracked.constant(false);

    UT_Array<SolveReal> newBubbleRegionVolumes;
    newBubbleRegionVolumes.setSize(bubbleRegionCount);
    newBubbleRegionVolumes.constant(0);

    const int threadCount = UT_Thread::getNumProcessors();

    if (doTrackBubbleIDs)
    {
    	std::cout << "\n// Track bubble IDs" << std::endl;

    	// Get particle geometry to pull bubble ID data from.
    	SIM_GeometryCopy *particleGeometry = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    	SIM_GeometryAutoWriteLock autoLockParticles(particleGeometry, SIM_DATA_ID_PRESERVE);
    	GU_Detail *particles = &autoLockParticles.getGdp();

    	assert(particles != nullptr);

    	GA_RWHandleI bubbleIDHandle(particles, GA_ATTRIB_POINT, "bubbleID");

    	if (!bubbleIDHandle.isValid())
    	{
    	    // Set default value for particles to UNINITIALIZED_PARTICLE. Any bubbles with UNINITIALIZED_PARTICLE
    	    // particle labels will not be removed from the system.
    	    particles->addIntTuple(GA_ATTRIB_POINT, "bubbleID", 1, GA_Defaults(UNINITIALIZED_PARTICLE));
    	    bubbleIDHandle = GA_RWHandleI(particles, GA_ATTRIB_POINT, "bubbleID");
    	    bubbleIDHandle.bumpDataId();
    	}

    	// Grab geometry for the previous timestep's bubble volume.
    	exint oldBubbleRegionCount = 0;
    	SIM_GeometryCopy *bubbleVolumeGeometry = nullptr;
    	GU_Detail *bubbleVolumeDetail = nullptr;

    	UT_UniquePtr<SIM_GeometryAutoWriteLock> autoLockBubbleVolume;

    	if (doTrackBubbleVolumes)
    	{
    	    bubbleVolumeGeometry = getOrCreateGeometry(obj, "bubbleVolumeGeometry");
    	    autoLockBubbleVolume = UTmakeUnique<SIM_GeometryAutoWriteLock>(bubbleVolumeGeometry, SIM_DATA_ID_PRESERVE);
    	    bubbleVolumeDetail = &autoLockBubbleVolume->getGdp();

    	    assert(bubbleVolumeDetail != nullptr);

    	    oldBubbleRegionCount = bubbleVolumeDetail->getNumPoints();
    	    std::cout << "  Debug check. Tracked bubble count: " << oldBubbleRegionCount << std::endl;
    	}

    	UT_Array<UT_Array<bool>> oldToNewBubbleMap;
	oldToNewBubbleMap.setSize(oldBubbleRegionCount);

    	for (exint oldBubbleRegion = 0; oldBubbleRegion < oldBubbleRegionCount; ++oldBubbleRegion)
	{
    	    oldToNewBubbleMap[oldBubbleRegion].setSize(bubbleRegionCount);
    	    oldToNewBubbleMap[oldBubbleRegion].constant(false);
	}

    	// Track bubble IDs by filling in array. If a particle with a valid bubble ID is adjacent to a
    	// newly detected bubble, indicate that the new bubble is indeed tracked.
    	{
    	    UT_PerfMonAutoSolveEvent event(this, "Track bubble IDs");
    	    std::cout << "  Track bubble IDs" << std::endl;

    	    trackBubbleIDs(isBubbleTracked,
			    oldToNewBubbleMap,
			    bubbleIDHandle,
			    *particles,
			    bubbleRegionIndices,
			    liquidSurface);

    	    // Update bubble ID on particles for particles adjacent to bubbles that were tracked
    	    // from the previous timestep.
	    updateBubbleIDs(bubbleIDHandle, isBubbleTracked, *particles, bubbleRegionIndices, liquidSurface);
	    bubbleIDHandle.bumpDataId();
    	}

    	if (doTrackBubbleVolumes)
    	{
    	    std::cout << "\n// Build bubble volume graph" << std::endl;
    	    UT_PerfMonAutoSolveEvent event(this, "Build bubble volume graph");

    	    assert(bubbleVolumeDetail != nullptr);
    	    {
		    GA_RWHandleF bubbleRestVolumeHandle(bubbleVolumeDetail, GA_ATTRIB_POINT, "bubbleRestVolume");

		    if (!bubbleRestVolumeHandle.isValid())
		    {
		    bubbleVolumeDetail->addFloatTuple(GA_ATTRIB_POINT, "bubbleRestVolume", 1, GA_Defaults(-1));
		    bubbleRestVolumeHandle = GA_RWHandleF(bubbleVolumeDetail, GA_ATTRIB_POINT, "bubbleRestVolume");
		    bubbleRestVolumeHandle.bumpDataId();
		    }

		    UT_Array<bool> isBubbleUninitialized;
		    isBubbleUninitialized.setSize(bubbleRegionCount);
		    isBubbleUninitialized.constant(false);

		    distributeOldBubbleVolumes(newBubbleRegionVolumes,
						isBubbleUninitialized,
						bubbleRestVolumeHandle,
						oldToNewBubbleMap,
						isBubbleTracked,
						bubbleRegionSizes);

		    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
		    {
		    if (isBubbleUninitialized[bubbleRegion])
			newBubbleRegionVolumes[bubbleRegion] = SolveReal(bubbleRegionSizes[bubbleRegion]);
		    }

		    // Correct for any global volume drift
		GA_RWHandleF globalBubbleRestVolumeHandle(bubbleVolumeDetail, GA_ATTRIB_POINT, "globalBubbleRestVolume");

		if (!globalBubbleRestVolumeHandle.isValid())
		    {
		    bubbleVolumeDetail->addFloatTuple(GA_ATTRIB_POINT, "globalBubbleRestVolume", 1, GA_Defaults(-1));
		    globalBubbleRestVolumeHandle = GA_RWHandleF(bubbleVolumeDetail, GA_ATTRIB_POINT, "globalBubbleRestVolume");
		    globalBubbleRestVolumeHandle.bumpDataId();
		    }

		SolveReal globalBubbleVolume;
		if (oldBubbleRegionCount > 0)
		    globalBubbleVolume = globalBubbleRestVolumeHandle.get(0);
		else
		{
		    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
			globalBubbleVolume += SolveReal(newBubbleRegionVolumes[bubbleRegion]);
		}

		    if (globalBubbleVolume > 0)
		{
		    SolveReal currentRestVolume = 0;
		    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
			currentRestVolume += SolveReal(newBubbleRegionVolumes[bubbleRegion]);

		    SolveReal volumeDiff = globalBubbleVolume - currentRestVolume;

		    std::cout << "  Bubble region correction for global rest volume: " << volumeDiff << std::endl;

		    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
			newBubbleRegionVolumes[bubbleRegion] += volumeDiff * newBubbleRegionVolumes[bubbleRegion] / currentRestVolume;

		    // Debug check
		    currentRestVolume = 0;
		    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
			currentRestVolume += SolveReal(newBubbleRegionVolumes[bubbleRegion]);

		    volumeDiff = globalBubbleVolume - currentRestVolume;

		    std::cout << "  Post bubble region correction for global rest volume: " << volumeDiff << std::endl;
		}
	    }

	    {
		// If volume tracking is enabled, the rest volume list has to be remapped.
		// After the remapping, we need to update particle IDs and the bubble volume detail.
		UT_PerfMonAutoSolveEvent event(this, "Update bubble volumes");

		bubbleVolumeDetail->clear();

		GA_RWHandleF bubbleRestVolumeHandle(bubbleVolumeDetail, GA_ATTRIB_POINT, "bubbleRestVolume");
		if (!bubbleRestVolumeHandle.isValid())
		{
		    bubbleVolumeDetail->addFloatTuple(GA_ATTRIB_POINT, "bubbleRestVolume", 1, GA_Defaults(-1));
		    bubbleRestVolumeHandle = GA_RWHandleF(bubbleVolumeDetail, GA_ATTRIB_POINT, "bubbleRestVolume");
		    bubbleRestVolumeHandle.bumpDataId();
		}

		// Correct for any global volume drift
		GA_RWHandleF globalBubbleRestVolumeHandle(bubbleVolumeDetail, GA_ATTRIB_POINT, "globalBubbleRestVolume");

		if (!globalBubbleRestVolumeHandle.isValid())
		{
		    bubbleVolumeDetail->addFloatTuple(GA_ATTRIB_POINT, "globalBubbleRestVolume", 1, GA_Defaults(-1));
		    globalBubbleRestVolumeHandle = GA_RWHandleF(bubbleVolumeDetail, GA_ATTRIB_POINT, "globalBubbleRestVolume");
		    globalBubbleRestVolumeHandle.bumpDataId();
		}

		SolveReal currentRestVolume = 0;
		for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
		    currentRestVolume += SolveReal(newBubbleRegionVolumes[bubbleRegion]);

		GA_Offset pointOffset = bubbleVolumeDetail->appendPointBlock(bubbleRegionCount);

		for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion, ++pointOffset)
		{
		    if (isBubbleTracked[bubbleRegion])
			bubbleRestVolumeHandle.set(pointOffset, newBubbleRegionVolumes[bubbleRegion]);

		    globalBubbleRestVolumeHandle.set(pointOffset, currentRestVolume);
		}
	    }
    	}
    }

    ////////////////////////////////////////////
    //
    // Delete bubbles that have no liquid boundary.
    // Presumably this means the air is ambient and
    // outside of a closed container.
    //
    ////////////////////////////////////////////

    UT_Array<bool> isBubbleRemoved;
    isBubbleRemoved.setSize(bubbleRegionCount);
    isBubbleRemoved.constant(false);

    {
        std::cout << "\n// Removing bubbles with no liquid boundary" << std::endl;
        for (exint combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
        {
            if (!isLiquidInCombinedRegion[combinedRegion])
            {
                int bubblesInRegionCount = 0;
                for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
                {
                    if (isBubbleInCombinedRegion[combinedRegion][bubbleRegion])
                    {
			++bubblesInRegionCount;

			if (doTrackBubbleIDs)
			    assert(isBubbleRemoved[bubbleRegion]);
			else
			{
			    isBubbleRemoved[bubbleRegion] = true;
			    std::cout << "  Bubble " << bubbleRegion << " in combined region " << combinedRegion << " has been flagged for removal" << std::endl;
			}
                    }
                }

                // Debug check
                // There cannot be two bubbles in a single
                // region without a liquid region separating them.
                assert(bubblesInRegionCount == 1);
            }
        }
    }

    if (getUseOpenBoundaries())
    {
        std::cout << "\n// Remove bubbles at open boundaries" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Remove bubbles at open boundaries");
	
	UT_Array<bool> isBubbleAtOpenBoundary;
	isBubbleAtOpenBoundary.setSize(bubbleRegionCount);
	isBubbleAtOpenBoundary.constant(0);

	removeBubblesAtOpenBoundaries(isBubbleAtOpenBoundary,
					bubbleRegionIndices,
					cutCellWeights);

	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (isBubbleAtOpenBoundary[bubbleRegion])
	    {
		std::cout << "  Bubble region " << bubbleRegion << " at open boundary removed" << std::endl;
		isBubbleRemoved[bubbleRegion] = true;
	    }
	}
    }

    ////////////////////////////////////////////
    //
    // Build divergence for each liquid cell
    //
    ////////////////////////////////////////////

    const exint totalDOFCount = liquidCellCount + bubbleRegionCount;

    Vector rhsVector = Vector::Zero(totalDOFCount);

    {
    	std::cout << "\n// Build liquid RHS" << std::endl;
    	UT_PerfMonAutoSolveEvent event(this, "Build liquid right-hand side");

    	buildLiquidRHS(rhsVector,
			materialCellLabels,
			liquidCellIndices,
			bubbleRegionIndices,
			combinedRegionIndices,
			*velocity,
    			solidVelocity,
    			cutCellWeights);
    }

    ////////////////////////////////////////////
    //
    // Build divergence for each bubble region
    //
    ////////////////////////////////////////////

    {
    	std::cout << "\n// Build bubble RHS" << std::endl;
    	UT_PerfMonAutoSolveEvent event(this, "Build bubble RHS");

	UT_Array<UT_Array<SolveReal>> parallelBubbleRegionDivergence;
	parallelBubbleRegionDivergence.setSize(threadCount);

   	for (int thread = 0; thread < threadCount; ++thread)
	{
	    parallelBubbleRegionDivergence[thread].setSize(bubbleRegionCount);
	    parallelBubbleRegionDivergence[thread].constant(0);
	}

    	buildBubbleRHS(parallelBubbleRegionDivergence,
			materialCellLabels,
			liquidCellIndices,
			bubbleRegionIndices,
			combinedRegionIndices,
			cutCellWeights,
			isBubbleRemoved,
			*velocity,
    			solidVelocity,
			bubbleRegionCount);

    	for (int thread = 0; thread < threadCount; ++thread)
    	    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
		rhsVector(bubbleRegion + liquidCellCount) += parallelBubbleRegionDivergence[thread][bubbleRegion];
    }

    ////////////////////////////////////////////
    //
    // Remove the null space introduced by latent
    // divergence in the system
    //
    ////////////////////////////////////////////

    {
	UT_Array<bool> doesRegionHaveRemovedBubble;
	doesRegionHaveRemovedBubble.setSize(combinedRegionCount);
	doesRegionHaveRemovedBubble.constant(0);

	for (exint combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
	    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    {
		if (isBubbleInCombinedRegion[combinedRegion][bubbleRegion])
		{
		    if (isBubbleRemoved[bubbleRegion])
			doesRegionHaveRemovedBubble[combinedRegion] = true;
		}
	    }

	UT_Array<exint> combinedRegionLiquidCellCount;
	combinedRegionLiquidCellCount.setSize(combinedRegionCount);
	combinedRegionLiquidCellCount.constant(0);

	UT_Array<exint> combinedRegionBubbleCellCount;
	combinedRegionBubbleCellCount.setSize(combinedRegionCount);
	combinedRegionBubbleCellCount.constant(0);

	{
	    std::cout << "\n// Project null space" << std::endl;
	    UT_PerfMonAutoSolveEvent event(this, "Project null space");

	    // Regions with no removed bubbles must have a zero average RHS.
	    UT_Array<SolveReal> combinedRegionDivergence;
	    combinedRegionDivergence.setSize(combinedRegionCount);
	    combinedRegionDivergence.constant(0);

	    for (exint combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
	    {
		if (!doesRegionHaveRemovedBubble[combinedRegion])
		{
		    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
		    {
			if (isBubbleInCombinedRegion[combinedRegion][bubbleRegion])
			{
			    combinedRegionDivergence[combinedRegion] += rhsVector(bubbleRegion + liquidCellCount);
			    combinedRegionBubbleCellCount[combinedRegion] += bubbleRegionSizes[bubbleRegion];
			}
		    }
		}
	    }

	    UT_Array<UT_Array<SolveReal>> parallelLiquidRegionDivergence;
	    parallelLiquidRegionDivergence.setSize(threadCount);

	    UT_Array<UT_Array<exint>> parallelLiquidRegionCellCount;
	    parallelLiquidRegionCellCount.setSize(threadCount);

	    for (int thread = 0; thread < threadCount; ++thread)
	    {
		parallelLiquidRegionDivergence[thread].setSize(combinedRegionCount);
		parallelLiquidRegionDivergence[thread].constant(0);

		parallelLiquidRegionCellCount[thread].setSize(combinedRegionCount);
		parallelLiquidRegionCellCount[thread].constant(0);
	    }

	    buildLiquidRegionDivergence(parallelLiquidRegionDivergence,
					parallelLiquidRegionCellCount,
					combinedRegionIndices,
					liquidCellIndices,
					doesRegionHaveRemovedBubble,
					rhsVector);

	    for (int thread = 0; thread < threadCount; ++thread)
		for (exint combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
		{
		    if (!doesRegionHaveRemovedBubble[combinedRegion])
		    {
			combinedRegionDivergence[combinedRegion] += parallelLiquidRegionDivergence[thread][combinedRegion];
			combinedRegionLiquidCellCount[combinedRegion] += parallelLiquidRegionCellCount[thread][combinedRegion];
		    }
		    else
		    {
			assert(parallelLiquidRegionDivergence[thread][combinedRegion] == 0);
			assert(parallelLiquidRegionCellCount[thread][combinedRegion] == 0);
		    }
		}

	    // Build averages per-cell divergence
	    for (exint combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
	    {
		std::cout << "  Region: " << combinedRegion << ". Divergence: " << combinedRegionDivergence[combinedRegion] << std::endl;
		
		if (!doesRegionHaveRemovedBubble[combinedRegion])
		{
    		    assert((combinedRegionLiquidCellCount[combinedRegion] + combinedRegionBubbleCellCount[combinedRegion]) > 0);
		    combinedRegionDivergence[combinedRegion] /= SolveReal(combinedRegionLiquidCellCount[combinedRegion] + combinedRegionBubbleCellCount[combinedRegion]);
		}
		else assert(combinedRegionDivergence[combinedRegion] == 0);
	    }

	    // Subtract average divergence from bubbles
	    for (exint combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
	    {
		if (!doesRegionHaveRemovedBubble[combinedRegion])
		{
		    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
		    {
			if (isBubbleInCombinedRegion[combinedRegion][bubbleRegion])
			    rhsVector(bubbleRegion + liquidCellCount) -= combinedRegionDivergence[combinedRegion] * SolveReal(bubbleRegionSizes[bubbleRegion]);
		    }
		}
	    }

	    // Subtract divergence from liquid cells
	    removeDivergenceFromLiquidCells(rhsVector,
					    combinedRegionDivergence,
					    combinedRegionIndices,
					    liquidCellIndices,
					    doesRegionHaveRemovedBubble);

	    std::cout << "  Average divergence in rhs vector after latent divergence removal: " << rhsVector.sum() / SolveReal(rhsVector.rows()) << std::endl;
	}

	////////////////////////////////////////////
	//
	// Now apply the divergence per bubble region
	// accounting for untracked bubbles and volume correction.
	//
	////////////////////////////////////////////

	if (doTrackBubbleIDs)
	{
	    std::cout << "\n// Add bubble tracking to RHS" << std::endl;
	    UT_PerfMonAutoSolveEvent event(this, "Add bubble tracking to RHS");

	    // Untracked bubbles must be killed off so let's add divergence
	    // to drive the volume to zero. We do this using dV/dt = div(u).

	    UT_Array<SolveReal> addedBubbleDivergence;
	    addedBubbleDivergence.setSize(bubbleRegionCount);
	    addedBubbleDivergence.constant(0);

	    bool doRemoveDivergence = false;
	    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    {
		if (!isBubbleTracked[bubbleRegion])
		{
		    // We should be scaling the divergence by cell volume but we have implicitly
		    // divided out a dx^2 factor in the linear system so we should only scale by dx.
		    SolveReal localDivergence = -SolveReal(bubbleRegionSizes[bubbleRegion]) * dx / dt;
		    std::cout << "  Untracked bubble: " << bubbleRegion << " of size " << bubbleRegionSizes[bubbleRegion]
				<< " is set to collapse with divergence " << localDivergence << std::endl;
		
		    addedBubbleDivergence[bubbleRegion] += localDivergence;
		    doRemoveDivergence = true;
		}
	    }

	    doRemoveDivergence |= doTrackBubbleVolumes;

	    if (doTrackBubbleVolumes)
	    {
		SolveReal correctionScale = .1;

		for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
		{
		    if (isBubbleTracked[bubbleRegion] && !isBubbleRemoved[bubbleRegion])
		    {
			std::cout << " Bubble: " << bubbleRegion << std::endl;
			std::cout << " Current bubble size: " << bubbleRegionSizes[bubbleRegion] << std::endl;
			std::cout << " Rest volume: " << newBubbleRegionVolumes[bubbleRegion] << std::endl;

			if (getCorrectVolumeDrift())
			{
			    SolveReal localDivergence = (newBubbleRegionVolumes[bubbleRegion] - SolveReal(bubbleRegionSizes[bubbleRegion])) * dx / dt;
			    std::cout << " Set divergence: " << localDivergence << std::endl;

			    addedBubbleDivergence[bubbleRegion] += correctionScale * localDivergence;
			}
		    }
		}
	    }

	    if (doRemoveDivergence)
	    {
		UT_Array<SolveReal> combinedRegionDivergence;
		combinedRegionDivergence.setSize(combinedRegionCount);
		combinedRegionDivergence.constant(0);

		// Build region divergence list
		for (exint combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
		{
		    if (!doesRegionHaveRemovedBubble[combinedRegion])
		    {
			for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
			{
			    if (isBubbleInCombinedRegion[combinedRegion][bubbleRegion])
			    {
				rhsVector(bubbleRegion + liquidCellCount) += addedBubbleDivergence[bubbleRegion];
				combinedRegionDivergence[combinedRegion] += addedBubbleDivergence[bubbleRegion];
			    }
			}

			combinedRegionDivergence[combinedRegion] /= SolveReal(combinedRegionLiquidCellCount[combinedRegion]);
		    }
		}

		// Subtract divergence from liquid cells
		removeDivergenceFromLiquidCells(rhsVector,
						combinedRegionDivergence,
						combinedRegionIndices,
						liquidCellIndices,
						doesRegionHaveRemovedBubble);

		std::cout << "  Average divergence in rhs vector after tracking divergence removal: " << rhsVector.sum() / SolveReal(rhsVector.rows()) << std::endl;
	    }
	}
    }

    ////////////////////////////////////////////
    //
    // Build matrix rows for liquid degrees of freedom.
    // Add non-zeros for entries for adjacent bubble regions.
    //
    ////////////////////////////////////////////

    std::vector<Eigen::Triplet<SolveReal>> sparseElements;
    sparseElements.reserve(totalDOFCount * 7);

    {
        std::cout << "\n// Build liquid matrix rows" << std::endl;
        UT_PerfMonAutoSolveEvent event(this, "Build liquid matrix rows");

	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelSparseElements(threadCount);

    	buildLiquidRows(parallelSparseElements,
			liquidCellIndices,
			bubbleRegionIndices,
    			liquidSurface,
    			cutCellWeights,
    			isBubbleRemoved,
			liquidDensity,
			liquidCellCount);

	exint elementListSize = 0;
	for (int thread = 0; thread < threadCount; ++thread)
	    elementListSize += parallelSparseElements[thread].size();

	sparseElements.reserve(elementListSize);
    	for (int thread = 0; thread < threadCount; ++thread)
    	    sparseElements.insert(sparseElements.end(), parallelSparseElements[thread].begin(), parallelSparseElements[thread].end());
    }

    ////////////////////////////////////////////
    //
    // Build bubble rows
    //
    ////////////////////////////////////////////

    UT_Array<SolveReal> bubbleRegionDiagonals;
    bubbleRegionDiagonals.setSize(bubbleRegionCount);
    bubbleRegionDiagonals.constant(0);

    {
    	std::cout << "\n// Build bubble matrix rows" << std::endl;
    	UT_PerfMonAutoSolveEvent event(this, "Build bubble matrix rows");

        std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelSparseElements(threadCount);

	UT_Array<UT_Array<SolveReal>> parallelBubbleRegionDiagonals;
	parallelBubbleRegionDiagonals.setSize(threadCount);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    parallelBubbleRegionDiagonals[thread].setSize(bubbleRegionCount);
	    parallelBubbleRegionDiagonals[thread].constant(0);
	}

    	buildBubbleRows(parallelSparseElements,
    			parallelBubbleRegionDiagonals,
			liquidCellIndices,
			bubbleRegionIndices,
    			liquidSurface,
    			cutCellWeights,
    			isBubbleRemoved,
			liquidDensity,
			liquidCellCount);

	exint elementListSize = 0;
	for (int thread = 0; thread < threadCount; ++thread)
	    elementListSize += parallelSparseElements[thread].size();

	sparseElements.reserve(elementListSize);
    	for (int thread = 0; thread < threadCount; ++thread)
	    sparseElements.insert(sparseElements.end(), parallelSparseElements[thread].begin(), parallelSparseElements[thread].end());

	// Collect diagonals per thread
	for (int thread = 0; thread < threadCount; ++thread)
	    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    {
		if (!isBubbleRemoved[bubbleRegion])
		{
		    if (parallelBubbleRegionDiagonals[thread][bubbleRegion] != 0)
		    {
			assert(parallelBubbleRegionDiagonals[thread][bubbleRegion] > 0);
			bubbleRegionDiagonals[bubbleRegion] += parallelBubbleRegionDiagonals[thread][bubbleRegion];
		    }
		}
		else assert(parallelBubbleRegionDiagonals[thread][bubbleRegion] == 0);
	    }

	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (!isBubbleRemoved[bubbleRegion])
		sparseElements.push_back(Eigen::Triplet<SolveReal>(bubbleRegion + liquidCellCount, bubbleRegion + liquidCellCount, bubbleRegionDiagonals[bubbleRegion]));
	}

	// Set 1 along diagonal for removed bubble
	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (isBubbleRemoved[bubbleRegion])
	    {
		assert(rhsVector(bubbleRegion + liquidCellCount) == 0);
		sparseElements.push_back(Eigen::Triplet<SolveReal>(bubbleRegion + liquidCellCount, bubbleRegion + liquidCellCount, 1));

		bubbleRegionDiagonals[bubbleRegion] = 1;
		assert(rhsVector(bubbleRegion + liquidCellCount) == 0);
	    }
	}
    }

    Vector solutionVector = Vector::Zero(totalDOFCount);

    ////////////////////////////////////////////
    //
    // Set initial guess for liquid pressure
    //
    ////////////////////////////////////////////

    if (getUseOldPressure())
    {
	std::cout << "\n// Apply warm start pressure" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Apply warm start pressure");

	HDK::Utilities::applyInitialGuess(solutionVector,
					    *pressure,
					    liquidCellIndices);

	UT_Array<UT_Array<SolveReal>> parallelAccumulatedPressure;
	parallelAccumulatedPressure.setSize(threadCount);

   	for (int thread = 0; thread < threadCount; ++thread)
	{
	    parallelAccumulatedPressure[thread].setSize(bubbleRegionCount);
	    parallelAccumulatedPressure[thread].constant(0);
	}

	computeBubblePressure(parallelAccumulatedPressure,
				*pressure,
				bubbleRegionIndices,
				isBubbleRemoved,
				bubbleRegionCount);

	UT_Array<SolveReal> accumulatedPressure;
	accumulatedPressure.setSize(bubbleRegionCount);
	accumulatedPressure.constant(0);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    {
		if (!isBubbleRemoved[bubbleRegion])
		    accumulatedPressure[bubbleRegion] += parallelAccumulatedPressure[thread][bubbleRegion];
	    }
	}

	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (!isBubbleRemoved[bubbleRegion])
	    {
		assert(bubbleRegionSizes[bubbleRegion] > 0);
		solutionVector(bubbleRegion + liquidCellCount) = accumulatedPressure[bubbleRegion] / SolveReal(bubbleRegionSizes[bubbleRegion]);
	    }
	}
    }

    ////////////////////////////////////////////
    //
    // Solve for pressure.
    //
    ////////////////////////////////////////////

    {
	std::cout << "\n// Solve liquid system" << std::endl;

 	UT_PerfMonAutoSolveEvent event(this, "Solve linear system");

	Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> sparseMatrix(totalDOFCount, totalDOFCount);
	sparseMatrix.setFromTriplets(sparseElements.begin(), sparseElements.end());
	sparseMatrix.makeCompressed();
	
#if !defined(NDEBUG)

	std::cout << "  Checking symmetry of system" << std::endl;

	for (exint k = 0; k < sparseMatrix.outerSize(); ++k)
	    for (typename Eigen::SparseMatrix<SolveReal, Eigen::RowMajor>::InnerIterator it(sparseMatrix, k); it; ++it)
	    {
		if (!SYSisEqual(sparseMatrix.coeff(it.row(), it.col()), sparseMatrix.coeff(it.col(), it.row())))
		{
		    std::cout << "Value at row " << it.row() << ", col " << it.col() << " is " << sparseMatrix.coeff(it.row(), it.col()) << std::endl;
		    std::cout << "Value at row " << it.col() << ", col " << it.row() << " is " << sparseMatrix.coeff(it.col(), it.row()) << std::endl;
		    assert(false);
		}
	    }
    
	// Check that it's finite
	std::cout << "  Check for NaNs\n" << std::endl;

	for (exint k = 0; k < sparseMatrix.outerSize(); ++k)
	    for (typename Eigen::SparseMatrix<SolveReal, Eigen::RowMajor>::InnerIterator it(sparseMatrix, k); it; ++it)
	    {
		if (!std::isfinite(it.value()))
		    assert(false);
	    }
#endif

	if (getUseMGPreconditioner())
	{
	    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

	    // Build MG Preconditioner for liquid-only regions.
	    UT_VoxelArray<int> mgDomainCellLabels;
	    std::array<UT_VoxelArray<SolveReal>, 3> mgBoundaryWeights;
	    
	    UT_Vector3I mgExpandedOffset;
	    int mgLevels;

	    {
		std::cout << "\n// Build Multigrid Preconditioner" << std::endl;
		UT_PerfMonAutoSolveEvent event(this, "Build Multigrid Preconditioner");

		UT_VoxelArray<int> baseDomainCellLabels;
		baseDomainCellLabels.size(liquidCellIndices.getVoxelRes()[0],
					    liquidCellIndices.getVoxelRes()[1],
					    liquidCellIndices.getVoxelRes()[2]);

		baseDomainCellLabels.constant(MGCellLabels::EXTERIOR_CELL);

		buildMGDomainLabels(baseDomainCellLabels,
				    materialCellLabels,
				    liquidCellIndices);

		// Build boundary weights
		std::array<UT_VoxelArray<SolveReal>, 3> baseBoundaryWeights;

		for (int axis : {0,1,2})
		{
		    UT_Vector3I size = baseDomainCellLabels.getVoxelRes();
		    ++size[axis];

		    baseBoundaryWeights[axis].size(size[0], size[1], size[2]);
		    baseBoundaryWeights[axis].constant(0);

		    buildMGBoundaryWeights(baseBoundaryWeights[axis],
					    *validFaces->getField(axis),
					    materialCellLabels,
					    bubbleRegionIndices,
					    baseDomainCellLabels,
					    liquidSurface,
					    *cutCellWeights[axis],
					    liquidDensity,
					    axis);
		}
		
		// Build expanded domain
		auto isExteriorCell = [](const int value) { return value == MGCellLabels::EXTERIOR_CELL; };
		auto isInteriorCell = [](const int value) { return value == MGCellLabels::INTERIOR_CELL; };
		auto isDirichletCell = [](const int value) { return value == MGCellLabels::DIRICHLET_CELL; };

		std::pair<UT_Vector3I, int> mgSettings = HDK::GeometricMultigridOperators::buildExpandedCellLabels(mgDomainCellLabels, baseDomainCellLabels, isExteriorCell, isInteriorCell, isDirichletCell);

		mgExpandedOffset = mgSettings.first;
		mgLevels = mgSettings.second;
		
		// Build expanded boundary weights
		for (int axis : {0,1,2})
		{
		    UT_Vector3I size = mgDomainCellLabels.getVoxelRes();
		    ++size[axis];

		    mgBoundaryWeights[axis].size(size[0], size[1], size[2]);
		    mgBoundaryWeights[axis].constant(0);

		    HDK::GeometricMultigridOperators::buildExpandedBoundaryWeights(mgBoundaryWeights[axis], baseBoundaryWeights[axis], mgDomainCellLabels, mgExpandedOffset, axis);
		}

		// Build boundary cells
		HDK::GeometricMultigridOperators::setBoundaryCellLabels(mgDomainCellLabels, mgBoundaryWeights);

	        assert(HDK::GeometricMultigridOperators::unitTestBoundaryCells<SolveReal>(mgDomainCellLabels, &mgBoundaryWeights));
		assert(HDK::GeometricMultigridOperators::unitTestExteriorCells(mgDomainCellLabels));
	    }

	    assert(unitTestMGLabels(mgDomainCellLabels, materialCellLabels, liquidCellIndices, bubbleRegionIndices, mgExpandedOffset));

	    HDK::GeometricMultigridPoissonSolver mgPreconditioner(mgDomainCellLabels,
								    mgBoundaryWeights,
								    mgLevels,
								    true /* use Gauss-Seidel smoother */);

	    mgLevels = mgPreconditioner.getMGLevels();

	    //
	    // Initialize bubble smoother
	    //

	    // Build band of liquid cells at the bubble boundary

	    auto cellCompare = [&](const UT_Vector3I &cellA, const UT_Vector3I &cellB)
	    {
		// Compare tile number first
		int tileNumberA = materialCellLabels.field()->indexToLinearTile(cellA[0], cellA[1], cellA[2]);
		int tileNumberB = materialCellLabels.field()->indexToLinearTile(cellB[0], cellB[1], cellB[2]);

		if (tileNumberA < tileNumberB)
		    return true;
		else if (tileNumberA == tileNumberB)
		{
		    if (cellA[2] < cellB[2])
			return true;
		    else if (cellA[2] == cellB[2])
		    {
			if (cellA[1] < cellB[1])
			    return true;
			else if (cellA[1] == cellB[1] &&
				    cellA[0] < cellB[0])
			return true;
		    }
		}

		return false;
	    };

	    UT_Array<UT_Vector3I> liquidSmootherCells;
	    UT_Array<UT_Vector3I> liquidCopyCells;
	    UT_Array<UT_Vector3I> bubbleBoundaryCells;
	    
	    {
		SIM_RawIndexField visitedLiquidCells;
		visitedLiquidCells.match(liquidSurface);
		visitedLiquidCells.makeConstant(UNVISITED_CELL);

		{
		    std::cout << "\n// Build layer of liquid cells in smoother" << std::endl;

		    UT_Array<UT_Array<UT_Vector3I>> parallelLiquidBandCells;
		    parallelLiquidBandCells.setSize(threadCount);

		    // Build initial list of exterior cells
		    buildInitialLiquidSmootherLayer(parallelLiquidBandCells, materialCellLabels, cutCellWeights);

		    UT_Array<UT_Vector3I> newLiquidCells;

		    UT_Array<bool> isTileOccupiedList;
		    isTileOccupiedList.setSize(visitedLiquidCells.field()->numTiles());

		    exint liquidSmootherBand = 6;
		    for (int layer = 0; layer < liquidSmootherBand; ++layer)
		    {
			// Combine parallel liquid band cells
			exint listSize = 0;
			for (int thread = 0; thread < threadCount; ++thread)
			    listSize += parallelLiquidBandCells[thread].size();

			newLiquidCells.clear();
			newLiquidCells.bumpCapacity(listSize);

			for (int thread = 0; thread < threadCount; ++thread)
			{
			    newLiquidCells.concat(parallelLiquidBandCells[thread]);
			    parallelLiquidBandCells[thread].clear();
			}

			UTparallelSort(newLiquidCells.begin(), newLiquidCells.end(), cellCompare);

			std::cout << "  Liquid smoother layer size: " << newLiquidCells.size() << std::endl;

			// Build tile list
			isTileOccupiedList.constant(false);
			HDK::Utilities::findOccupiedIndexTiles(isTileOccupiedList,
								newLiquidCells,
								visitedLiquidCells);

			HDK::Utilities::uncompressTiles(visitedLiquidCells, isTileOccupiedList);

			setVisitedLiquidCells(visitedLiquidCells, newLiquidCells);

			if (layer < liquidSmootherBand - 1)
			{
			    // Build next exterior band list
			    buildNextLiquidSmootherLayer(parallelLiquidBandCells,
							    newLiquidCells,
							    materialCellLabels,
							    visitedLiquidCells,
							    cutCellWeights);
			}
		    }
		}

		{
		    std::cout << "\n// Build list of bubble smoother cells and liquid copy cells" << std::endl;

		    UT_Array<UT_Array<UT_Vector3I>> parallelLiquidSmootherCells;
		    parallelLiquidSmootherCells.setSize(threadCount);

		    UT_Array<UT_Array<UT_Vector3I>> parallelLiquidCopyCells;
		    parallelLiquidCopyCells.setSize(threadCount);

		    UT_Array<UT_Array<UT_Vector3I>> parallelBubbleBoundaryCells;
		    parallelBubbleBoundaryCells.setSize(threadCount);

		    // Set smoother cells
		    buildSmootherCells(parallelLiquidSmootherCells,
					parallelLiquidCopyCells,
					parallelBubbleBoundaryCells,
					materialCellLabels,
					visitedLiquidCells,
					cutCellWeights);

		    exint smootherListSize = 0;
		    exint copyListSize = 0;
		    exint bubbleBoundaryListSize = 0;
		    
		    for (int thread = 0; thread < threadCount; ++thread)
		    {
			smootherListSize += parallelLiquidSmootherCells[thread].size();
			copyListSize += parallelLiquidCopyCells[thread].size();
			bubbleBoundaryListSize += parallelBubbleBoundaryCells[thread].size();
		    }
		    
		    liquidSmootherCells.bumpCapacity(smootherListSize);
		    liquidCopyCells.bumpCapacity(copyListSize);
		    bubbleBoundaryCells.bumpCapacity(copyListSize);

		    for (int thread = 0; thread < threadCount; ++thread)
		    {
			liquidSmootherCells.concat(parallelLiquidSmootherCells[thread]);
			liquidCopyCells.concat(parallelLiquidCopyCells[thread]);
			bubbleBoundaryCells.concat(parallelBubbleBoundaryCells[thread]);
		    }

		    UTparallelSort(liquidSmootherCells.begin(), liquidSmootherCells.end(), cellCompare);
		    UTparallelSort(liquidCopyCells.begin(), liquidCopyCells.end(), cellCompare);
		    UTparallelSort(bubbleBoundaryCells.begin(), bubbleBoundaryCells.end(), cellCompare);
		}
	    }

	    UT_VoxelArray<SolveReal> mgSourceGrid;
	    mgSourceGrid.size(mgDomainCellLabels.getVoxelRes()[0], mgDomainCellLabels.getVoxelRes()[1], mgDomainCellLabels.getVoxelRes()[2]);

	    UT_VoxelArray<SolveReal> mgDestinationGrid;
	    mgDestinationGrid.size(mgDomainCellLabels.getVoxelRes()[0], mgDomainCellLabels.getVoxelRes()[1], mgDomainCellLabels.getVoxelRes()[2]);

	    UT_VoxelArray<SolveReal> smootherSourceGrid;
	    smootherSourceGrid.size(liquidSurface.getVoxelRes()[0], liquidSurface.getVoxelRes()[1], liquidSurface.getVoxelRes()[2]);

	    UT_VoxelArray<SolveReal> smootherDestinationGrid;
	    smootherDestinationGrid.size(liquidSurface.getVoxelRes()[0], liquidSurface.getVoxelRes()[1], liquidSurface.getVoxelRes()[2]);

	    UT_Array<SolveReal> bubbleRegionSources;
	    bubbleRegionSources.setSize(bubbleRegionCount);
	    bubbleRegionSources.constant(0);

	    UT_Array<SolveReal> bubbleRegionDestinations;
	    bubbleRegionDestinations.setSize(bubbleRegionCount);
	    bubbleRegionDestinations.constant(0);

	    UT_Array<SolveReal> tempSmootherDestinationValues;
	    tempSmootherDestinationValues.setSize(liquidSmootherCells.size());

	    // Build tile list
	    UT_Array<bool> isSmootherTileOccupied;
	    isSmootherTileOccupied.setSize(smootherDestinationGrid.numTiles());
	    isSmootherTileOccupied.constant(false);
	    
	    HDK::Utilities::findOccupiedIndexTiles(isSmootherTileOccupied,
						    liquidCopyCells,
						    materialCellLabels);

	    auto MultigridPreconditioner = [&](Vector &destinationVector, const Vector &sourceVector)
	    {
		int bubbleSmootherIterations = 10;
		int preconditionerIterations = 1;

		assert(destinationVector.rows() == sourceVector.rows() &&
			sourceVector.rows() == totalDOFCount);

		UT_StopWatch timer;
		timer.clear();
		timer.start();

		// Reset and uncompress smoother grids
		smootherSourceGrid.constant(0);
		smootherDestinationGrid.constant(0);

		bubbleRegionDestinations.constant(0);

		bubbleRegionSources.constant(0);
		for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
		{
		    if (!isBubbleRemoved[bubbleRegion])
			bubbleRegionSources[bubbleRegion] = sourceVector(liquidCellCount + bubbleRegion);
		}

		HDK::GeometricMultigridOperators::uncompressTiles(smootherSourceGrid, isSmootherTileOccupied);
		HDK::GeometricMultigridOperators::uncompressTiles(smootherDestinationGrid, isSmootherTileOccupied);

		// Copy source vector to smoother grid
		copySourceToSmootherGrid(smootherSourceGrid,
					    sourceVector,
					    liquidCellIndices,
					    liquidSmootherCells);

		time = timer.stop();
		std::cout << "  Copy source to smoother grid. Compute time: " << time << std::endl;
		timer.clear();
		timer.start();

		for (int smootherIteration = 0; smootherIteration < bubbleSmootherIterations; ++smootherIteration)
		{
		    applyBubbleSmoother(bubbleRegionDestinations,
					bubbleRegionSources,
					materialCellLabels,
					bubbleRegionIndices,
					liquidSurface,
					cutCellWeights,
					smootherDestinationGrid,
					bubbleRegionDiagonals,
					bubbleBoundaryCells,
					isBubbleRemoved,
					liquidDensity);

		    tempSmootherDestinationValues.constant(0);

		    applyLiquidSmoother(tempSmootherDestinationValues,
					materialCellLabels,
					bubbleRegionIndices,
					liquidSurface,
					cutCellWeights,
					mgDomainCellLabels,
					smootherDestinationGrid,
					smootherSourceGrid,
					bubbleRegionDestinations,
					liquidSmootherCells,
					mgExpandedOffset,
					isBubbleRemoved,
					liquidDensity);

		    copyTempDestinationToGrid(smootherDestinationGrid,
						tempSmootherDestinationValues,
						liquidSmootherCells);
		}

		time = timer.stop();
		std::cout << "  Jacobi smoothing for bubble and boudary liquid cells. Compute time: " << time << std::endl;
		timer.clear();
		timer.start();

		mgSourceGrid.constant(0);
		mgDestinationGrid.constant(0);

		// Uncompress destination within the bubble smoother
		{
		    UT_Array<bool> isMGTileOccupied;
		    isMGTileOccupied.setSize(mgDestinationGrid.numTiles());
		    isMGTileOccupied.constant(false);

		    UT_Interrupt *boss = UTgetInterrupt();

		    UTparallelForLightItems(UT_BlockedRange<exint>(0, liquidSmootherCells.size()), [&](const UT_BlockedRange<exint> &range)
		    {
			if (boss->opInterrupt())
			    return;

			for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
			{
			    UT_Vector3I cell = liquidSmootherCells[cellIndex];
			    UT_Vector3I expandedCell = cell + mgExpandedOffset;
			    int tileNumber = mgDestinationGrid.indexToLinearTile(expandedCell[0], expandedCell[1], expandedCell[2]);

			    if (!isMGTileOccupied[tileNumber])
				isMGTileOccupied[tileNumber] = true;
			}
		    });

		    HDK::GeometricMultigridOperators::uncompressTiles(mgDestinationGrid, isMGTileOccupied);
		    HDK::GeometricMultigridOperators::uncompressTiles(mgSourceGrid, isMGTileOccupied);
		}

		time = timer.stop();
		std::cout << "  Uncompress MG grids. Compute time: " << time << std::endl;
		timer.clear();
		timer.start();
	    
		copySourceToMGGrid(mgSourceGrid,
				    materialCellLabels,
				    liquidCellIndices,
				    mgDomainCellLabels,
				    sourceVector,
				    mgExpandedOffset,
				    liquidDensity);

		time = timer.stop();
		std::cout << "  Copy source to MG grid. Compute time: " << time << std::endl;
		timer.clear();
		timer.start();

		// TODO: add parameter to set preconditioner inner iterations
		for (int precondIter = 0; precondIter < preconditionerIterations; ++precondIter)
		{
		    //
		    // Apply multigrid v-cycle
		    //

		    applyBubbleSmoother(bubbleRegionDestinations,
					bubbleRegionSources,
					materialCellLabels,
					bubbleRegionIndices,
					liquidSurface,
					cutCellWeights,
					smootherDestinationGrid,
					bubbleRegionDiagonals,
					bubbleBoundaryCells,
					isBubbleRemoved,
					liquidDensity);

		    applyDirichletToMGGrid(mgSourceGrid,
					    mgDestinationGrid,
					    materialCellLabels,
					    liquidCellIndices,
					    bubbleRegionIndices,
					    liquidSurface,
					    cutCellWeights,
					    mgDomainCellLabels,
					    smootherDestinationGrid,
					    sourceVector,
					    bubbleRegionDestinations,
					    liquidSmootherCells,
					    mgExpandedOffset,
					    isBubbleRemoved,
					    liquidDensity);

		    time = timer.stop();
		    std::cout << "  Copy to MG and apply Dirichlet boundaries. Compute time: " << time << std::endl;
		    timer.clear();
		    timer.start();

		    // Apply v-cycle
		    mgPreconditioner.applyVCycle(mgDestinationGrid, mgSourceGrid, true);

		    time = timer.stop();
		    std::cout << "  Apply v-cycle time: " << time << std::endl;
		    timer.clear();
		    timer.start();
		    

		    copyMGToSmoother(smootherDestinationGrid,
					materialCellLabels,
					mgDomainCellLabels,
					mgDestinationGrid,
					liquidCopyCells,
					mgExpandedOffset);

		    time = timer.stop();
		    std::cout << "  Copy MG grid to smoother grid. Compute time: " << time << std::endl;
		    timer.clear();
		    timer.start();

		    for (int smootherIteration = 0; smootherIteration < bubbleSmootherIterations; ++smootherIteration)
		    {
			applyBubbleSmoother(bubbleRegionDestinations,
					    bubbleRegionSources,
					    materialCellLabels,
					    bubbleRegionIndices,
					    liquidSurface,
					    cutCellWeights,
					    smootherDestinationGrid,
					    bubbleRegionDiagonals,
					    bubbleBoundaryCells,
					    isBubbleRemoved,
					    liquidDensity);

			tempSmootherDestinationValues.constant(0);

			applyLiquidSmoother(tempSmootherDestinationValues,
					    materialCellLabels,
					    bubbleRegionIndices,
					    liquidSurface,
					    cutCellWeights,
					    mgDomainCellLabels,
					    smootherDestinationGrid,
					    smootherSourceGrid,
					    bubbleRegionDestinations,
					    liquidSmootherCells,
					    mgExpandedOffset,
					    isBubbleRemoved,
					    liquidDensity);

			copyTempDestinationToGrid(smootherDestinationGrid,
						    tempSmootherDestinationValues,
						    liquidSmootherCells);
		    }
		}
		    
		applyBubbleSmoother(bubbleRegionDestinations,
				    bubbleRegionSources,
				    materialCellLabels,
				    bubbleRegionIndices,
				    liquidSurface,
				    cutCellWeights,
				    smootherDestinationGrid,
				    bubbleRegionDiagonals,
				    bubbleBoundaryCells,
				    isBubbleRemoved,
				    liquidDensity);

		// Transfer bubble results back to main vector
		for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
		{
		    if (!isBubbleRemoved[bubbleRegion])
			destinationVector(liquidCellCount + bubbleRegion) = bubbleRegionDestinations[bubbleRegion];
		    else assert(bubbleRegionDestinations[bubbleRegion] == 0);
		}

		// Transfer results back to main vector
		copyMGToDestinationVector(destinationVector,
					    materialCellLabels,
					    liquidCellIndices,
					    mgDomainCellLabels,
					    mgDestinationGrid,
					    mgExpandedOffset);

		time = timer.stop();
		std::cout << "  Copy MG grid to output vector. Compute time: " << time << std::endl;
		timer.clear();
		timer.start();

		copySmootherToDestinationVector(destinationVector,
						liquidCellIndices,
						smootherDestinationGrid,
						liquidSmootherCells);

		time = timer.stop();
		std::cout << "  Copy smoother grid to output time: " << time << std::endl;
		timer.clear();
		timer.start();

		return destinationVector;
	    };

	    HDK::Utilities::solveConjugateGradient(solutionVector,
						    rhsVector,
						    [&sparseMatrix](Vector &destinationVector, const Vector &sourceVector) { destinationVector.noalias() = sparseMatrix * sourceVector; },
						    MultigridPreconditioner,				    
						    getTolerance(),
						    getMaxIterations());
	}
	else
	{
	    // Build diagonal preconditioner
	    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> diagonalPrecondMatrix(totalDOFCount, totalDOFCount);

	    {
		tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<SolveReal>>> parallelDiagonalPrecondElements;

		tbb::parallel_for(tbb::blocked_range<exint>(0, sparseMatrix.outerSize(), 1000), [&](const tbb::blocked_range<exint> &range)
		{
		    auto &localDiagonalPrecondElements = parallelDiagonalPrecondElements.local();
		    for (exint i = range.begin(); i != range.end(); ++i)
			for (typename Eigen::SparseMatrix<SolveReal, Eigen::RowMajor>::InnerIterator it(sparseMatrix, i); it; ++it)
			{
			    if (it.row() == it.col())
				localDiagonalPrecondElements.push_back(Eigen::Triplet<SolveReal>(it.row(), it.row(), 1. / it.value()));
			}
		});

		exint listSize = 0;
		parallelDiagonalPrecondElements.combine_each([&](const std::vector<Eigen::Triplet<SolveReal>> &localElements)
		{
		    listSize += localElements.size();
		});

		std::vector<Eigen::Triplet<SolveReal>> diagonalPrecondElements;
		diagonalPrecondElements.reserve(listSize);

		parallelDiagonalPrecondElements.combine_each([&](const std::vector<Eigen::Triplet<SolveReal>> &localElements)
		{
		    diagonalPrecondElements.insert(diagonalPrecondElements.end(), localElements.begin(), localElements.end());
		});

		diagonalPrecondMatrix.setFromTriplets(diagonalPrecondElements.begin(), diagonalPrecondElements.end());
		diagonalPrecondMatrix.makeCompressed();
	    }

	    HDK::Utilities::solveConjugateGradient(solutionVector,
						    rhsVector,
						    [&sparseMatrix](Vector &destinationVector, const Vector &sourceVector) { destinationVector.noalias() = sparseMatrix * sourceVector; },
						    [&diagonalPrecondMatrix](Vector &destinationVector, const Vector &sourceVector) { destinationVector.noalias() = diagonalPrecondMatrix * sourceVector; },
						    getTolerance(),
						    getMaxIterations());
	}

	Vector residual = rhsVector - sparseMatrix * solutionVector;
	std::cout << "  Re-computed Relative L2 Error: " << std::sqrt(residual.squaredNorm() / rhsVector.squaredNorm()) << std::endl;
    }

    ////////////////////////////////////////////
    //
    // Apply pressure gradient to the velocity field
    //
    ////////////////////////////////////////////

    {
	std::cout << "\n// Apply solution to pressure" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Apply solution to pressure");

	pressure->makeConstant(0);

	// Apply solution from liquid cells.
	applySolutionToPressure(*pressure, liquidCellIndices, solutionVector, 0);

	// Apply solution from bubbles.
	applySolutionToPressure(*pressure, bubbleRegionIndices, solutionVector, liquidCellCount);
    }

    {
	std::cout << "\n//  Update velocity from pressure gradient" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Update velocity from pressure gradient");

	for (int axis : {0,1,2})
	{
	    applyPressureGradient(*velocity->getField(axis),
				    *validFaces->getField(axis),
				    *pressure,
				    liquidSurface,
				    materialCellLabels,
				    bubbleRegionIndices,
				    liquidDensity,
				    axis);
	}
    }

    ////////////////////////////////////////////
    //
    // Debug check to verify that 
    // 1. each bubble region region is divergence free
    // 2. each liquid cell is divergence free
    //
    ////////////////////////////////////////////

    {
	UT_PerfMonAutoSolveEvent event(this, "Verify divergence-free constraint");
	std::cout << "\n// Verify divergence-free constraint" << std::endl;

	UT_Array<SolveReal> parallelAccumulatedLiquidDivergence;
	parallelAccumulatedLiquidDivergence.setSize(threadCount);
	parallelAccumulatedLiquidDivergence.constant(0);

	UT_Array<SolveReal> parallelMaxLiquidDivergence;
	parallelMaxLiquidDivergence.setSize(threadCount);
	parallelMaxLiquidDivergence.constant(0);

	UT_Array<SolveReal> parallelLiquidCellCount;
	parallelLiquidCellCount.setSize(threadCount);
	parallelLiquidCellCount.constant(0);

	UT_Array<UT_Array<SolveReal>> parallelBubbleRegionDivergence;
	parallelBubbleRegionDivergence.setSize(threadCount);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    parallelBubbleRegionDivergence[thread].setSize(bubbleRegionCount);
	    parallelBubbleRegionDivergence[thread].constant(0);
	}

	computeResultingDivergence(parallelAccumulatedLiquidDivergence,
				    parallelMaxLiquidDivergence,
				    parallelLiquidCellCount,
				    parallelBubbleRegionDivergence,
				    materialCellLabels,
				    bubbleRegionIndices,
				    *velocity,
				    cutCellWeights,
				    solidVelocity,
				    isBubbleTracked,
				    isBubbleRemoved);

	SolveReal accumulatedLiquidDivergence = 0;
	SolveReal maxLiquidDivergence = 0;
	SolveReal liquidCellCount = 0;

	UT_Array<SolveReal> bubbleRegionDivergence;
	bubbleRegionDivergence.setSize(bubbleRegionCount);
	bubbleRegionDivergence.constant(0);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    accumulatedLiquidDivergence += parallelAccumulatedLiquidDivergence[thread];
	    maxLiquidDivergence = std::max(maxLiquidDivergence, parallelMaxLiquidDivergence[thread]);
	    liquidCellCount += parallelLiquidCellCount[thread];

	    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    {
		if (isBubbleTracked[bubbleRegion] && !isBubbleRemoved[bubbleRegion])
		    bubbleRegionDivergence[bubbleRegion] += parallelBubbleRegionDivergence[thread][bubbleRegion];
	    }		
	}

	std::cout << "    Max liquid divergence: " << maxLiquidDivergence << std::endl;
	std::cout << "    Accumulated liquid divergence: " << accumulatedLiquidDivergence << std::endl;
	std::cout << "    Average liquid divergence: " << accumulatedLiquidDivergence / liquidCellCount << std::endl;

	std::cout << "    Bubble region divergence: " << std::endl;
	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (isBubbleTracked[bubbleRegion] && !isBubbleRemoved[bubbleRegion])
		std::cout << "      Bubble region: " << bubbleRegion << ". Divergence: " << bubbleRegionDivergence[bubbleRegion] << ". Average divergence: " << bubbleRegionDivergence[bubbleRegion] / SolveReal(bubbleRegionSizes[bubbleRegion]) << std::endl;
	}
    }

    pressureField->pubHandleModification();
    velocity->pubHandleModification();
    validFaces->pubHandleModification();

    return true;
}

void
HDK_ConstraintBubblePressureSolver::setMaterialCellLabels(SIM_RawIndexField &materialCellLabels) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(materialCellLabels.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());

	if (boss->opInterrupt())
	    return;    

	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (vit.isTileConstant())
	    {
		auto label = vit.getValue();
		
		MaterialLabels cellLabel;

		if (label == HDK::Utilities::LIQUID_CELL)
		    cellLabel = MaterialLabels::LIQUID_CELL;
		else if (label == HDK::Utilities::AIR_CELL)
		    cellLabel = MaterialLabels::BUBBLE_CELL;
		else
		{
		    assert(label == HDK::Utilities::SOLID_CELL);
		    cellLabel = MaterialLabels::SOLID_CELL;
		}

		vit.getTile()->makeConstant(cellLabel);
	    }
	    else
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    auto label = vit.getValue();
		    if (label == HDK::Utilities::LIQUID_CELL)
			vit.setValue(MaterialLabels::LIQUID_CELL);
		    else if (label == HDK::Utilities::AIR_CELL)
			vit.setValue(MaterialLabels::BUBBLE_CELL);
		    else
		    {
			assert(label == HDK::Utilities::SOLID_CELL);
			vit.setValue(MaterialLabels::SOLID_CELL);
		    }
		}
	    }
	}
    });

    materialCellLabels.fieldNC()->collapseAllTiles();
}

void
HDK_ConstraintBubblePressureSolver::buildValidFaces(SIM_VectorField &validFaces,
					    const SIM_RawIndexField &materialCellLabels,
					    const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    UT_Array<bool> isTileOccupiedList;
    for (int axis : {0,1,2})
    {
	validFaces.getField(axis)->makeConstant(HDK::Utilities::INVALID_FACE);

	isTileOccupiedList.clear();
	isTileOccupiedList.setSize(validFaces.getField(axis)->field()->numTiles());
	isTileOccupiedList.constant(false);

	HDK::Utilities::findOccupiedFaceTiles(isTileOccupiedList,
						*validFaces.getField(axis),
						materialCellLabels,
						[](const exint label){ return label == MaterialLabels::LIQUID_CELL; },
						axis);

	HDK::Utilities::uncompressTiles(*validFaces.getField(axis), isTileOccupiedList);

	HDK::Utilities::classifyValidFaces(*validFaces.getField(axis),
					    materialCellLabels,
					    *cutCellWeights[axis],
					    [](const exint label){ return label == MaterialLabels::LIQUID_CELL; },
					    axis);
    }
}

exint
HDK_ConstraintBubblePressureSolver::buildLiquidCellIndices(SIM_RawIndexField &liquidCellIndices,
							    const SIM_RawIndexField &materialCellLabels) const
{
    using SIM::FieldUtils::setFieldValue;

    liquidCellIndices.match(materialCellLabels);
    liquidCellIndices.makeConstant(HDK::Utilities::UNLABELLED_CELL);

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(materialCellLabels.field());

    UT_VoxelTileIteratorI vitt;

    exint liquidCellCount = 0;

    // Build liquid cell indices
    UT_Interrupt *boss = UTgetInterrupt();
    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
	if (boss->opInterrupt())
	    break;

	if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::LIQUID_CELL)
	{
	    vitt.setTile(vit);

	    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
	    {
		if (vitt.getValue() == MaterialLabels::LIQUID_CELL)
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
		    setFieldValue(liquidCellIndices, cell, liquidCellCount++);
		}
	    }
	}
    }

    return liquidCellCount;
}


void
HDK_ConstraintBubblePressureSolver::integrateSurfacePressure(SIM_RawField &velocity,
								const SIM_RawField &validFaces,
								const SIM_RawField &surfacePressure,
								const SIM_RawField &liquidSurface,
								const SIM_RawIndexField &liquidCellIndices,
								const SolveReal liquidDensity,
								const SolveReal dt,
								const int axis) const
{
    using SIM::FieldUtils::faceToCellMap;
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    const SolveReal invDx = 1. / velocity.getVoxelSize().maxComponent();

    UTparallelForEachNumber(validFaces.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorF vit;
	vit.setConstArray(validFaces.field());

	if (boss->opInterrupt())
	    return;

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	    {
		vit.myTileStart = tileNumber;
		vit.myTileEnd = tileNumber + 1;
		vit.rewind();

		if (!vit.isTileConstant() || vit.getValue() == HDK::Utilities::VALID_FACE)
		{
		    for (; !vit.atEnd(); vit.advance())
		    {
			if (vit.getValue() == HDK::Utilities::VALID_FACE)
			{
			    UT_Vector3I face(vit.x(), vit.y(), vit.z());

			    UT_Vector3I backwardCell = faceToCellMap(face, axis, 0);
			    UT_Vector3I forwardCell = faceToCellMap(face, axis, 1);

			    if (backwardCell[axis] >= 0 && forwardCell[axis] < liquidSurface.getVoxelRes()[axis])
			    {
				assert(getFieldValue(liquidCellIndices, backwardCell) >= 0 ||
					getFieldValue(liquidCellIndices, forwardCell) >= 0);

				if (getFieldValue(liquidCellIndices, backwardCell) == HDK::Utilities::UNLABELLED_CELL ||
				    getFieldValue(liquidCellIndices, forwardCell) == HDK::Utilities::UNLABELLED_CELL)
				{
				    SolveReal phi0 = getFieldValue(liquidSurface, backwardCell);
				    SolveReal phi1 = getFieldValue(liquidSurface, forwardCell);

				    // Compute ghost fluid weight
				    SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
				    theta = SYSclamp(theta, .01, 1.);

				    SolveReal localPressure;
				    SolveReal sign;

				    SolveReal backwardPressure = getFieldValue(surfacePressure, backwardCell);
				    SolveReal forwardPressure = getFieldValue(surfacePressure, forwardCell);

				    if (getFieldValue(liquidCellIndices, backwardCell) == HDK::Utilities::UNLABELLED_CELL)
				    {
					localPressure = SYSlerp(forwardPressure, backwardPressure, theta);
					sign = -1;
				    }
				    else
				    {
					localPressure = SYSlerp(backwardPressure, forwardPressure, theta);
					sign = 1;
				    }

				    localPressure *= invDx * dt;

				    localPressure /= (theta * SolveReal(liquidDensity));

				    setFieldValue(velocity, face, getFieldValue(velocity, face) - sign * localPressure);
				}
			    }
			}
		    }
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::computeBubbleRegionSize(UT_Array<UT_Array<exint>> &parallelBubbleRegionSizes,
						    const SIM_RawIndexField &bubbleRegionIndices) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm computeBubbleRegionSizeAlgorithm;
    computeBubbleRegionSizeAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(bubbleRegionIndices.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<exint> &localBubbleRegionSizes = parallelBubbleRegionSizes[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (vit.isTileConstant())
	    {
		exint bubbleRegion = vit.getValue();
		if (bubbleRegion >= 0)
		    localBubbleRegionSizes[bubbleRegion] += SolveReal(TILESIZE * TILESIZE * TILESIZE);
	    }
	    else
	    {
		vitt.setTile(vit);
		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    exint bubbleRegion = vitt.getValue();

		    if (bubbleRegion >= 0)
			++localBubbleRegionSizes[bubbleRegion];
		}
	    }
	}

	return 0;
    });    
}

exint
HDK_ConstraintBubblePressureSolver::buildBubbleRegionIndices(SIM_RawIndexField &bubbleRegionIndices,
								UT_Array<exint> &bubbleRegionSizes,
								const SIM_RawIndexField &materialCellLabels,
								std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    std::cout << "  Build connected bubble regions" << std::endl;

    SIM_VolumetricConnectedComponentBuilder interiorRegionBuilder(bubbleRegionIndices, materialCellLabels, cutCellWeights.data());
    
    exint bubbleRegionCount = interiorRegionBuilder.buildConnectedComponents([](const exint label)
    {
	return label == MaterialLabels::BUBBLE_CELL;
    });

    HDK::Utilities::overwriteIndices(bubbleRegionIndices,
					SIM_VolumetricConnectedComponentBuilder::INACTIVE_REGION,
					HDK::Utilities::UNLABELLED_CELL);

    std::cout << "  Compute bubble sizes for " << bubbleRegionCount << "-many bubbles." << std::endl;

    const int threadCount = UT_Thread::getNumProcessors();
    
    UT_Array<UT_Array<exint>> parallelBubbleRegionSizes;
    parallelBubbleRegionSizes.setSize(threadCount);

    for (int thread = 0; thread < threadCount; ++thread)
    {
	parallelBubbleRegionSizes[thread].setSize(bubbleRegionCount);
	parallelBubbleRegionSizes[thread].constant(0);
    }

    computeBubbleRegionSize(parallelBubbleRegionSizes, bubbleRegionIndices);

    bubbleRegionSizes.clear();
    bubbleRegionSizes.setSize(bubbleRegionCount);
    bubbleRegionSizes.constant(0);

    for (int thread = 0; thread < threadCount; ++thread)
    {
	for (int bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    bubbleRegionSizes[bubbleRegion] += parallelBubbleRegionSizes[thread][bubbleRegion];
    }

    return bubbleRegionCount;
}

exint
HDK_ConstraintBubblePressureSolver::buildCombinedRegionIndices(SIM_RawIndexField &combinedRegionIndices,
								const SIM_RawIndexField &materialCellLabels,
								std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    combinedRegionIndices.match(materialCellLabels);

    SIM_VolumetricConnectedComponentBuilder combinedRegionRegionBuilder(combinedRegionIndices, materialCellLabels, cutCellWeights.data());

    exint combinedRegionCount = combinedRegionRegionBuilder.buildConnectedComponents([](const exint label)
    {
	return label == MaterialLabels::BUBBLE_CELL || label == MaterialLabels::LIQUID_CELL;
    });

    HDK::Utilities::overwriteIndices(combinedRegionIndices,
					SIM_VolumetricConnectedComponentBuilder::INACTIVE_REGION,
					HDK::Utilities::UNLABELLED_CELL);

    return combinedRegionCount;
}

void
HDK_ConstraintBubblePressureSolver::buildRegionBooleans(UT_Array<UT_Array<bool>> &isBubbleInCombinedRegion,
							UT_Array<bool> &isLiquidInCombinedRegion,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &liquidCellIndices,
							const SIM_RawIndexField &bubbleRegionIndices,
							const SIM_RawIndexField &combinedRegionIndices,
							const exint combinedRegionCount,
							const exint bubbleRegionCount) const
{
    using SIM::FieldUtils::getFieldValue;

    isBubbleInCombinedRegion.setSize(combinedRegionCount);

    for (int combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
    {
	isBubbleInCombinedRegion[combinedRegion].setSize(bubbleRegionCount);
	isBubbleInCombinedRegion[combinedRegion].constant(false);
    }

    isLiquidInCombinedRegion.setSize(combinedRegionCount);
    isLiquidInCombinedRegion.constant(false);

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelFor(UT_BlockedRange<int>(0, combinedRegionIndices.field()->numTiles()), [&](const UT_BlockedRange<int> &range)
    {
        if (boss->opInterrupt())
	    return;

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(combinedRegionIndices.field());

	UT_Array<UT_Array<bool>> localIsBubbleInCombinedRegion;
	localIsBubbleInCombinedRegion.setSize(combinedRegionCount);

	for (int combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
	{
	    localIsBubbleInCombinedRegion[combinedRegion].setSize(bubbleRegionCount);
	    localIsBubbleInCombinedRegion[combinedRegion].constant(false);
	}

	UT_Array<bool> localIsLiquidInCombinedRegion;
	localIsLiquidInCombinedRegion.setSize(combinedRegionCount);
	localIsLiquidInCombinedRegion.constant(false);

	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (!vit.isTileConstant() || vit.getValue() >= 0)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    UT_Vector3I cell(vit.x(), vit.y(), vit.z());

		    exint combinedRegionIndex = vit.getValue();
		    if (combinedRegionIndex >= 0)
		    {
			assert(combinedRegionIndex < combinedRegionCount);

			exint bubbleRegionIndex = getFieldValue(bubbleRegionIndices, cell);
			exint liquidCellIndex = getFieldValue(liquidCellIndices, cell);

			if (bubbleRegionIndex >= 0)
			{
			    assert(liquidCellIndex == HDK::Utilities::UNLABELLED_CELL);
			    assert(bubbleRegionIndex < bubbleRegionCount);
			    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::BUBBLE_CELL);

			    localIsBubbleInCombinedRegion[combinedRegionIndex][bubbleRegionIndex] = true;
			}
			else if (liquidCellIndex >= 0)
			{
			    // Make sure that the bubble region index isn't some wrong negative number.
			    assert(bubbleRegionIndex == HDK::Utilities::UNLABELLED_CELL);
			    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::LIQUID_CELL);

			    localIsLiquidInCombinedRegion[combinedRegionIndex] = true;
			}
			else assert(false);
		    }
		    else
		    {
			assert(getFieldValue(bubbleRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL &&
				getFieldValue(liquidCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL &&
				getFieldValue(materialCellLabels, cell) == MaterialLabels::SOLID_CELL);
		    }
		}
	    }
	}

	assert(isBubbleInCombinedRegion.size() == localIsBubbleInCombinedRegion.size());
	for (int combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
	{
	    assert(isBubbleInCombinedRegion[combinedRegion].size() == localIsBubbleInCombinedRegion[combinedRegion].size());

	    for (int bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
		if (localIsBubbleInCombinedRegion[combinedRegion][bubbleRegion]) isBubbleInCombinedRegion[combinedRegion][bubbleRegion] = true;
	}

	assert(localIsLiquidInCombinedRegion.size() == isLiquidInCombinedRegion.size());
	for (int combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
	    if (localIsLiquidInCombinedRegion[combinedRegion]) isLiquidInCombinedRegion[combinedRegion] = true;
    });

#if !defined(NDEBUG)
    // Debug check. A bubble cannot exist in two regions at once.
    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
    {
	int combinedRegionCountCheck = 0;
	std::cout << "  Bubble " << bubbleRegion << " is in region: ";

	for (exint combinedRegion = 0; combinedRegion < combinedRegionCount; ++combinedRegion)
	{
	if (isBubbleInCombinedRegion[combinedRegion][bubbleRegion])
	{
	    std::cout << combinedRegion << " " << std::endl;
	    ++combinedRegionCountCheck;

	    if (isLiquidInCombinedRegion[combinedRegion] == 0)
		std::cout << "    No liquid in region with bubble" << std::endl;
	}
	}

	std::cout << std::endl;
	assert(combinedRegionCountCheck == 1);
    }
#endif

}

void
HDK_ConstraintBubblePressureSolver::trackBubbleIDs(UT_Array<bool> &isBubbleTrackedList,
					    UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
					    const GA_RWHandleI &bubbleIDHandle,
					    const GU_Detail &particles,
					    const SIM_RawIndexField &bubbleRegionIndices,
					    const SIM_RawField &liquidSurface) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::cellToCellMap;

    UT_Vector3I voxelRes = bubbleRegionIndices.getVoxelRes();

    SolveReal dx = liquidSurface.getVoxelSize().maxComponent();
    SolveReal bubbleDistance = -2 * dx;

    bool doMapBubbles = oldToNewBubbleMap.size() > 0;

    UTparallelForLightItems(GA_SplittableRange(particles.getPointRange()), [&](const GA_SplittableRange &range)
    {
    	UT_Array<bool> localIsBubbleTrackedList;
	localIsBubbleTrackedList.setSize(isBubbleTrackedList.size());
	localIsBubbleTrackedList.constant(false);

    	UT_Array<UT_Array<bool>> localOldToNewBubbleMap;

    	if (doMapBubbles)
    	{
    	    localOldToNewBubbleMap.setSize(oldToNewBubbleMap.size());

    	    for (int oldBubbleRegion = 0; oldBubbleRegion < oldToNewBubbleMap.size(); ++oldBubbleRegion)
	    {
		localOldToNewBubbleMap[oldBubbleRegion].setSize(oldToNewBubbleMap[oldBubbleRegion].size());
		localOldToNewBubbleMap[oldBubbleRegion].constant(false);
	    }
    	}

    	GA_Offset startOffset, endOffset;
    	for (GA_Iterator it(range); it.blockAdvance(startOffset, endOffset); )
    	{
    	    GA_PageNum pageNumber = GAgetPageNum(startOffset);

    	    // If the entire page is constant and unlabeled then we can skip the whole page since none of them
    	    // actually track anything.
    	    if (bubbleIDHandle.isPageConstant(pageNumber) && bubbleIDHandle.get(startOffset) == UNLABELLED_PARTICLE)
                continue;
    	    
    	    for (GA_Offset pointOffset = startOffset; pointOffset < endOffset; ++pointOffset)
    	    {
		exint bubbleID = bubbleIDHandle.get(pointOffset);

		if (bubbleID == UNLABELLED_PARTICLE)
		    continue;

		// If the particle has a valid or uninitialized bubble ID and it's close
		// to the liquid surface, check for a bubble label nearby.
		UT_Vector3 point = particles.getPos3(pointOffset);

		UT_Vector3i cell;
		liquidSurface.posToIndex(point, cell[0], cell[1], cell[2]);

                // Make sure that particle actually falls inside the grid domain
                if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 ||
                    cell[0] >= voxelRes[0] || cell[1] >= voxelRes[1] || cell[2] >= voxelRes[2])
                    continue;

    		// A particle deep in the liquid can't have an adjacent new bubble
    		// so we can safely skip it.
    		if (getFieldValue(liquidSurface, cell) < bubbleDistance)
    		    continue;

    		exint newBubbleRegion = getFieldValue(bubbleRegionIndices, cell);
    		if (newBubbleRegion >= 0)
    		{
    		    localIsBubbleTrackedList[newBubbleRegion] = true;

    		    // bubble ID can be uninitialized. If so, we can't map anything.
    		    if (doMapBubbles && bubbleID >= 0) 
			localOldToNewBubbleMap[bubbleID][newBubbleRegion] = true;
    		}

    		for (int axis : {0,1,2})
    		    for (int direction : {0,1})
    		    {
			UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

			if (adjacentCell[axis] < 0 || adjacentCell[axis] >= voxelRes[axis])
			    continue;

			exint newBubbleRegion = getFieldValue(bubbleRegionIndices, adjacentCell);

			if (newBubbleRegion >= 0)
			{
			    localIsBubbleTrackedList[newBubbleRegion] = true;

			    if (doMapBubbles && bubbleID >= 0)
				localOldToNewBubbleMap[bubbleID][newBubbleRegion] = true;
			}
    		    }
    	    }
    	}

    	for (exint newBubbleRegion = 0; newBubbleRegion < isBubbleTrackedList.size(); ++newBubbleRegion)
    	{
    	    if (localIsBubbleTrackedList[newBubbleRegion])
                isBubbleTrackedList[newBubbleRegion] = true;
    	}

	if (doMapBubbles)
	{
	    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldToNewBubbleMap.size(); ++oldBubbleRegion)
	    {
		assert(oldToNewBubbleMap[oldBubbleRegion].size() == isBubbleTrackedList.size());
		for (int newBubbleRegion = 0; newBubbleRegion < oldToNewBubbleMap[oldBubbleRegion].size(); ++newBubbleRegion)
		{
		    if (localOldToNewBubbleMap[oldBubbleRegion][newBubbleRegion])
			oldToNewBubbleMap[oldBubbleRegion][newBubbleRegion] = true;
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::updateBubbleIDs(GA_RWHandleI &bubbleIDHandle,
					    const UT_Array<bool> &isBubbleTrackedList,
					    const GU_Detail &particles,
					    const SIM_RawIndexField &bubbleRegionIndices,
					    const SIM_RawField &liquidSurface) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::cellToCellMap;

    UT_Vector3I voxelRes = bubbleRegionIndices.getVoxelRes();

    SolveReal dx = liquidSurface.getVoxelSize().maxComponent();
    SolveReal bubbleDistance = -2 * dx;

    UTparallelFor(GA_SplittableRange(particles.getPointRange()), [&](const GA_SplittableRange &range)
    {
	GA_Offset startOffset, endOffset;
	for (GA_Iterator it(range); it.blockAdvance(startOffset, endOffset); )
	{
	    for (GA_Offset pointOffset = startOffset; pointOffset < endOffset; ++pointOffset)
	    {
		UT_Vector3 point = particles.getPos3(pointOffset);
		exint bubbleID = bubbleIDHandle.get(pointOffset);

		UT_Vector3i cell;
		liquidSurface.posToIndex(point, cell[0], cell[1], cell[2]);

		if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 ||
		    cell[0] >= voxelRes[0] || cell[1] >= voxelRes[1] || cell[2] >= voxelRes[2])
		{
		    if (bubbleID != UNLABELLED_PARTICLE)
			bubbleIDHandle.set(pointOffset, UNLABELLED_PARTICLE);
		    continue;
		}
		else if (getFieldValue(liquidSurface, cell) < bubbleDistance)
		{
		    if (bubbleID != UNLABELLED_PARTICLE)
			bubbleIDHandle.set(pointOffset, UNLABELLED_PARTICLE);
		    continue;
		}

		// If the particle is near a bubble, we want to associate the particle with the bubble.
		// Ideally we would just pick nearest bubble but that could be slower than desired.
		// If the bubble is too small to have other particles map to it then it is probably not worth tracking anyway.

		exint newBubbleID = UNLABELLED_PARTICLE;
		exint bubbleRegion = getFieldValue(bubbleRegionIndices, cell);

		// If the bubble region isn't currently being tracked, we want it to remain untracked.
		if (bubbleRegion >= 0 && isBubbleTrackedList[bubbleRegion])
		    newBubbleID = bubbleRegion;
		else
		{
		    bool bubbleFound = false;

		    for (int axis = 0; axis < 3 && !bubbleFound; ++axis)
			for (int direction : {0,1})
			{
			    UT_Vector3I adjacentCell = cellToCellMap(UT_Vector3I(cell), axis, direction);

			    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= voxelRes[axis])
				continue;

			    exint bubbleRegion = getFieldValue(bubbleRegionIndices, adjacentCell);

			    if (bubbleRegion >= 0 && isBubbleTrackedList[bubbleRegion])
			    {
				newBubbleID = bubbleRegion;
				bubbleFound = true;
				break;	
			    }
			}
		}

		if (newBubbleID != UNLABELLED_PARTICLE || bubbleID != UNLABELLED_PARTICLE)
		    bubbleIDHandle.set(pointOffset, newBubbleID);
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::distributeOldBubbleVolumes(UT_Array<SolveReal> &newBubbleVolumes,
								UT_Array<bool> &isBubbleUninitializedList,
								const GA_RWHandleF &bubbleVolumeHandle,
								const UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
								const UT_Array<bool> &isBubbleTrackedList,
								const UT_Array<exint> &newBubbleRegionSizes) const
{
    constexpr exint UNVISITED_CELL = -2;
    constexpr exint VISITED_CELL = -1;

    std::cout << "    Distribute old bubble masses" << std::endl;

    // In order to update the new bubble volumes/masses based on the old bubble volumes,
    // we need to build a mapping between the two. We can handle this with a bi-partite graph.
    
    exint oldBubbleRegionCount = oldToNewBubbleMap.size();
    exint newBubbleRegionCount = newBubbleRegionSizes.size();

    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldBubbleRegionCount; ++oldBubbleRegion)
	assert(newBubbleRegionCount == oldToNewBubbleMap[oldBubbleRegion].size());

    UT_Array<exint> oldBubbleComponents;
    oldBubbleComponents.setSize(oldBubbleRegionCount);
    oldBubbleComponents.constant(UNVISITED_CELL);

    UT_Array<exint> newBubbleComponents;
    newBubbleComponents.setSize(newBubbleRegionCount);
    newBubbleComponents.constant(UNVISITED_CELL);

    UT_Array<exint> oldBubblesToVisit;

    exint connectedComponentCount = 0;

    // Build connected components
    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldBubbleRegionCount; ++oldBubbleRegion)
    {
	// Start at an unvisited old bubble region to begin building the connected components in the bipartite graph.
	if (oldBubbleComponents[oldBubbleRegion] == UNVISITED_CELL)
	{
	    oldBubblesToVisit.clear();
	    oldBubblesToVisit.append(oldBubbleRegion);

	    oldBubbleComponents[oldBubbleRegion] = VISITED_CELL;

	    int oldToNewEdgeCount = 0;
	    while(!oldBubblesToVisit.isEmpty())
	    {
		exint localOldBubbleRegion = oldBubblesToVisit.last();
		oldBubblesToVisit.removeLast();

		assert(oldBubbleComponents[localOldBubbleRegion] == VISITED_CELL);

		oldBubbleComponents[localOldBubbleRegion] = connectedComponentCount;

		// Sweep through the set of new bubble regions and check for an edge from the old bubble to
		// the current new bubble. If an edge exists, continue with search.
		for (exint newBubbleRegion = 0; newBubbleRegion < newBubbleRegionCount; ++newBubbleRegion)
		{
		    if (oldToNewBubbleMap[localOldBubbleRegion][newBubbleRegion])
		    {
			assert(isBubbleTrackedList[newBubbleRegion]);

			if (newBubbleComponents[newBubbleRegion] == UNVISITED_CELL)
			{
			    // Loop over the set of old bubble and check for an edge from this new bubble back to an old
			    // bubble. If that old bubble has not yet been visited, add it to the queue.
			    for (exint adjacentOldBubbleRegion = 0; adjacentOldBubbleRegion < oldBubbleRegionCount; ++adjacentOldBubbleRegion)
			    {
				if (oldToNewBubbleMap[adjacentOldBubbleRegion][newBubbleRegion])
				{
				    if (oldBubbleComponents[adjacentOldBubbleRegion] == UNVISITED_CELL)
				    {
					oldBubblesToVisit.append(adjacentOldBubbleRegion);
					oldBubbleComponents[adjacentOldBubbleRegion] = VISITED_CELL;
				    }
				    else if (oldBubbleComponents[adjacentOldBubbleRegion] != VISITED_CELL)
					assert(oldBubbleComponents[adjacentOldBubbleRegion] == connectedComponentCount);
				}
			    }

			    newBubbleComponents[newBubbleRegion] = connectedComponentCount;
			    ++oldToNewEdgeCount;
			}
			else assert(newBubbleComponents[newBubbleRegion] == connectedComponentCount);
		    }
		}
	    }

	    if (oldToNewEdgeCount == 0)
		std::cout << "Old bubble: " << oldBubbleRegion << " has been abandoned. Volume: " << bubbleVolumeHandle.get(oldBubbleRegion);
	    
	    ++connectedComponentCount;
	}
    }

    // Build total volumes for old and new bubbles in each connected component
    UT_Array<SolveReal> combinedOldVolumes;
    combinedOldVolumes.setSize(connectedComponentCount);
    combinedOldVolumes.constant(0);

    UT_Array<SolveReal> combinedNewVolumes;
    combinedNewVolumes.setSize(connectedComponentCount);
    combinedNewVolumes.constant(0);

    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldBubbleRegionCount; ++oldBubbleRegion)
    {
	exint localComponent = oldBubbleComponents[oldBubbleRegion];

	assert(localComponent >= 0);

	combinedOldVolumes[localComponent] += bubbleVolumeHandle.get(oldBubbleRegion);
    }

    for (exint newBubbleRegion = 0; newBubbleRegion < newBubbleRegionCount; ++newBubbleRegion)
    {
	exint localComponent = newBubbleComponents[newBubbleRegion];

	if (localComponent >= 0)
	    combinedNewVolumes[localComponent] += SolveReal(newBubbleRegionSizes[newBubbleRegion]);
	else
	{
	    assert(localComponent == UNVISITED_CELL);

	    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldBubbleRegionCount; ++oldBubbleRegion)
		assert(!oldToNewBubbleMap[oldBubbleRegion][newBubbleRegion]);
	}
    }

    // Distribute old rest volumes to new bubbles as a fraction of total new bubble volume.
    newBubbleVolumes.setSize(newBubbleRegionCount);
    newBubbleVolumes.constant(0);

    for (exint newBubbleRegion = 0; newBubbleRegion < newBubbleRegionCount; ++newBubbleRegion)
    {
	exint localComponent = newBubbleComponents[newBubbleRegion];

	if (localComponent >= 0)
	{
	    SolveReal newBubbleVolume = SolveReal(newBubbleRegionSizes[newBubbleRegion]);

	    assert(combinedNewVolumes[localComponent] > 0);
	    newBubbleVolumes[newBubbleRegion] = combinedOldVolumes[localComponent] * newBubbleVolume / combinedNewVolumes[localComponent]; 
	}
	else
	{
	    std::cout << "    New bubble region: " << newBubbleRegion << " has no old bubble connection and cannot preserve rest volume" << std::endl;
	    isBubbleUninitializedList[newBubbleRegion] = true;
	}
    }
}

void
HDK_ConstraintBubblePressureSolver::removeBubblesAtOpenBoundaries(UT_Array<bool> &isBubbleRemoved,
							    const SIM_RawIndexField &bubbleRegionIndices,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(bubbleRegionIndices.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(bubbleRegionIndices.field());
        
	if (boss->opInterrupt())
	    return;
	
	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (!vit.isTileConstant() || vit.getValue() >= 0)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    exint bubbleRegion = vit.getValue();
		    if (bubbleRegion >= 0 && !isBubbleRemoved[bubbleRegion])
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				if (getFieldValue(*cutCellWeights[axis], face) > 0)
				{
				    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= bubbleRegionIndices.getVoxelRes()[axis])
					isBubbleRemoved[bubbleRegion] = true;
				}
			    }
		    }
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::buildLiquidRHS(Vector &rhsVector,
					    const SIM_RawIndexField &materialCellLabels,
					    const SIM_RawIndexField &liquidCellIndices,
					    const SIM_RawIndexField &bubbleRegionIndices,
					    const SIM_RawIndexField &combinedRegionIndices,
					    const SIM_VectorField &velocity,
					    const SIM_VectorField *solidVelocity,
					    const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::cellToFaceMap;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(liquidCellIndices.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(liquidCellIndices.field());
        
	if (boss->opInterrupt())
	    return;    

	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (!vit.isTileConstant())
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    exint liquidIndex = vit.getValue();
		    if (liquidIndex >= 0)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::LIQUID_CELL);
			assert(getFieldValue(bubbleRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			assert(getFieldValue(combinedRegionIndices, cell) >= 0);

			SolveReal divergence = 0;

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				SolveReal sign = (direction % 2 == 0) ? 1. : -1.;
				SolveReal weight = getFieldValue(*cutCellWeights[axis], face);
		    
				if (weight > 0)
				    divergence += sign * weight * getFieldValue(*velocity.getField(axis), face);

				if (solidVelocity != nullptr && weight < 1.)
				{
				    UT_Vector3 point;
				    velocity.getField(axis)->indexToPos(face[0], face[1], face[2], point);

				    divergence += sign * (1. - weight) * solidVelocity->getField(axis)->getValue(point);
				}
			    }

			rhsVector(liquidIndex) = divergence;
		    }
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::buildBubbleRHS(UT_Array<UT_Array<SolveReal>> &parallelBubbleRegionDivergence,
					    const SIM_RawIndexField &materialCellLabels,
					    const SIM_RawIndexField &liquidCellIndices,
					    const SIM_RawIndexField &bubbleRegionIndices,
					    const SIM_RawIndexField &combinedRegionIndices,
					    const std::array<const SIM_RawField *, 3> &cutCellWeights,
					    const UT_Array<bool> &isBubbleRemoved,
					    const SIM_VectorField &velocity,
					    const SIM_VectorField *solidVelocity,
					    const exint bubbleRegionCount) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::cellToCellMap;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildBubbleRHSAlgorithm;
    buildBubbleRHSAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(bubbleRegionIndices.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<SolveReal> &localBubbleRegionDivergence = parallelBubbleRegionDivergence[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() || vit.getValue() >= 0)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    const exint bubbleRegion = vitt.getValue();

		    if (bubbleRegion >= 0 && !isBubbleRemoved[bubbleRegion])
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::BUBBLE_CELL);
			assert(getFieldValue(liquidCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			assert(getFieldValue(combinedRegionIndices, cell) >= 0);

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				// We're assuming a solid static boundary so there is no right-hand-side velocity.
				if (adjacentCell[axis] < 0 || adjacentCell[axis] >= bubbleRegionIndices.getVoxelRes()[axis])
				    continue;

				SolveReal sign = (direction == 0) ? 1. : -1.;
				UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				SolveReal weight = getFieldValue(*cutCellWeights[axis], face);

				int adjacentLiquidIndex = getFieldValue(liquidCellIndices, adjacentCell);

				if (adjacentLiquidIndex >= 0 && weight > 0.)
				{
				    assert(getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::LIQUID_CELL);
				    assert(getFieldValue(bubbleRegionIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);
				    assert(getFieldValue(combinedRegionIndices, adjacentCell) >= 0);

				    localBubbleRegionDivergence[bubbleRegion] += sign * weight * getFieldValue(*velocity.getField(axis), face);
				}

				if (solidVelocity != nullptr && weight < 1.)
				{
				    UT_Vector3 point;
				    velocity.getField(axis)->indexToPos(face[0], face[1], face[2], point);

				    localBubbleRegionDivergence[bubbleRegion] += sign * (1. - weight) * solidVelocity->getField(axis)->getValue(point);
				}
			    }
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_ConstraintBubblePressureSolver::buildLiquidRegionDivergence(UT_Array<UT_Array<SolveReal>> &parallelLiquidRegionDivergence,
							UT_Array<UT_Array<exint>> &parallelLiquidRegionCellCount,
							const SIM_RawIndexField &combinedRegionIndices,
							const SIM_RawIndexField &liquidCellIndices,
							const UT_Array<bool> &doesRegionHaveRemovedBubble,
							const Vector &rhsVector) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildLiquidRegionDivergenceAlgorithm;
    buildLiquidRegionDivergenceAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(liquidCellIndices.field());
	vit.splitByTile(info);
        
	UT_VoxelTileIteratorI vitt;

	UT_Array<SolveReal> &localLiquidRegionDivergence = parallelLiquidRegionDivergence[info.job()];
	UT_Array<exint> &localLiquidRegionCellCount = parallelLiquidRegionCellCount[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant())
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    exint liquidIndex = vitt.getValue();
		    if (liquidIndex >= 0)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			exint combinedRegion = getFieldValue(combinedRegionIndices, cell);
			assert(combinedRegion >= 0);

			if (!doesRegionHaveRemovedBubble[combinedRegion])
			{
			    localLiquidRegionDivergence[combinedRegion] += rhsVector(liquidIndex);
			    ++localLiquidRegionCellCount[combinedRegion];
			}
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_ConstraintBubblePressureSolver::removeDivergenceFromLiquidCells(Vector &rhsVector,
							    const UT_Array<SolveReal> &regionDivergence,
							    const SIM_RawIndexField &combinedRegionIndices,
							    const SIM_RawIndexField &liquidCellIndices,
							    const UT_Array<bool> &doesRegionHaveRemovedBubble) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(liquidCellIndices.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(liquidCellIndices.field());
        
	if (boss->opInterrupt())
	    return;    

	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (!vit.isTileConstant())
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    exint liquidIndex = vit.getValue();
		    if (liquidIndex >= 0)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());

			exint combinedRegion = getFieldValue(combinedRegionIndices, cell);
			assert(combinedRegion >= 0);

			if (!doesRegionHaveRemovedBubble[combinedRegion])
			    rhsVector(liquidIndex) -= regionDivergence[combinedRegion];
		    }
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::buildLiquidRows(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelSparseElements,
						    const SIM_RawIndexField &liquidCellIndices,
						    const SIM_RawIndexField &bubbleRegionIndices,
						    const SIM_RawField &liquidSurface,
						    const std::array<const SIM_RawField *, 3> &cutCellWeights,
						    const UT_Array<bool> &isBubbleRemoved,
						    const SolveReal liquidDensity,
						    const exint liquidCellCount) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildLiquidRowsAlgorithm;
    buildLiquidRowsAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(liquidCellIndices.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	std::vector<Eigen::Triplet<SolveReal>> &localSparseElements = parallelSparseElements[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;;

	    if (!vit.isTileConstant())
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    exint liquidIndex = vitt.getValue();
		    if (liquidIndex >= 0)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
			// Build non-zeros for each liquid voxel row
			SolveReal diagonal = 0.;

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				if (adjacentCell[axis] < 0 || adjacentCell[axis] >= liquidCellIndices.getVoxelRes()[axis])
				    continue;

				UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				SolveReal weight = getFieldValue(*cutCellWeights[axis], face);
				if (weight > 0)
				{
				    exint adjacentLiquidIndex = getFieldValue(liquidCellIndices, adjacentCell);

				    if (adjacentLiquidIndex >= 0)
				    {
					weight /= liquidDensity;
					diagonal += weight;

					localSparseElements.push_back(Eigen::Triplet<SolveReal>(liquidIndex, adjacentLiquidIndex, -weight));
					assert(getFieldValue(bubbleRegionIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);
				    }
				    else
				    {
					SolveReal phi0 = getFieldValue(liquidSurface, cell);
					SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

					SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
					theta = SYSclamp(theta, .01, 1.);

					SolveReal localDensity = theta * liquidDensity;

					exint bubbleRegion = getFieldValue(bubbleRegionIndices, adjacentCell);
					assert(bubbleRegion >= 0);

					weight /= localDensity;

					if (!isBubbleRemoved[bubbleRegion])
					    localSparseElements.push_back(Eigen::Triplet<SolveReal>(liquidIndex, bubbleRegion + liquidCellCount, -weight));
				
					diagonal += weight;
				    }
				}
			    }

			    assert(diagonal > 0.);
			    localSparseElements.push_back(Eigen::Triplet<SolveReal>(liquidIndex, liquidIndex, diagonal));
		    }
		}
	    }
	}

	return 0;
    });    
}

void
HDK_ConstraintBubblePressureSolver::buildBubbleRows(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelSparseElements,
						    UT_Array<UT_Array<SolveReal>> &parallelBubbleRegionDiagonals,
						    const SIM_RawIndexField &liquidCellIndices,
						    const SIM_RawIndexField &bubbleRegionIndices,
						    const SIM_RawField &liquidSurface,
						    const std::array<const SIM_RawField *, 3> &cutCellWeights,
						    const UT_Array<bool> &isBubbleRemoved,
						    const SolveReal liquidDensity,
						    const exint liquidCellCount) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;

    UT_ThreadedAlgorithm buildBubbleRowsAlgorithm;
    buildBubbleRowsAlgorithm.run([&](const UT_JobInfo &info)
    {
	std::vector<Eigen::Triplet<SolveReal>> &localSparseElements = parallelSparseElements[info.job()];
	UT_Array<SolveReal> &localBubbleRegionDiagonals = parallelBubbleRegionDiagonals[info.job()];

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(bubbleRegionIndices.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Interrupt *boss = UTgetInterrupt();

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;;

	    if (!vit.isTileConstant() || (vit.getValue() >= 0 && !isBubbleRemoved[vit.getValue()]))
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    exint bubbleRegion = vitt.getValue();

		    if (bubbleRegion >= 0 && !isBubbleRemoved[bubbleRegion])
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			assert(getFieldValue(liquidCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);

		    	// Build air-liquid pressure gradient entries.
			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				if (adjacentCell[axis] < 0 || adjacentCell[axis] >= bubbleRegionIndices.getVoxelRes()[axis])
				    continue;

				UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				SolveReal weight = getFieldValue(*cutCellWeights[axis], face);
	                    	
				if (weight > 0.)
				{
				    // Build row entries and liquid-air flux.
				    exint adjacentLiquidIndex = getFieldValue(liquidCellIndices, adjacentCell);

				    if (adjacentLiquidIndex >= 0)
				    {
					assert(getFieldValue(bubbleRegionIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);

					SolveReal phi0 = getFieldValue(liquidSurface, cell);
					SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

					assert(phi0 > 0);

					SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
					theta = SYSclamp(theta, .01, 1.);

					SolveReal localDensity = theta * liquidDensity;

					weight /= localDensity;

					localBubbleRegionDiagonals[bubbleRegion] += weight;

					localSparseElements.push_back(Eigen::Triplet<SolveReal>(bubbleRegion + liquidCellCount, adjacentLiquidIndex, -weight));
				    }
				    else assert(getFieldValue(bubbleRegionIndices, adjacentCell) == bubbleRegion);
				}
			    }
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_ConstraintBubblePressureSolver::computeBubblePressure(UT_Array<UT_Array<SolveReal>> &parallelAccumulatedPressure,
						    const SIM_RawField &pressure,
						    const SIM_RawIndexField &bubbleRegionIndices,
						    const UT_Array<bool> &isBubbleRemoved,
						    const exint bubbleRegionCount) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm computeBubblePressureAlgorithm;
    computeBubblePressureAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(bubbleRegionIndices.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<SolveReal> &localAccumulatedPressure = parallelAccumulatedPressure[info.job()];
    
	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;
	    
	    if (!vit.isTileConstant() || (vit.getValue() >= 0 && !isBubbleRemoved[vit.getValue()]))
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    exint bubbleRegion = vitt.getValue();
		    if (bubbleRegion >= 0 && !isBubbleRemoved[bubbleRegion])
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
			localAccumulatedPressure[bubbleRegion] += getFieldValue(pressure, cell);
		    }
		}
	    }
	}

	return 0;
    });
}

//
// Multigrid operations
//

void
HDK_ConstraintBubblePressureSolver::buildMGDomainLabels(UT_VoxelArray<int> &mgDomainCellLabels,
						const SIM_RawIndexField &materialCellLabels,
						const SIM_RawIndexField &liquidCellIndices) const
{
    using SIM::FieldUtils::getFieldValue;
    
    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(materialCellLabels.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());

	if (boss->opInterrupt())
	    return;    

	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    // TODO: handle constant tiles
	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::LIQUID_CELL ||
		vit.getValue() == MaterialLabels::BUBBLE_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    UT_Vector3I cell(vit.x(), vit.y(), vit.z());

		    if (vit.getValue() == MaterialLabels::LIQUID_CELL)
		    {
			assert(getFieldValue(liquidCellIndices, cell) >= 0);			    
			mgDomainCellLabels.setValue(cell, MGCellLabels::INTERIOR_CELL);
		    }
		    else if (vit.getValue() == MaterialLabels::BUBBLE_CELL)
		    {
			assert(getFieldValue(liquidCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			mgDomainCellLabels.setValue(cell, MGCellLabels::DIRICHLET_CELL);
		    }
		    else
		    {
			assert(mgDomainCellLabels(cell) == MGCellLabels::EXTERIOR_CELL);
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::SOLID_CELL);
			assert(getFieldValue(liquidCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
		    }
		}
	    }
	}
    });

    mgDomainCellLabels.collapseAllTiles();
}

void
HDK_ConstraintBubblePressureSolver::buildMGBoundaryWeights(UT_VoxelArray<SolveReal> &boundaryWeights,
							    const SIM_RawField &validFaces,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &bubbleRegionIndices,
							    const UT_VoxelArray<int> &domainCellLabels,
							    const SIM_RawField &liquidSurface,
							    const SIM_RawField &cutCellWeights,
							    const SolveReal liquidDensity,
							    const int axis) const
{
    using SIM::FieldUtils::faceToCellMap;
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(validFaces.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorF vit;
	vit.setConstArray(validFaces.field());

	if (boss->opInterrupt())
	    return;
     
	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (!vit.isTileConstant() ||
		vit.getValue() == HDK::Utilities::VALID_FACE)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == HDK::Utilities::VALID_FACE)
		    {
			UT_Vector3I face(vit.x(), vit.y(), vit.z());

			SolveReal weight = getFieldValue(cutCellWeights, face);
			assert(weight > 0);

			UT_Vector3I backwardCell = faceToCellMap(face, axis, 0);
			UT_Vector3I forwardCell = faceToCellMap(face, axis, 1);

			auto backwardCellLabel = domainCellLabels(backwardCell);
			auto forwardCellLabel = domainCellLabels(forwardCell);

			if ((backwardCellLabel == MGCellLabels::INTERIOR_CELL && forwardCellLabel == MGCellLabels::DIRICHLET_CELL) ||
			    (backwardCellLabel == MGCellLabels::DIRICHLET_CELL && forwardCellLabel == MGCellLabels::INTERIOR_CELL))
			{
			    assert(getFieldValue(materialCellLabels, backwardCell) == MaterialLabels::BUBBLE_CELL ||
				    getFieldValue(materialCellLabels, forwardCell) == MaterialLabels::BUBBLE_CELL);

			    SolveReal phi0 = getFieldValue(liquidSurface, backwardCell);
			    SolveReal phi1 = getFieldValue(liquidSurface, forwardCell);

			    SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
			    theta = SYSclamp(theta, .01, 1.);

			    SolveReal interpolatedDensity = theta * liquidDensity;

			    weight *= liquidDensity / interpolatedDensity;
			}

			boundaryWeights.setValue(face, weight);
		    }
		}
	    }
	}
    });
}

bool
HDK_ConstraintBubblePressureSolver::unitTestMGLabels(const UT_VoxelArray<int> &mgDomainCellLabels,
						const SIM_RawIndexField &materialCellLabels,
						const SIM_RawIndexField &liquidCellIndices,
						const SIM_RawIndexField &bubbleRegionIndices,
						const UT_Vector3I &mgExpandedOffset) const
{
    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    bool passedTest = true;

    UTparallelForEachNumber(materialCellLabels.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());

	if (boss->opInterrupt())
	    return;
 
	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    for (; !vit.atEnd(); vit.advance())
	    {
		if (!passedTest)
		    return;

		UT_Vector3I cell(vit.x(), vit.y(), vit.z());
		UT_Vector3I expandedCell = cell + mgExpandedOffset;

		if (vit.getValue() == MaterialLabels::LIQUID_CELL)
		{
		    if (mgDomainCellLabels(expandedCell) != MGCellLabels::INTERIOR_CELL &&
			mgDomainCellLabels(expandedCell) != MGCellLabels::BOUNDARY_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(liquidCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(bubbleRegionIndices, cell) != HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		}
		else if (vit.getValue() == MaterialLabels::BUBBLE_CELL)
		{
		    if (mgDomainCellLabels(expandedCell) != MGCellLabels::DIRICHLET_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(liquidCellIndices, cell) != HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(bubbleRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		}
		else
		{
		    if (vit.getValue() != MaterialLabels::SOLID_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (mgDomainCellLabels(expandedCell) != MGCellLabels::EXTERIOR_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(liquidCellIndices, cell) != HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(bubbleRegionIndices, cell) != HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		}
	    }
	}
    });

    return passedTest;
}

void
HDK_ConstraintBubblePressureSolver::buildInitialLiquidSmootherLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidBandCells,
							    const SIM_RawIndexField &materialCellLabels,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInitialLiquidSmootherLayerAlgorithm;
    buildInitialLiquidSmootherLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<UT_Vector3I> &localLiquidBandCells = parallelLiquidBandCells[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::LIQUID_CELL)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			bool isBoundaryCell = false;
			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
				    continue;

				UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				if (getFieldValue(*cutCellWeights[axis], face) > 0 && 
				    getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::BUBBLE_CELL)
				    isBoundaryCell = true;
			}

			if (isBoundaryCell)
			    localLiquidBandCells.append(cell);
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_ConstraintBubblePressureSolver::setVisitedLiquidCells(SIM_RawIndexField &visitedLiquidCells,
						    const UT_Array<UT_Vector3I> &newLiquidCells) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForLightItems(UT_BlockedRange<exint>(0, newLiquidCells.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;

	exint startIndex = range.begin();
    
	if (startIndex > 0)
	{
	    while (startIndex != range.end() && newLiquidCells[startIndex] == newLiquidCells[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	for (exint cellIndex = startIndex; cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = newLiquidCells[cellIndex];

	    if (cell == oldCell)
		continue;
	
	    oldCell = cell;

	    assert(getFieldValue(visitedLiquidCells, cell) == UNVISITED_CELL);
	    setFieldValue(visitedLiquidCells, cell, VISITED_CELL);
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::buildNextLiquidSmootherLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidLayerCells,
							    const UT_Array<UT_Vector3I> &oldLiquidLayerCells,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &visitedLiquidCells,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildNextLiquidSmootherLayerAlgorithm;
    buildNextLiquidSmootherLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_Array<UT_Vector3I> &localLiquidLayerCells = parallelLiquidLayerCells[info.job()];

	exint startIndex, endIndex;
	const exint elementSize = oldLiquidLayerCells.entries();
	info.divideWork(elementSize, startIndex, endIndex);

	if (boss->opInterrupt())
	    return 0;
    
	if (startIndex > 0)
	{
	    while (startIndex != endIndex && oldLiquidLayerCells[startIndex - 1] == oldLiquidLayerCells[startIndex])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	const exint localEndIndex = endIndex;
	for (exint i = startIndex; i < localEndIndex; ++i)
	{
	    UT_Vector3I cell = oldLiquidLayerCells[i];

	    if (oldCell == cell)
		continue;

	    oldCell = cell;

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::LIQUID_CELL);
	    assert(getFieldValue(visitedLiquidCells, cell) == VISITED_CELL);

	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
	    
		    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
			continue;

		    UT_Vector3I face = cellToFaceMap(cell, axis, direction);
		    if (getFieldValue(*cutCellWeights[axis], face) > 0)
		    {
			if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::LIQUID_CELL &&
			    getFieldValue(visitedLiquidCells, adjacentCell) == UNVISITED_CELL)
			    localLiquidLayerCells.append(adjacentCell);
		    }
		}
	}

	return 0;
    });
}

void
HDK_ConstraintBubblePressureSolver::buildSmootherCells(UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidSmootherCells,
						UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidCopyCells,
						UT_Array<UT_Array<UT_Vector3I>> &parallelBubbleBoundaryCells,
						const SIM_RawIndexField &materialCellLabels,
						const SIM_RawIndexField &visitedLiquidCells,
						const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildBubbleSmootherCellsAlgorithm;
    buildBubbleSmootherCellsAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<UT_Vector3I> &localLiquidSmootherCells = parallelLiquidSmootherCells[info.job()];
	UT_Array<UT_Vector3I> &localLiquidCopyCells = parallelLiquidCopyCells[info.job()];
	UT_Array<UT_Vector3I> &localBubbleBoundaryCells = parallelBubbleBoundaryCells[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
		    if (vitt.getValue() == MaterialLabels::LIQUID_CELL)
		    {
			if (getFieldValue(visitedLiquidCells, cell) == VISITED_CELL)
			{
			    localLiquidSmootherCells.append(cell);
			    localLiquidCopyCells.append(cell);

			    bool isBubbleBoundaryCell = false;
			    for (int axis : {0,1,2})
				for (int direction : {0,1})
				{
				    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
					continue;

				    UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				    if (getFieldValue(*cutCellWeights[axis], face) > 0)
				    {
					if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::BUBBLE_CELL)
					{
					    assert(getFieldValue(visitedLiquidCells, adjacentCell) == UNVISITED_CELL);
					    isBubbleBoundaryCell = true;
					}
				    }
				}

			    if (isBubbleBoundaryCell)
				localBubbleBoundaryCells.append(cell);
			}
			else
			{
			    assert(getFieldValue(visitedLiquidCells, cell) == UNVISITED_CELL);
			    
			    // Add adjacent liquid cells to copy list
			    bool isSmootherBoundaryCell = false;
			    for (int axis : {0,1,2})
				for (int direction : {0,1})
				{
				    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
					continue;

				    UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				    if (getFieldValue(*cutCellWeights[axis], face) > 0)
				    {
					assert(getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::LIQUID_CELL);
					
					if (getFieldValue(visitedLiquidCells, adjacentCell) == VISITED_CELL)
					    isSmootherBoundaryCell = true;
				    }
				}

			    if (isSmootherBoundaryCell)
				localLiquidCopyCells.append(cell);
			}			
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_ConstraintBubblePressureSolver::copySourceToSmootherGrid(UT_VoxelArray<SolveReal> &smootherSourceGrid,
							const Vector &sourceVector,
							const SIM_RawIndexField &liquidCellIndices,
							const UT_Array<UT_Vector3I> &bubbleSmootherCells) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForLightItems(UT_BlockedRange<exint>(0, bubbleSmootherCells.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;

	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = bubbleSmootherCells[cellIndex];
	    exint index = getFieldValue(liquidCellIndices, cell);
	    assert(index >= 0);
	    smootherSourceGrid.setValue(cell, sourceVector(index));
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::buildBubbleLaplacian(UT_Array<UT_Array<SolveReal>> &parallelBubbleRegionLaplacian,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &bubbleRegionIndices,
							    const SIM_RawField &liquidSurface,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights,
							    const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
							    const UT_Array<UT_Vector3I> &bubbleBoundaryCells,
							    const UT_Array<bool> &isBubbleRemoved,
							    const SolveReal liquidDensity) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();
    
    UT_ThreadedAlgorithm buildBubbleLaplacianAlgorithm;
    buildBubbleLaplacianAlgorithm.run([&](const UT_JobInfo &info)
    {
	exint start, end;
	exint elements = bubbleBoundaryCells.size();

	UT_Array<SolveReal> &localBubbleLaplacian = parallelBubbleRegionLaplacian[info.job()];

	info.divideWork(elements, start, end);

	if (boss->opInterrupt())
	    return 0;
	
	for (exint i = start, iend = end; i < iend; ++i)
	{
	    UT_Vector3I cell = bubbleBoundaryCells[i];

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::LIQUID_CELL);

#if !defined(NDEBUG)
	    bool hasBubbleNeighbour = false;
#endif
	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

		    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
			continue;

		    UT_Vector3I face = cellToFaceMap(cell, axis, direction);

		    SolveReal weight = getFieldValue(*cutCellWeights[axis], face);

		    if (weight > 0)
		    {
			if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::BUBBLE_CELL)
			{
#if !defined(NDEBUG)
			    hasBubbleNeighbour = true;
#endif

			    exint bubbleRegion = getFieldValue(bubbleRegionIndices, adjacentCell);
			    assert(bubbleRegion >= 0);

			    if (isBubbleRemoved[bubbleRegion])
				continue;

			    // TODO: pre-compute weights
			    SolveReal phi0 = getFieldValue(liquidSurface, cell);
			    SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

			    SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
			    theta = SYSclamp(theta, .01, 1.);

			    SolveReal localDensity = theta * liquidDensity;

			    weight /= localDensity;

			    localBubbleLaplacian[bubbleRegion] -= weight * smootherDestinationGrid(cell);
			}
		    }
		}

	    assert(hasBubbleNeighbour);
	}

	return 0;
     });
}

void
HDK_ConstraintBubblePressureSolver::applyBubbleSmoother(UT_Array<SolveReal> &bubbleRegionDestinations,
							const UT_Array<SolveReal> &bubbleRegionSources,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &bubbleRegionIndices,
							const SIM_RawField &liquidSurface,
							const std::array<const SIM_RawField *, 3> &cutCellWeights,
							const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
							const UT_Array<SolveReal> &bubbleRegionDiagonals,
							const UT_Array<UT_Vector3I> &bubbleBoundaryCells,
							const UT_Array<bool> &isBubbleRemoved,
							const SolveReal liquidDensity) const
{
    int threadCount = UT_Thread::getNumProcessors();
    UT_Array<UT_Array<SolveReal>> parallelBubbleRegionalLaplacian;
    parallelBubbleRegionalLaplacian.setSize(threadCount);

    const exint bubbleRegionCount = bubbleRegionDestinations.size();

    for (int thread = 0; thread < threadCount; ++thread)
    {
	parallelBubbleRegionalLaplacian[thread].setSize(bubbleRegionCount);
	parallelBubbleRegionalLaplacian[thread].constant(0);
    }

    buildBubbleLaplacian(parallelBubbleRegionalLaplacian,
			    materialCellLabels,
			    bubbleRegionIndices,
			    liquidSurface,
			    cutCellWeights,
			    smootherDestinationGrid,
			    bubbleBoundaryCells,
			    isBubbleRemoved,
			    liquidDensity);

    UT_Array<SolveReal> bubbleRegionLaplacian;
    bubbleRegionLaplacian.setSize(bubbleRegionCount);
    bubbleRegionLaplacian.constant(0);

    // Collect Laplacian
    for (int thread = 0; thread < threadCount; ++thread)
	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (!isBubbleRemoved[bubbleRegion])
		bubbleRegionLaplacian[bubbleRegion] += parallelBubbleRegionalLaplacian[thread][bubbleRegion];
	    else assert(parallelBubbleRegionalLaplacian[thread][bubbleRegion] == 0);
	}

    constexpr SolveReal dampedWeight = 2./ 3.;

    // Update bubble destination list
    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
    {
	if (!isBubbleRemoved[bubbleRegion])
	{
	    bubbleRegionLaplacian[bubbleRegion] += bubbleRegionDiagonals[bubbleRegion] * bubbleRegionDestinations[bubbleRegion];
	    assert(bubbleRegionDiagonals[bubbleRegion] > 0);

	    SolveReal residual = bubbleRegionSources[bubbleRegion] - bubbleRegionLaplacian[bubbleRegion];
	    residual /= bubbleRegionDiagonals[bubbleRegion];

	    bubbleRegionDestinations[bubbleRegion] += dampedWeight * residual;
	}
    }
}

void
HDK_ConstraintBubblePressureSolver::applyLiquidSmoother(UT_Array<SolveReal> &tempSmootherDestinationValues,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &bubbleRegionIndices,
							const SIM_RawField &liquidSurface,
							const std::array<const SIM_RawField *, 3> &cutCellWeights,
							const UT_VoxelArray<int> &mgDomainCellLabels,
							const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
							const UT_VoxelArray<SolveReal> &smootherSourceGrid,
							const UT_Array<SolveReal> &bubbleRegionDestinations,
							const UT_Array<UT_Vector3I> &bubbleSmootherCells,
							const UT_Vector3I &mgExpandedOffset,
							const UT_Array<bool> &isBubbleRemoved,
							const SolveReal liquidDensity) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

    UT_Interrupt *boss = UTgetInterrupt();
    
    constexpr SolveReal dampedWeight = 2. / 3.;

    UTparallelFor(UT_BlockedRange<exint>(0, bubbleSmootherCells.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;
	
	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = bubbleSmootherCells[cellIndex];

	    SolveReal diagonal = 0;
	    SolveReal laplacian = 0;

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::LIQUID_CELL);

	    UT_Vector3I expandedCell = cell + mgExpandedOffset;
	    if (mgDomainCellLabels(expandedCell) == MGCellLabels::INTERIOR_CELL)
	    {
		for (int axis : {0,1,2})
		    for (int direction : {0,1})
		    {
			UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
			laplacian -= smootherDestinationGrid(adjacentCell);
		    }

		laplacian /= liquidDensity;
		diagonal = 6. / liquidDensity;
	    }
	    else
	    {
		assert(mgDomainCellLabels(expandedCell) == MGCellLabels::BOUNDARY_CELL);

		for (int axis : {0,1,2})
		    for (int direction : {0,1})
		    {
			UT_Vector3I face = cellToFaceMap(cell, axis, direction);
			SolveReal weight = getFieldValue(*cutCellWeights[axis], face);

			if (weight > 0)
			{
			    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

			    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
				continue;

			    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::LIQUID_CELL)
			    {
				weight /= liquidDensity;
				diagonal += weight;
				laplacian -= weight * smootherDestinationGrid(adjacentCell);
			    }
			    else
			    {
				assert(getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::BUBBLE_CELL);

				SolveReal phi0 = getFieldValue(liquidSurface, cell);
				SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

				SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
				theta = SYSclamp(theta, .01, 1.);

				SolveReal localDensity = theta * liquidDensity;

				exint bubbleRegion = getFieldValue(bubbleRegionIndices, adjacentCell);
				assert(bubbleRegion >= 0);

				weight /= localDensity;

				if (!isBubbleRemoved[bubbleRegion])
				    laplacian -= weight * bubbleRegionDestinations[bubbleRegion];

				diagonal += weight;
			    }			    
			}
		    }
	    }

	    laplacian += diagonal * smootherDestinationGrid(cell);
	    assert(diagonal > 0);
	    SolveReal residual = smootherSourceGrid(cell) - laplacian;
	    residual /= diagonal;

	    tempSmootherDestinationValues[cellIndex] = smootherDestinationGrid(cell) + dampedWeight * residual;
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::copyTempDestinationToGrid(UT_VoxelArray<SolveReal> &smootherDestinationGrid,
							const UT_Array<SolveReal> &tempSmootherDestinationValues,
							const UT_Array<UT_Vector3I> &bubbleSmootherCells) const
{
    UT_Interrupt *boss = UTgetInterrupt();
    
    UTparallelForLightItems(UT_BlockedRange<exint>(0, bubbleSmootherCells.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;
	
	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = bubbleSmootherCells[cellIndex];
	    smootherDestinationGrid.setValue(cell, tempSmootherDestinationValues[cellIndex]);
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::copySourceToMGGrid(UT_VoxelArray<SolveReal> &mgSourceGrid,
						const SIM_RawIndexField &materialCellLabels,
						const SIM_RawIndexField &liquidCellIndices,
						const UT_VoxelArray<int> &mgDomainCellLabels,
						const Vector &sourceVector,
						const UT_Vector3I &mgExpandedOffset,
						const SolveReal liquidDensity) const
{
    using SIM::FieldUtils::getFieldValue;
    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

    UT_Interrupt *boss = UTgetInterrupt();
    UTparallelForEachNumber(mgDomainCellLabels.numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIterator<int> vit;
	vit.setConstArray(&mgDomainCellLabels);
        
	if (boss->opInterrupt())
	    return;

	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (!vit.isTileConstant() ||
		vit.getValue() == MGCellLabels::INTERIOR_CELL ||
		vit.getValue() == MGCellLabels::BOUNDARY_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == MGCellLabels::INTERIOR_CELL ||
			vit.getValue() == MGCellLabels::BOUNDARY_CELL)
		    {
			UT_Vector3I expandedCell(vit.x(), vit.y(), vit.z());
			UT_Vector3I cell = expandedCell - mgExpandedOffset;

			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::LIQUID_CELL);

			exint liquidIndex = getFieldValue(liquidCellIndices, cell);
			assert(liquidIndex >= 0);

			mgSourceGrid.setValue(expandedCell, liquidDensity * sourceVector(liquidIndex));
		    }
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::applyDirichletToMGGrid(UT_VoxelArray<SolveReal> &mgSourceGrid,
							    UT_VoxelArray<SolveReal> &mgDestinationGrid,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &liquidCellIndices,
							    const SIM_RawIndexField &bubbleRegionIndices,
							    const SIM_RawField &liquidSurface,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights,
							    const UT_VoxelArray<int> &mgDomainCellLabels,
							    const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
							    const Vector &sourceVector,
							    const UT_Array<SolveReal> &bubbleRegionDestinations,
							    const UT_Array<UT_Vector3I> &bubbleSmootherCells,
							    const UT_Vector3I &mgExpandedOffset,
							    const UT_Array<bool> &isBubbleRemoved,
							    const SolveReal liquidDensity) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;
    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

    UT_Interrupt *boss = UTgetInterrupt();
    
    UTparallelFor(UT_BlockedRange<exint>(0, bubbleSmootherCells.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;
	
	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = bubbleSmootherCells[cellIndex];
	    UT_Vector3I expandedCell = cell + mgExpandedOffset;
	    
	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::LIQUID_CELL);

	    // Copy smoother values over
	    mgDestinationGrid.setValue(expandedCell, smootherDestinationGrid(cell));

	    // Apply Dirichlet to source grid
	    if (mgDomainCellLabels(expandedCell) == MGCellLabels::BOUNDARY_CELL)
	    {
		SolveReal laplacian = 0;

		for (int axis : {0,1,2})
		    for (int direction : {0,1})
		    {
			UT_Vector3I face = cellToFaceMap(cell, axis, direction);
			SolveReal weight = getFieldValue(*cutCellWeights[axis], face);

			if (weight > 0)
			{
			    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

			    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
				continue;

			    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::BUBBLE_CELL)
			    {
				SolveReal phi0 = getFieldValue(liquidSurface, cell);
				SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

				SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
				theta = SYSclamp(theta, .01, 1.);

				SolveReal localDensity = theta * liquidDensity;

				exint bubbleRegion = getFieldValue(bubbleRegionIndices, adjacentCell);
				assert(bubbleRegion >= 0);

				weight *= liquidDensity / localDensity;

				if (!isBubbleRemoved[bubbleRegion])
				    laplacian -= weight * bubbleRegionDestinations[bubbleRegion];
			    }
			}
		    }

		if (laplacian != 0)
		{
		    exint liquidIndex = getFieldValue(liquidCellIndices, cell);
		    assert(liquidIndex >= 0);
		    mgSourceGrid.setValue(expandedCell, liquidDensity * sourceVector(liquidIndex) - laplacian);
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::copyMGToSmoother(UT_VoxelArray<SolveReal> &smootherDestinationGrid,
						const SIM_RawIndexField &materialCellLabels,
						const UT_VoxelArray<int> &mgDomainCellLabels,
						const UT_VoxelArray<SolveReal> &mgDestinationGrid,
						const UT_Array<UT_Vector3I> &liquidCopyCells,
						const UT_Vector3I &mgExpandedOffset) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;
    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

    UT_Interrupt *boss = UTgetInterrupt();
    
    UTparallelFor(UT_BlockedRange<exint>(0, liquidCopyCells.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;
	
	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = liquidCopyCells[cellIndex];
	    UT_Vector3I expandedCell = cell + mgExpandedOffset;

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::LIQUID_CELL);
	    assert(mgDomainCellLabels(expandedCell) == MGCellLabels::INTERIOR_CELL ||
		    mgDomainCellLabels(expandedCell) == MGCellLabels::BOUNDARY_CELL);

	    smootherDestinationGrid.setValue(cell, mgDestinationGrid(expandedCell));
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::copyMGToDestinationVector(Vector &destinationVector,
    							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &liquidCellIndices,
							const UT_VoxelArray<int> &mgDomainCellLabels,
							const UT_VoxelArray<SolveReal> &mgDestinationGrid,
							const UT_Vector3I &mgExpandedOffset) const
{
    using SIM::FieldUtils::getFieldValue;
    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

    UT_Interrupt *boss = UTgetInterrupt();
    UTparallelForEachNumber(materialCellLabels.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
        
	if (boss->opInterrupt())
	    return;

	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::LIQUID_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == MaterialLabels::LIQUID_CELL)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());
			UT_Vector3I expandedCell = cell + mgExpandedOffset;

			exint liquidIndex = getFieldValue(liquidCellIndices, cell);
			assert(liquidIndex >= 0);

			assert(mgDomainCellLabels(expandedCell) == MGCellLabels::BOUNDARY_CELL ||
				mgDomainCellLabels(expandedCell) == MGCellLabels::INTERIOR_CELL);

			destinationVector(liquidIndex) = mgDestinationGrid(expandedCell);
		    }
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::copySmootherToDestinationVector(Vector &destinationVector,
							    const SIM_RawIndexField &liquidCellIndices,
							    const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
							    const UT_Array<UT_Vector3I> &bubbleSmootherCells) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();
    
    UTparallelFor(UT_BlockedRange<exint>(0, bubbleSmootherCells.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;
	
	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = bubbleSmootherCells[cellIndex];

	    exint liquidIndex = getFieldValue(liquidCellIndices, cell);
	    assert(liquidIndex >= 0);

	    destinationVector(liquidIndex) = smootherDestinationGrid(cell);
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::applySolutionToPressure(SIM_RawField &pressure,
						    const SIM_RawIndexField &indices,
						    const Vector &solution,
						    const exint offset) const
{
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(indices.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(indices.field());

	if (boss->opInterrupt())
	    return;

	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (vit.isTileConstant())
	    {
		exint row = vit.getValue();
		if (row >= 0)
		{
		    // TODO: remove once free surface simulation is working
		    assert(offset > 0); 
		    int tileNumber = vit.getLinearTileNum();

		    UT_Vector3i tileIndex;
		    indices.field()->linearTileToXYZ(tileNumber, tileIndex[0], tileIndex[1], tileIndex[2]);

		    pressure.fieldNC()->getTile(tileIndex[0], tileIndex[1], tileIndex[2])->makeConstant(solution(row + offset));
		}
	    }
	    else
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    exint row = vit.getValue();

		    if (row >= 0)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());
			setFieldValue(pressure, cell, solution(row + offset));
		    }
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::applyPressureGradient(SIM_RawField &velocity,
							    const SIM_RawField &validFaces,
							    const SIM_RawField &pressure,
							    const SIM_RawField &liquidSurface,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &bubbleRegionIndices,
							    const SolveReal &liquidDensity,
							    const int axis) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;
    using SIM::FieldUtils::faceToCellMap;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(validFaces.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorF vit;
	vit.setConstArray(validFaces.field());

    	if (boss->opInterrupt())
	    return;

	UT_Vector3I faceRes;
	faceRes[0] = validFaces.getXRes();
	faceRes[1] = validFaces.getYRes();
	faceRes[2] = validFaces.getZRes();
	
	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (!vit.isTileConstant() || vit.getValue() == HDK::Utilities::VALID_FACE)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == HDK::Utilities::VALID_FACE)
		    {
			UT_Vector3I face(vit.x(), vit.y(), vit.z());

			UT_Vector3I backwardCell = faceToCellMap(face, axis, 0);
			UT_Vector3I forwardCell = faceToCellMap(face, axis, 1);

			assert(backwardCell[axis] >= 0 && forwardCell[axis] < materialCellLabels.getVoxelRes()[axis]);

			auto backwardMaterial = getFieldValue(materialCellLabels, backwardCell);
			auto forwardMaterial = getFieldValue(materialCellLabels, forwardCell);

			assert(backwardMaterial == MaterialLabels::LIQUID_CELL || forwardMaterial == MaterialLabels::LIQUID_CELL);

			SolveReal localDensity;
			if (backwardMaterial == MaterialLabels::BUBBLE_CELL || forwardMaterial == MaterialLabels::BUBBLE_CELL)
			{
			    SolveReal backwardPhi = getFieldValue(liquidSurface, backwardCell);
			    SolveReal forwardPhi = getFieldValue(liquidSurface, forwardCell);

			    SolveReal theta = HDK::Utilities::computeGhostFluidWeight(backwardPhi, forwardPhi);
			    theta = SYSclamp(theta, .01, 1.);

			    localDensity = theta * liquidDensity;
			}
			else
			    localDensity = liquidDensity;

			// Pressure includes both liquid and constraint bubble pressure
			SolveReal gradient = getFieldValue(pressure, forwardCell) - getFieldValue(pressure, backwardCell);

			setFieldValue(velocity, face, getFieldValue(velocity, face) - gradient / localDensity);
		    }
		}
	    }
	}
    });
}

void
HDK_ConstraintBubblePressureSolver::computeResultingDivergence(UT_Array<SolveReal> &parallelAccumulatedLiquidDivergence,
							UT_Array<SolveReal> &parallelMaxLiquidDivergence,
							UT_Array<SolveReal> &parallelLiquidCellCount,
							UT_Array<UT_Array<SolveReal>> &parallelBubbleRegionDivergence,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &bubbleRegionIndices,
							const SIM_VectorField &velocity,
							const std::array<const SIM_RawField *, 3> &cutCellWeights,
							const SIM_VectorField *solidVelocity,
							const UT_Array<bool> &isBubbleTracked,
							const UT_Array<bool> &isBubbleRemoved) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm computeResultingDivergenceAlgorithm;
    computeResultingDivergenceAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	SolveReal &localAccumulatedLiquidDivergence = parallelAccumulatedLiquidDivergence[info.job()];
	SolveReal &localMaxLiquidDivergence = parallelMaxLiquidDivergence[info.job()];
	SolveReal &localLiquidCellCount = parallelLiquidCellCount[info.job()];

	UT_Array<SolveReal> &localBubbleRegionDivergence = parallelBubbleRegionDivergence[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::LIQUID_CELL ||
		vit.getValue() == MaterialLabels::BUBBLE_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::LIQUID_CELL)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			SolveReal divergence = 0;

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I face = cellToFaceMap(cell, axis, direction);
				SolveReal weight = getFieldValue(*cutCellWeights[axis], face);

				SolveReal sign = (direction == 0) ? -1 : 1;

				if (weight > 0)
				    divergence += sign * weight * getFieldValue(*velocity.getField(axis), face);
				if (solidVelocity != nullptr && weight < 1)
				{
				    UT_Vector3 point;
				    velocity.getField(axis)->indexToPos(face[0], face[1], face[2], point);

				    divergence += sign * (1. - weight) * solidVelocity->getField(axis)->getValue(point);
				}
			    }

			localAccumulatedLiquidDivergence += divergence;
			localMaxLiquidDivergence = std::max(localMaxLiquidDivergence, divergence);
			++localLiquidCellCount;
		    }
		    else if (vitt.getValue() == MaterialLabels::BUBBLE_CELL)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			exint bubbleRegion = getFieldValue(bubbleRegionIndices, cell);
			assert(bubbleRegion >= 0);		

			if (!isBubbleTracked[bubbleRegion] || isBubbleRemoved[bubbleRegion])
			    continue;

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				// We're assuming a solid static boundary so there is no divergence contribution.
				if (adjacentCell[axis] < 0 || adjacentCell[axis] >= bubbleRegionIndices.getVoxelRes()[axis])
				    continue;

				SolveReal sign = (direction == 0) ? 1. : -1.;
				UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				SolveReal weight = getFieldValue(*cutCellWeights[axis], face);

				if (weight > 0. && getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::LIQUID_CELL)
				    localBubbleRegionDivergence[bubbleRegion] += sign * weight * getFieldValue(*velocity.getField(axis), face);

				if (solidVelocity != nullptr && weight < 1.)
				{
				    UT_Vector3 point;
				    velocity.getField(axis)->indexToPos(face[0], face[1], face[2], point);

				    localBubbleRegionDivergence[bubbleRegion] += sign * (1. - weight) * solidVelocity->getField(axis)->getValue(point);
				}
			    }
		    }
		}
	    }
	}

	return 0;
    });
}