#include "HDK_TwoPhasePressureSolver.h"

#include <Eigen/LU>

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
#include <UT/UT_PerfMonAutoEvent.h>

#include <UT/UT_StopWatch.h>

void initializeSIM(void *)
{
    IMPLEMENT_DATAFACTORY(HDK_TwoPhasePressureSolver);
}

// Standard constructor, note that BaseClass was crated by the
// DECLARE_DATAFACTORY and provides an easy way to chain through
// the class hierarchy.
HDK_TwoPhasePressureSolver::HDK_TwoPhasePressureSolver(const SIM_DataFactory *factory)
    : BaseClass(factory)
{
}

HDK_TwoPhasePressureSolver::~HDK_TwoPhasePressureSolver()
{
}

const SIM_DopDescription* HDK_TwoPhasePressureSolver::getDopDescription()
{
    static PRM_Name	theSurfaceName(GAS_NAME_SURFACE, "Surface Field");
    static PRM_Default 	theSurfaceDefault(0,"surface");

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

    static PRM_Name    theUseOldPressureName("useOldPressure", "Use Old Pressure as an Initial Guess");

    static PRM_Name 	theGoalDivergenceName("goalDivergence", "Goal Divergence");
    static PRM_Default 	theGoalDivergenceNameDefault(0, "goaldivergence");

    static PRM_Name 	theDensityName(GAS_NAME_DENSITY, "Liquid Density Field");
    static PRM_Default 	theDensityDefault(0, "massdensity");

    static PRM_Name     theBubbleDensityName("bubbleDensity", "Bubble Density");
    static PRM_Default  theBubbleDensityDefault(1000);

    static PRM_Name 	theValidFacesName("validFaces", "Valid Faces Field");

    static PRM_Name 	theToleranceName(SIM_NAME_TOLERANCE, "Error Tolerance");
    static PRM_Default 	theToleranceDefault(1e-5);

    static PRM_Name 	theMaxIterations("maxIterations", "Max Solver Iterations");
    static PRM_Default 	theMaxIterationsDefault(2500);

    static PRM_Name	theTrackBubbleIDsName("trackBubbleIDs", "Track Bubble IDs");

    static PRM_Name	theTrackBubbleVolumesName("trackBubbleVolumes", "Track Bubble Volumes");
    static PRM_Conditional theTrackBubbleVolumesDisable("{ trackBubbleIDs == 0 }");
    
    static PRM_Name	theCorrectVolumeDriftName("correctVolumeDrift", "Correct Volume Drift");
    static PRM_Conditional theCorrectVolumeDriftDisable("{ trackBubbleVolumes == 0 }");

    static PRM_Name	theGeometryName(GAS_NAME_GEOMETRY, "Bubble ID Tracking Geometry");
    static PRM_Default  theGeometryNameDefault(0, "Geometry");

    static PRM_Name	thePrintCellsName("printCells", "Print Cell Activity");
    static PRM_Name	theOnlyPrintCellsName("onlyPrintCells", "Only Print Cells");

    static PRM_Name    thePrintedCellsGeometryName("cellActivityGeometry", "Cell Activity Geometry");
    static PRM_Default thePrintedCellsGeometryDefault(0, "CellActivityGeometry");

    static PRM_Name theBubbleVolumeGeometryName("bubbleVolumeGeometry", "Bubble Volume Tracking Geometry");
    static PRM_Default  theBubbleVolumeGeometryNameDefault(0, "BubbleVolume");

    static PRM_Template	theTemplates[] =
    {
        PRM_Template(PRM_STRING, 1, &theSurfaceName, &theSurfaceDefault),
    	PRM_Template(PRM_STRING, 1, &theVelocityName, &theVelocityDefault),

    	PRM_Template(PRM_STRING, 1, &theSolidSurfaceName, &theSolidSurfaceDefault),

    	PRM_Template(PRM_STRING, 1, &theSolidVelocityName, &theSolidVelocityDefault),
    	PRM_Template(PRM_STRING, 1, &theCutCellWeightsName, &theCutCellWeightsDefault),
    
	PRM_Template(PRM_STRING, 1, &thePressureName, &thePressureNameDefault),

	PRM_Template(PRM_TOGGLE, 1, &theUseOldPressureName, PRMoneDefaults),

	PRM_Template(PRM_STRING, 1, &theGoalDivergenceName, &theGoalDivergenceNameDefault),

	PRM_Template(PRM_STRING, 1, &theDensityName, &theDensityDefault),

	PRM_Template(PRM_FLT, 1, &theBubbleDensityName, &theBubbleDensityDefault),
	
        PRM_Template(PRM_STRING, 1, &theValidFacesName),

    	PRM_Template(PRM_FLT, 1, &theToleranceName, &theToleranceDefault),
    	PRM_Template(PRM_FLT, 1, &theMaxIterations, &theMaxIterationsDefault),

        PRM_Template(PRM_TOGGLE, 1, &theTrackBubbleIDsName, PRMzeroDefaults),

    	PRM_Template(PRM_TOGGLE, 1, &theTrackBubbleVolumesName, PRMzeroDefaults,
                        0, 0, 0, 0, 1, 0, &theTrackBubbleVolumesDisable),

    	PRM_Template(PRM_TOGGLE, 1, &theCorrectVolumeDriftName, PRMzeroDefaults,
                        0, 0, 0, 0, 1, 0, &theCorrectVolumeDriftDisable),

        PRM_Template(PRM_STRING, 1, &theGeometryName, &theGeometryNameDefault),

    	PRM_Template(PRM_STRING, 1, &theBubbleVolumeGeometryName, &theBubbleVolumeGeometryNameDefault),

	PRM_Template(PRM_TOGGLE, 1, &thePrintCellsName, PRMzeroDefaults),
	PRM_Template(PRM_TOGGLE, 1, &theOnlyPrintCellsName, PRMzeroDefaults),

	PRM_Template(PRM_STRING, 1, &thePrintedCellsGeometryName,
				    &thePrintedCellsGeometryDefault),

    	PRM_Template()
    };

    static SIM_DopDescription theDopDescription(true,
						"HDK_TwoPhasePressureSolver",
						"HDK Two Phase Pressure Solver",
						"$OS",
						classname(),
						theTemplates);

    setGasDescription(theDopDescription);

    return &theDopDescription;
}

bool
HDK_TwoPhasePressureSolver::solveGasSubclass(SIM_Engine &engine,
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

    // Load goal divergence
    const SIM_ScalarField *goalDivergence = getConstScalarField(obj, "goalDivergence");
    
    if (goalDivergence != nullptr && !goalDivergence->getField()->isAligned(&liquidSurface))
    {
	addError(obj, SIM_MESSAGE, "Goal divergence field sampling needs to match the surface field", UT_ERROR_ABORT);
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
    const SolveReal bubbleDensity = getBubbleDensity();

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
    // Build unique indices for each active fluid cell
    //
    ////////////////////////////////////////////

    exint activeCellCount = 0;
    SIM_RawIndexField activeCellIndices;
    {
	UT_PerfMonAutoSolveEvent event(this, "Build active cell indices");
    	std::cout << "\n// Build active cell indices" << std::endl;

	activeCellCount = buildActiveCellIndices(activeCellIndices, materialCellLabels);
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
    // the entire volume.
    //
    ////////////////////////////////////////////

    SIM_RawIndexField activeRegionIndices;
    exint activeRegionCount;

    {
    	std::cout << "\n// Building active connected components" << std::endl;
    	UT_PerfMonAutoSolveEvent event(this, "Build active connected components");

	activeRegionCount = buildActiveRegionIndices(activeRegionIndices, materialCellLabels, cutCellWeights);

	std::cout << "  Active region count: " << activeRegionCount << std::endl;
    }

    ////////////////////////////////////////////
    //
    //
    ////////////////////////////////////////////

    bool doTrackBubbleVolumes = getTrackBubbleVolumes();
    bool doTrackBubbleIDs = doTrackBubbleVolumes ? true : getTrackBubbleIDs();

    UT_Array<bool> areBubblesTracked;
    UT_Array<SolveReal> newBubbleRestVolumes;

    const int threadCount = UT_Thread::getNumProcessors();

    if (doTrackBubbleIDs)
    {
	UT_Array<UT_Array<bool>> oldToNewBubbleMap;

	buildBubbleIDs(areBubblesTracked,
			oldToNewBubbleMap,
			materialCellLabels,
			obj,
			liquidSurface,
			bubbleRegionIndices,
			cutCellWeights,
			bubbleRegionCount,
			doTrackBubbleVolumes);

	if (doTrackBubbleVolumes)
	{
	    buildBubbleRestVolumes(newBubbleRestVolumes,
				    obj,
				    oldToNewBubbleMap,
				    areBubblesTracked,
				    bubbleRegionSizes,
				    bubbleRegionCount);
	}
    }

    ////////////////////////////////////////////
    //
    // Print geometry for display and debugging
    //
    ////////////////////////////////////////////

    if (getPrintCells())
    {
	std::cout << "// Printing cell geometry" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Printing cell geometry");

	printCellGeometry(obj,
			    activeCellIndices,
			    materialCellLabels,
			    velocity->getOrig(),
			    dx);

        if (getOnlyPrintCells())
	    return true;
    }

    ////////////////////////////////////////////
    //
    // Build rhs for active cells
    //
    ////////////////////////////////////////////

    Vector rhsVector = Vector::Zero(activeCellCount);

    {
	std::cout << "\n// Build RHS" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Build RHS");

	buildRHS(rhsVector,
		    materialCellLabels,
		    activeCellIndices,
		    *velocity,
		    cutCellWeights,
		    solidVelocity,
		    dx);
    }

    ////////////////////////////////////////////
    //
    // Apply goal divergence for liquid cells
    //
    ////////////////////////////////////////////
    
    if (goalDivergence != nullptr)
    {
	std::cout << "\n// Apply goal divergence" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Apply goal divergence");

	applyGoalDivergence(rhsVector,
			    materialCellLabels,
			    activeCellIndices,
			    *goalDivergence->getField());
    }

    ////////////////////////////////////////////
    //
    // Remove the null space introduced by latent
    // divergence in the system
    //
    ////////////////////////////////////////////

    {
	std::cout << "\n// Project background null space" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Project background null space");

	removeBackgroundNullSpace(rhsVector, activeCellIndices, activeRegionIndices, activeRegionCount);
    }

    if (doTrackBubbleIDs)
    {
	std::cout << "\n// Apply divergence for untracked bubbles and volume correction" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Apply divergence for untracked bubbles and volume correction");

	applyBubbleCorrection(rhsVector,
				materialCellLabels,
				activeCellIndices,
				activeRegionIndices,
				bubbleRegionIndices,
				newBubbleRestVolumes,
				bubbleRegionSizes,
				areBubblesTracked,
				dx, dt,
				doTrackBubbleVolumes,
				getCorrectVolumeDrift(),
				activeRegionCount,
				bubbleRegionCount);
    }

    ////////////////////////////////////////////
    //
    // Build poisson sparse elements for active cells
    //
    ////////////////////////////////////////////

    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> poissonMatrix(activeCellCount, activeCellCount);
    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> diagonalPrecondMatrix(activeCellCount, activeCellCount);

    {
	std::cout << "\n// Build Poisson system" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Build Poisson system");

	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelPoissonElements(threadCount);
	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelDiagonalPrecondElements(threadCount);

	buildLinearSystem(parallelPoissonElements,
			    parallelDiagonalPrecondElements,
			    materialCellLabels,
			    activeCellIndices,
			    cutCellWeights,
			    liquidSurface,
			    liquidDensity,
			    bubbleDensity);

	// Compile poisson system elements
	{
	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelPoissonElements[thread].size();

	    std::vector<Eigen::Triplet<SolveReal>> poissonElements;
	    poissonElements.reserve(listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
		poissonElements.insert(poissonElements.end(), parallelPoissonElements[thread].begin(), parallelPoissonElements[thread].end());

	    poissonMatrix.setFromTriplets(poissonElements.begin(), poissonElements.end());
	    poissonMatrix.makeCompressed();
	}

	// Compile diagonal preconditioner elements
	{
	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelDiagonalPrecondElements[thread].size();

	    std::vector<Eigen::Triplet<SolveReal>> diagonalPrecondElements;
	    diagonalPrecondElements.reserve(listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
		diagonalPrecondElements.insert(diagonalPrecondElements.end(), parallelDiagonalPrecondElements[thread].begin(), parallelDiagonalPrecondElements[thread].end());

	    diagonalPrecondMatrix.setFromTriplets(diagonalPrecondElements.begin(), diagonalPrecondElements.end());
	    diagonalPrecondMatrix.makeCompressed();
	}
    }

    ////////////////////////////////////////////
    //
    // Set valid faces
    //
    ////////////////////////////////////////////

    {
	UT_PerfMonAutoSolveEvent event(this, "Build valid face flags");
	std::cout << "// Build valid face flags" << std::endl;

	buildValidFaces(*validFaces,
			materialCellLabels,
			cutCellWeights);
    }

    ////////////////////////////////////////////
    //
    // Set initial guess for liquid pressure
    //
    ////////////////////////////////////////////

    Vector solutionVector = Vector::Zero(activeCellCount);

    if (getUseOldPressure())
    {
	std::cout << "\n// Apply warm start pressure" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Apply warm start pressure");

	HDK::Utilities::applyInitialGuess(solutionVector,
					    *pressure,
					    activeCellIndices);
    }

    ////////////////////////////////////////////
    //
    // Solve linear system
    //
    ////////////////////////////////////////////

    {
	std::cout << "// Solve affine tiled linear system" << std::endl;
    	UT_PerfMonAutoSolveEvent event(this, "Solve affine tiled linear system");

	auto MatrixVectorMultiply = [&](Vector &destinationVector, const Vector &sourceVector)
	{
	    destinationVector.noalias() = poissonMatrix * sourceVector;
	};

	auto DiagonalPreconditioner = [&](Vector &destinationVector, const Vector &sourceVector)
	{
	    destinationVector.noalias() = diagonalPrecondMatrix * sourceVector;
	};

	HDK::Utilities::solveConjugateGradient(solutionVector,
						rhsVector,
						MatrixVectorMultiply,
						DiagonalPreconditioner,
						getTolerance(),
						getMaxIterations());

	Vector residualVector(rhsVector.rows());
	MatrixVectorMultiply(residualVector, solutionVector);
	residualVector = rhsVector - residualVector;

	std::cout << "Re-computed residual: " << std::sqrt(residualVector.squaredNorm() / rhsVector.squaredNorm()) << std::endl;
	std::cout << "L-infinity residual: " << residualVector.lpNorm<Eigen::Infinity>() << std::endl;
    }

    ////////////////////////////////////////////
    //
    // Copy pressure to grid
    //
    ////////////////////////////////////////////

    // TODO: consider extrapolating pressures to get better initial guesses

    {
	UT_PerfMonAutoSolveEvent event(this, "Copy pressure solution to grid");
	std::cout << "// Copy pressure solution to grid" << std::endl;

	pressure->makeConstant(0);

	copyPressureVectorToGrid(*pressure,
				    activeCellIndices,
				    solutionVector);
    }

    ////////////////////////////////////////////
    //
    // Apply solution to velocity field
    //
    ////////////////////////////////////////////

    {
	UT_PerfMonAutoSolveEvent event(this, "Update velocities with pressure and affine fields");
	std::cout << "// Update velocities with pressure and affine fields" << std::endl;

	for (int axis : {0,1,2})
	{
	    applySolutionToVelocity(*velocity->getField(axis),
				    *validFaces->getField(axis),
				    materialCellLabels,
				    activeCellIndices,
				    *pressure,
				    liquidSurface,
				    *cutCellWeights[axis],
				    liquidDensity,
				    bubbleDensity,
				    axis,
				    dx);
	}
    }

    ////////////////////////////////////////////
    //
    // Debug check to verify that 
    // 1. each interior region is divergence free
    // 2. each active cell is divergence free
    //
    ////////////////////////////////////////////

    {
	UT_PerfMonAutoSolveEvent event(this, "Verify divergence-free constraint");
	std::cout << "// Verify divergence-free constraint" << std::endl;

	std::array<SolveReal, 3> divergenceStats;
	computeResultingDivergence(divergenceStats,
				    materialCellLabels,
				    *velocity,
				    cutCellWeights,
				    solidVelocity);

	std::cout << "    Accumulated divergence: " << divergenceStats[0] << std::endl;
	std::cout << "    Average divergence: " << divergenceStats[0] / divergenceStats[2] << std::endl;
	std::cout << "    Max divergence: " << divergenceStats[1] << std::endl;
    }

    pressureField->pubHandleModification();
    velocity->pubHandleModification();
    validFaces->pubHandleModification();

    return true;
}

void
HDK_TwoPhasePressureSolver::setMaterialCellLabels(SIM_RawIndexField &materialCellLabels) const
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

exint
HDK_TwoPhasePressureSolver::buildActiveCellIndices(SIM_RawIndexField &activeCellIndices,
						    const SIM_RawIndexField &materialCellLabels) const
{
    using SIM::FieldUtils::setFieldValue;

    activeCellIndices.match(materialCellLabels);
    activeCellIndices.makeConstant(HDK::Utilities::UNLABELLED_CELL);

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(materialCellLabels.field());

    UT_VoxelTileIteratorI vitt;

    exint activeCellCount = 0;

    // Build liquid cell indices
    UT_Interrupt *boss = UTgetInterrupt();
    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
	if (boss->opInterrupt())
	    break;

	if (!vit.isTileConstant() ||
	    vit.getValue() == MaterialLabels::LIQUID_CELL ||
	    vit.getValue() == MaterialLabels::BUBBLE_CELL)
	{
	    vitt.setTile(vit);

	    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
	    {
		if (vitt.getValue() == MaterialLabels::LIQUID_CELL ||
		    vitt.getValue() == MaterialLabels::BUBBLE_CELL)
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
		    setFieldValue(activeCellIndices, cell, activeCellCount++);
		}
	    }
	}
    }

    return activeCellCount;
}

void
HDK_TwoPhasePressureSolver::computeBubbleRegionSize(UT_Array<UT_Array<exint>> &parallelBubbleRegionSizes,
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
HDK_TwoPhasePressureSolver::buildBubbleRegionIndices(SIM_RawIndexField &bubbleRegionIndices,
							UT_Array<exint> &bubbleRegionSizes,
							const SIM_RawIndexField &materialCellLabels,
							std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    std::cout << "  Build connected bubble regions" << std::endl;

    SIM_VolumetricConnectedComponentBuilder<> interiorRegionBuilder(bubbleRegionIndices, materialCellLabels, cutCellWeights.data());
    
    exint bubbleRegionCount = interiorRegionBuilder.buildConnectedComponents([](const exint label)
    {
	return label == MaterialLabels::BUBBLE_CELL;
    });

    HDK::Utilities::overwriteIndices(bubbleRegionIndices,
					SIM_VolumetricConnectedComponentBuilder<>::INACTIVE_REGION,
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
HDK_TwoPhasePressureSolver::buildActiveRegionIndices(SIM_RawIndexField &activeRegionIndices,
							const SIM_RawIndexField &materialCellLabels,
							std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    activeRegionIndices.match(materialCellLabels);

    SIM_VolumetricConnectedComponentBuilder<> activeRegionRegionBuilder(activeRegionIndices, materialCellLabels, cutCellWeights.data());

    exint activeRegionCount = activeRegionRegionBuilder.buildConnectedComponents([](const exint label)
    {
	return label == MaterialLabels::BUBBLE_CELL || label == MaterialLabels::LIQUID_CELL;
    });

    HDK::Utilities::overwriteIndices(activeRegionIndices,
					SIM_VolumetricConnectedComponentBuilder<>::INACTIVE_REGION,
					HDK::Utilities::UNLABELLED_CELL);

    return activeRegionCount;
}

void
HDK_TwoPhasePressureSolver::trackBubbleIDs(UT_Array<bool> &areBubblesTracked,
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
    	UT_Array<bool> localAreBubblesTracked;
	localAreBubblesTracked.setSize(areBubblesTracked.size());
	localAreBubblesTracked.constant(false);

    	UT_Array<UT_Array<bool>> localOldToNewBubbleMap;

   	if (doMapBubbles)
    	{
    	    localOldToNewBubbleMap.setSize(oldToNewBubbleMap.size());

    	    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldToNewBubbleMap.size(); ++oldBubbleRegion)
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
    		    localAreBubblesTracked[newBubbleRegion] = true;

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
			    localAreBubblesTracked[newBubbleRegion] = true;
			    
			    if (doMapBubbles && bubbleID >= 0)
				localOldToNewBubbleMap[bubbleID][newBubbleRegion] = true;
			}
    		    }
    	    }
    	}

    	for (exint newBubbleRegion = 0; newBubbleRegion < areBubblesTracked.size(); ++newBubbleRegion)
    	{
    	    if (localAreBubblesTracked[newBubbleRegion])
                areBubblesTracked[newBubbleRegion] = true;
    	}

	if (doMapBubbles)
	{
	    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldToNewBubbleMap.size(); ++oldBubbleRegion)
	    {
		assert(oldToNewBubbleMap[oldBubbleRegion].size() == areBubblesTracked.size());
		for (exint newBubbleRegion = 0; newBubbleRegion < oldToNewBubbleMap[oldBubbleRegion].size(); ++newBubbleRegion)
		{
		    if (localOldToNewBubbleMap[oldBubbleRegion][newBubbleRegion])
			oldToNewBubbleMap[oldBubbleRegion][newBubbleRegion] = true;
		}
	    }
	}
    });
}

void
HDK_TwoPhasePressureSolver::updateBubbleIDs(GA_RWHandleI &bubbleIDHandle,
						const UT_Array<bool> &areBubblesTracked,
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
		if (bubbleRegion >= 0 && areBubblesTracked[bubbleRegion])
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

			    if (bubbleRegion >= 0 && areBubblesTracked[bubbleRegion])
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
HDK_TwoPhasePressureSolver::buildBubbleIDs(UT_Array<bool> &areBubblesTracked,
    					    UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
					    SIM_RawIndexField &materialCellLabels,
					    SIM_Object *obj,
					    const SIM_RawField &liquidSurface,
					    const SIM_RawIndexField &bubbleRegionIndices,
					    std::array<const SIM_RawField *, 3> &cutCellWeights,
					    const exint bubbleRegionCount,
					    const bool doTrackBubbleVolumes)
{
    UT_PerfMonAutoSolveEvent event(this, "Track bubble IDs");
    std::cout << "  Track bubble IDs" << std::endl;

    areBubblesTracked.setSize(bubbleRegionCount);
    areBubblesTracked.constant(false);

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
    SIM_GeometryCopy *bubbleVolumeGeometry = nullptr;
    GU_Detail *bubbleVolumeDetail = nullptr;

    UT_UniquePtr<SIM_GeometryAutoWriteLock> autoLockBubbleVolume;

    exint oldBubbleRegionCount = 0;

    if (doTrackBubbleVolumes)
    {
	bubbleVolumeGeometry = getOrCreateGeometry(obj, "bubbleVolumeGeometry");
	autoLockBubbleVolume = UTmakeUnique<SIM_GeometryAutoWriteLock>(bubbleVolumeGeometry, SIM_DATA_ID_PRESERVE);
	bubbleVolumeDetail = &autoLockBubbleVolume->getGdp();

	assert(bubbleVolumeDetail != nullptr);

	oldBubbleRegionCount = bubbleVolumeDetail->getNumPoints();
	std::cout << "  Debug check. Tracked bubble count: " << oldBubbleRegionCount << std::endl;
    }

    oldToNewBubbleMap.setSize(oldBubbleRegionCount);

    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldBubbleRegionCount; ++oldBubbleRegion)
    {
	oldToNewBubbleMap[oldBubbleRegion].setSize(bubbleRegionCount);
	oldToNewBubbleMap[oldBubbleRegion].constant(false);
    }

    // Track bubble IDs by filling in array. If a particle with a valid bubble ID is adjacent to a
    // newly detected bubble, indicate that the new bubble is indeed tracked.
    trackBubbleIDs(areBubblesTracked,
		    oldToNewBubbleMap,
		    bubbleIDHandle,
		    *particles,
		    bubbleRegionIndices,
		    liquidSurface);

    // Update bubble ID on particles for particles adjacent to bubbles that were tracked
    // from the previous timestep.
    updateBubbleIDs(bubbleIDHandle, areBubblesTracked, *particles, bubbleRegionIndices, liquidSurface);
    bubbleIDHandle.bumpDataId();
}

void
HDK_TwoPhasePressureSolver::distributeOldBubbleVolumes(UT_Array<SolveReal> &newBubbleRestVolumes,
							UT_Array<bool> &areBubblesUninitialized,
							const GA_RWHandleF &bubbleVolumeHandle,
							const UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
							const UT_Array<bool> &areBubblesTracked,
							const UT_Array<exint> &newBubbleRegionSizes) const
{
    constexpr exint UNVISITED_BUBBLE = -2;
    constexpr exint VISITED_BUBBLE = -1;

    std::cout << "    Distribute old bubble rest volumes" << std::endl;

    // In order to update the new bubble volumes based on the old bubble volumes,
    // we need to build a mapping between the two. We can handle this with a bi-partite graph.
    
    exint oldBubbleRegionCount = oldToNewBubbleMap.size();
    exint newBubbleRegionCount = newBubbleRegionSizes.size();

    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldBubbleRegionCount; ++oldBubbleRegion)
	assert(newBubbleRegionCount == oldToNewBubbleMap[oldBubbleRegion].size());

    UT_Array<exint> oldBubbleComponents;
    oldBubbleComponents.setSize(oldBubbleRegionCount);
    oldBubbleComponents.constant(UNVISITED_BUBBLE);

    UT_Array<exint> newBubbleComponents;
    newBubbleComponents.setSize(newBubbleRegionCount);
    newBubbleComponents.constant(UNVISITED_BUBBLE);

    UT_Array<exint> oldBubblesToVisit;

    exint connectedComponentCount = 0;

    // Build connected components
    for (exint oldBubbleRegion = 0; oldBubbleRegion < oldBubbleRegionCount; ++oldBubbleRegion)
    {
	// Start at an unvisited old bubble region to begin building the connected components in the bipartite graph.
	if (oldBubbleComponents[oldBubbleRegion] == UNVISITED_BUBBLE)
	{
	    oldBubblesToVisit.clear();
	    oldBubblesToVisit.append(oldBubbleRegion);

	    oldBubbleComponents[oldBubbleRegion] = VISITED_BUBBLE;

	    int oldToNewEdgeCount = 0;
	    while(!oldBubblesToVisit.isEmpty())
	    {
		exint localOldBubbleRegion = oldBubblesToVisit.last();
		oldBubblesToVisit.removeLast();

		assert(oldBubbleComponents[localOldBubbleRegion] == VISITED_BUBBLE);

		oldBubbleComponents[localOldBubbleRegion] = connectedComponentCount;

		// Sweep through the set of new bubble regions and check for an edge from the old bubble to
		// the current new bubble. If an edge exists, continue with search.
		for (exint newBubbleRegion = 0; newBubbleRegion < newBubbleRegionCount; ++newBubbleRegion)
		{
		    if (oldToNewBubbleMap[localOldBubbleRegion][newBubbleRegion])
		    {
			assert(areBubblesTracked[newBubbleRegion]);

			if (newBubbleComponents[newBubbleRegion] == UNVISITED_BUBBLE)
			{
			    // Loop over the set of old bubble and check for an edge from this new bubble back to an old
			    // bubble. If that old bubble has not yet been visited, add it to the queue.
			    for (exint adjacentOldBubbleRegion = 0; adjacentOldBubbleRegion < oldBubbleRegionCount; ++adjacentOldBubbleRegion)
			    {
				if (oldToNewBubbleMap[adjacentOldBubbleRegion][newBubbleRegion])
				{
				    if (oldBubbleComponents[adjacentOldBubbleRegion] == UNVISITED_BUBBLE)
				    {
					oldBubblesToVisit.append(adjacentOldBubbleRegion);
					oldBubbleComponents[adjacentOldBubbleRegion] = VISITED_BUBBLE;
				    }
				    else if (oldBubbleComponents[adjacentOldBubbleRegion] != VISITED_BUBBLE)
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
		std::cout << "Old bubble: " << oldBubbleRegion << " has been abandoned. Volume: " << bubbleVolumeHandle.get(oldBubbleRegion) << std::endl;
	    
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
    newBubbleRestVolumes.setSize(newBubbleRegionCount);
    newBubbleRestVolumes.constant(0);

    for (exint newBubbleRegion = 0; newBubbleRegion < newBubbleRegionCount; ++newBubbleRegion)
    {
	exint localComponent = newBubbleComponents[newBubbleRegion];

	if (localComponent >= 0)
	{
	    SolveReal newBubbleVolume = SolveReal(newBubbleRegionSizes[newBubbleRegion]);

	    assert(combinedNewVolumes[localComponent] > 0);
	    newBubbleRestVolumes[newBubbleRegion] = combinedOldVolumes[localComponent] * newBubbleVolume / combinedNewVolumes[localComponent]; 
	}
	else
	{
	    std::cout << "    New bubble region: " << newBubbleRegion << " has no old bubble connection and cannot preserve rest volume" << std::endl;
	    areBubblesUninitialized[newBubbleRegion] = true;
	}
    }
}

void
HDK_TwoPhasePressureSolver::buildBubbleRestVolumes(UT_Array<SolveReal> &newBubbleRestVolumes,
						    SIM_Object *obj,
						    const UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
						    const UT_Array<bool> &areBubblesTracked,
						    const UT_Array<exint> &newBubbleRegionSizes,
						    const exint bubbleRegionCount)
{
    std::cout << "\n// Build bubble rest volumes" << std::endl;
    UT_PerfMonAutoSolveEvent event(this, "Build bubble rest volumes");

    SIM_GeometryCopy *bubbleVolumeGeometry = getOrCreateGeometry(obj, "bubbleVolumeGeometry");
    UT_UniquePtr<SIM_GeometryAutoWriteLock> autoLockBubbleVolume = UTmakeUnique<SIM_GeometryAutoWriteLock>(bubbleVolumeGeometry, SIM_DATA_ID_PRESERVE);
    GU_Detail *bubbleVolumeDetail = &autoLockBubbleVolume->getGdp();

    assert(bubbleVolumeDetail != nullptr);

    {
	GA_RWHandleF bubbleRestVolumeHandle(bubbleVolumeDetail, GA_ATTRIB_POINT, "bubbleRestVolume");

	if (!bubbleRestVolumeHandle.isValid())
	{
	    bubbleVolumeDetail->addFloatTuple(GA_ATTRIB_POINT, "bubbleRestVolume", 1, GA_Defaults(-1));
	    bubbleRestVolumeHandle = GA_RWHandleF(bubbleVolumeDetail, GA_ATTRIB_POINT, "bubbleRestVolume");
	    bubbleRestVolumeHandle.bumpDataId();
	}

	UT_Array<bool> areBubblesUninitialized;
	areBubblesUninitialized.setSize(bubbleRegionCount);
	areBubblesUninitialized.constant(false);

	newBubbleRestVolumes.setSize(bubbleRegionCount);
	newBubbleRestVolumes.constant(0);

	distributeOldBubbleVolumes(newBubbleRestVolumes,
				    areBubblesUninitialized,
				    bubbleRestVolumeHandle,
				    oldToNewBubbleMap,
				    areBubblesTracked,
				    newBubbleRegionSizes);

	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (areBubblesUninitialized[bubbleRegion])
		newBubbleRestVolumes[bubbleRegion] = SolveReal(newBubbleRegionSizes[bubbleRegion]);
	}

	// Correct for any global volume drift
	GA_RWHandleF globalBubbleRestVolumeHandle(bubbleVolumeDetail, GA_ATTRIB_POINT, "globalBubbleRestVolume");

	if (!globalBubbleRestVolumeHandle.isValid())
	{
	    bubbleVolumeDetail->addFloatTuple(GA_ATTRIB_POINT, "globalBubbleRestVolume", 1, GA_Defaults(-1));
	    globalBubbleRestVolumeHandle = GA_RWHandleF(bubbleVolumeDetail, GA_ATTRIB_POINT, "globalBubbleRestVolume");
	    globalBubbleRestVolumeHandle.bumpDataId();
	}

	const exint oldBubbleRegionCount = bubbleVolumeDetail->getNumPoints();
	assert(oldToNewBubbleMap.size() == oldBubbleRegionCount);

	SolveReal globalBubbleVolume;
	if (oldBubbleRegionCount > 0)
	    globalBubbleVolume = globalBubbleRestVolumeHandle.get(0);
	else
	{
	    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    {
		if (areBubblesUninitialized[bubbleRegion])
		    globalBubbleVolume += SolveReal(newBubbleRestVolumes[bubbleRegion]);
	    }
	}

	assert(globalBubbleVolume > 0);
	
	SolveReal currentRestVolume = 0;
	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (areBubblesTracked[bubbleRegion])
		currentRestVolume += SolveReal(newBubbleRestVolumes[bubbleRegion]);
	}

	SolveReal volumeDiff = globalBubbleVolume - currentRestVolume;
	std::cout << "  Bubble region correction for global rest volume: " << volumeDiff << std::endl;

	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (areBubblesTracked[bubbleRegion])
		newBubbleRestVolumes[bubbleRegion] *= globalBubbleVolume / currentRestVolume;
	}

	// Debug check
	currentRestVolume = 0;
	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (areBubblesTracked[bubbleRegion])
		currentRestVolume += SolveReal(newBubbleRestVolumes[bubbleRegion]);
	}

	volumeDiff = globalBubbleVolume - currentRestVolume;
	std::cout << "  Post bubble region correction for global rest volume: " << volumeDiff << std::endl;
    }

    {
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

	// Update bubble geometry with new volumes
	SolveReal currentRestVolume = 0;
	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	{
	    if (areBubblesTracked[bubbleRegion])
		currentRestVolume += SolveReal(newBubbleRestVolumes[bubbleRegion]);
	}

	GA_Offset pointOffset = bubbleVolumeDetail->appendPointBlock(bubbleRegionCount);

	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion, ++pointOffset)
	{
	    if (areBubblesTracked[bubbleRegion])
		bubbleRestVolumeHandle.set(pointOffset, newBubbleRestVolumes[bubbleRegion]);

	    globalBubbleRestVolumeHandle.set(pointOffset, currentRestVolume);
	}
    }
}

void
HDK_TwoPhasePressureSolver::printCellGeometry(SIM_Object *obj,
						const SIM_RawIndexField &activeCellIndices,
						const SIM_RawIndexField &materialCellLabels,
						const UT_Vector3 origin,
						const fpreal dx)
{
    using SIM::FieldUtils::getFieldValue;

    SIM_GeometryCopy *cellActivityGeo = getOrCreateGeometry(obj, "cellActivityGeometry");

    SIM_GeometryAutoWriteLock autoLockCellActivity(cellActivityGeo, SIM_DATA_ID_PRESERVE);
    GU_Detail *cellActivityDetail = &autoLockCellActivity.getGdp();

    cellActivityDetail->clear();
    GA_RWHandleF cellSizeHandle(cellActivityDetail, GA_ATTRIB_POINT, "pscale");

    if (!cellSizeHandle.isValid())
    {
	cellActivityDetail->addFloatTuple(GA_ATTRIB_POINT, "pscale", 1, GA_Defaults(0));
	cellSizeHandle = GA_RWHandleF(cellActivityDetail, GA_ATTRIB_POINT, "pscale");
	cellSizeHandle.bumpDataId();
    }

    GA_RWHandleI cellActivityHandle(cellActivityDetail, GA_ATTRIB_POINT, "cellactivity");
    if (!cellActivityHandle.isValid())
    {
	cellActivityDetail->addIntTuple(GA_ATTRIB_POINT, "cellactivity", 1, GA_Defaults(-1));
	cellActivityHandle = GA_RWHandleI(cellActivityDetail, GA_ATTRIB_POINT, "cellactivity");
	cellActivityHandle.bumpDataId();
    }

    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(materialCellLabels.field());

    UT_VoxelTileIteratorI vitt;
    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
	if (boss->opInterrupt())
	    break;
	
	if (!vit.isTileConstant() ||
	    vit.getValue() == MaterialLabels::LIQUID_CELL ||
	    vit.getValue() == MaterialLabels::BUBBLE_CELL)
	{
	    vitt.setTile(vit);

	    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
	    {
		if (vitt.getValue() == MaterialLabels::LIQUID_CELL ||
		    vitt.getValue() == MaterialLabels::BUBBLE_CELL)
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

		    GA_Offset pointOffset = cellActivityDetail->appendPoint();
		    cellSizeHandle.set(pointOffset, dx);

		    if (vitt.getValue() == MaterialLabels::LIQUID_CELL)
		    {
			assert(getFieldValue(activeCellIndices, cell) >= 0);
			cellActivityHandle.set(pointOffset, 0);
		    }
		    else if (vitt.getValue() == MaterialLabels::BUBBLE_CELL)
		    {
			assert(getFieldValue(activeCellIndices, cell) >= 0);
			cellActivityHandle.set(pointOffset, 1);
		    }
		    else assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);

		    UT_Vector3 cellCenter = origin + dx * UT_Vector3(vitt.x() + .5, vitt.y() + .5, vitt.z() + .5);
		    cellActivityDetail->setPos3(pointOffset, cellCenter);
		}
	    }
	}
    }

    cellSizeHandle.bumpDataId();
    cellActivityHandle.bumpDataId();
    cellActivityDetail->getAttributes().bumpAllDataIds(GA_ATTRIB_POINT);
}

void
HDK_TwoPhasePressureSolver::buildRHS(Vector &rhsVector,
					    const SIM_RawIndexField &materialCellLabels,
					    const SIM_RawIndexField &activeCellIndices,
					    const SIM_VectorField &velocity,
					    const std::array<const SIM_RawField *, 3> &cutCellWeights,
					    const SIM_VectorField *solidVelocity,
					    const SolveReal dx) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(activeCellIndices.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(activeCellIndices.field());

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
		    UT_Vector3I cell(vit.x(), vit.y(), vit.z());

		    exint index = vit.getValue();
		    if (index >= 0)
		    {
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::LIQUID_CELL ||
				getFieldValue(materialCellLabels, cell) == MaterialLabels::BUBBLE_CELL);

			SolveReal divergence = 0;

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I face = cellToFaceMap(cell, axis, direction);
				SolveReal weight = getFieldValue(*cutCellWeights[axis], face);

				SolveReal sign = (direction == 0) ? 1 : -1;

				if (weight > 0)
				    divergence += sign * weight * getFieldValue(*velocity.getField(axis), face);

				if (solidVelocity != nullptr && weight < 1)
				{
				    UT_Vector3 point;
				    velocity.getField(axis)->indexToPos(face[0], face[1], face[2], point);

				    divergence += sign * (1. - weight) * solidVelocity->getField(axis)->getValue(point);
				}
			    }

			rhsVector(index) = divergence;
		    }
		    else assert(getFieldValue(materialCellLabels, cell) != MaterialLabels::LIQUID_CELL &&
				getFieldValue(materialCellLabels, cell) != MaterialLabels::BUBBLE_CELL);
		}
	    }
	}
    });
}

void
HDK_TwoPhasePressureSolver::applyGoalDivergence(Vector &rhsVector,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_RawIndexField &activeCellIndices,
						    const SIM_RawField &goalDivergence) const
{
    using SIM::FieldUtils::getFieldValue;

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

			exint index = getFieldValue(activeCellIndices, cell);
			assert(index >= 0);

			rhsVector(index) += getFieldValue(goalDivergence, cell);
		    }
		}
	    }
	}
    });
}

void
HDK_TwoPhasePressureSolver::buildActiveRegionDivergence(UT_Array<UT_Array<SolveReal>> &parallelActiveRegionDivergence,
							    UT_Array<UT_Array<SolveReal>> &parallelActiveRegionCellCount,
							    const SIM_RawIndexField &activeCellIndices,
							    const SIM_RawIndexField &activeRegionIndices,
							    const Vector &rhsVector) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildActiveRegionDivergenceAlgorithm;
    buildActiveRegionDivergenceAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_Array<SolveReal> &localActiveRegionDivergence = parallelActiveRegionDivergence[info.job()];
	UT_Array<SolveReal> &localActiveRegionCellCount = parallelActiveRegionCellCount[info.job()];

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(activeCellIndices.field());
	vit.splitByTile(info);
	
	UT_VoxelTileIteratorI vitt;

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant())
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

		    exint index = vitt.getValue();
		    if (index >= 0)
		    {
			exint activeRegion = getFieldValue(activeRegionIndices, cell);
			assert(activeRegion >= 0);

			++localActiveRegionCellCount[activeRegion];
			localActiveRegionDivergence[activeRegion] += rhsVector(index);
		    }
		    else
		    { 
			assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			assert(getFieldValue(activeRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_TwoPhasePressureSolver::removeAverageDivergence(Vector &rhsVector,
							const SIM_RawIndexField &activeCellIndices,
							const SIM_RawIndexField &activeRegionIndices,
							const UT_Array<SolveReal> &averageActiveRegionDivergence) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(activeCellIndices.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(activeCellIndices.field());

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
		    UT_Vector3I cell(vit.x(), vit.y(), vit.z());

		    exint index = vit.getValue();
		    if (index >= 0)
		    {
			exint activeRegion = getFieldValue(activeRegionIndices, cell);
			assert(activeRegion >= 0);

			rhsVector(index) -= averageActiveRegionDivergence[activeRegion];
		    }
		    else
		    { 
			assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			assert(getFieldValue(activeRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
		    }
		}
	    }
	}
    });
}

void
HDK_TwoPhasePressureSolver::removeBackgroundNullSpace(Vector &rhsVector,
							const SIM_RawIndexField &activeCellIndices,
							const SIM_RawIndexField &activeRegionIndices,
							const exint activeRegionCount) const
{
    const int threadCount = UT_Thread::getNumProcessors();

    // Compute divergence for each active region
    UT_Array<UT_Array<SolveReal>> parallelActiveRegionDivergence;
    parallelActiveRegionDivergence.setSize(threadCount);

    UT_Array<UT_Array<SolveReal>> parallelActiveRegionCellCount;
    parallelActiveRegionCellCount.setSize(threadCount);

    for (int thread = 0; thread < threadCount; ++thread)
    {
	parallelActiveRegionDivergence[thread].setSize(activeRegionCount);
	parallelActiveRegionDivergence[thread].constant(0);

	parallelActiveRegionCellCount[thread].setSize(activeRegionCount);
	parallelActiveRegionCellCount[thread].constant(0);
    }

    buildActiveRegionDivergence(parallelActiveRegionDivergence,
				parallelActiveRegionCellCount,
				activeCellIndices,
				activeRegionIndices,
				rhsVector);

    UT_Array<SolveReal> averageActiveRegionDivergence;
    averageActiveRegionDivergence.setSize(activeRegionCount);
    averageActiveRegionDivergence.constant(0);

    UT_Array<SolveReal> activeRegionCellCount;
    activeRegionCellCount.setSize(activeRegionCount);
    activeRegionCellCount.constant(0);

    for (int thread = 0; thread < threadCount; ++thread)
	for (exint activeRegion = 0; activeRegion < activeRegionCount; ++activeRegion)
	{
	    averageActiveRegionDivergence[activeRegion] += parallelActiveRegionDivergence[thread][activeRegion];
	    activeRegionCellCount[activeRegion] += parallelActiveRegionCellCount[thread][activeRegion];
	}

    for (exint activeRegion = 0; activeRegion < activeRegionCount; ++activeRegion)
    {
	assert(activeRegionCellCount[activeRegion] > 0);
	averageActiveRegionDivergence[activeRegion] /= activeRegionCellCount[activeRegion];
    }

    removeAverageDivergence(rhsVector,
			    activeCellIndices,
			    activeRegionIndices,
			    averageActiveRegionDivergence);

    std::cout << "    Average divergence in rhs vector after latent divergence removal: " << rhsVector.sum() / SolveReal(rhsVector.rows()) << std::endl;
}

void
HDK_TwoPhasePressureSolver::setCorrectedBubbleDivergence(UT_Array<UT_Array<SolveReal>> &parallelCorrectedBubbleDivergence,
							    UT_Array<UT_Array<exint>> &parallelUncorrectedCellCount,
							    UT_Array<UT_Array<exint>> &parallelCorrectedBubbleCellCount,
							    Vector &rhsVector,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &activeCellIndices,
							    const SIM_RawIndexField &activeRegionIndices,
							    const SIM_RawIndexField &bubbleRegionIndices,
							    const UT_Array<bool> &areBubblesTracked,
							    const UT_Array<SolveReal> &newBubbleRestVolumes,
							    const UT_Array<exint> &bubbleRegionSizes,
							    const SolveReal dx,
							    const SolveReal dt,
							    const bool doCorrectVolumeDrift) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    constexpr SolveReal correctionScale = .1;

    const exint bubbleRegionCount = areBubblesTracked.size();

    UT_Array<SolveReal> fractionalBubbleVolumeCorrection;
    fractionalBubbleVolumeCorrection.setSize(bubbleRegionCount);
    fractionalBubbleVolumeCorrection.constant(0);

    if (doCorrectVolumeDrift)
    {
	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    fractionalBubbleVolumeCorrection[bubbleRegion] = (newBubbleRestVolumes[bubbleRegion] - SolveReal(bubbleRegionSizes[bubbleRegion])) / SolveReal(bubbleRegionSizes[bubbleRegion]);
    }

    UT_ThreadedAlgorithm setUntrackedDivergenceAlgorithm;
    setUntrackedDivergenceAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<SolveReal> &localCorrectedBubbleDivergence = parallelCorrectedBubbleDivergence[info.job()];
	UT_Array<exint> &localUncorrectedCellCount = parallelUncorrectedCellCount[info.job()];
	UT_Array<exint> &localCorrectedBubbleCellCount = parallelCorrectedBubbleCellCount[info.job()];

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
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

		    if (vitt.getValue() == MaterialLabels::LIQUID_CELL)
		    {
			exint activeRegion = getFieldValue(activeRegionIndices, cell);
			assert(activeRegion >= 0);

			++localUncorrectedCellCount[activeRegion];
		    }
		    else if (vitt.getValue() == MaterialLabels::BUBBLE_CELL)
		    {
			exint bubbleRegion = getFieldValue(bubbleRegionIndices, cell);
			assert(bubbleRegion >= 0);

			exint activeRegion = getFieldValue(activeRegionIndices, cell);
			assert(activeRegion >= 0);

			if (areBubblesTracked[bubbleRegion])
			{
			    if (doCorrectVolumeDrift)
			    {
				SolveReal localDivergence = dx * fractionalBubbleVolumeCorrection[bubbleRegion] / dt;

				localCorrectedBubbleDivergence[activeRegion] += correctionScale * localDivergence;

				exint index = getFieldValue(activeCellIndices, cell);
				assert(index >= 0);

				rhsVector(index) += correctionScale * localDivergence;

				++localCorrectedBubbleCellCount[bubbleRegion];
			    }
			    else
				++localUncorrectedCellCount[activeRegion];
			}
			else
			{
			    SolveReal localDivergence = -dx / dt;
			    localCorrectedBubbleDivergence[activeRegion] += localDivergence;

			    exint index = getFieldValue(activeCellIndices, cell);
			    assert(index >= 0);

			    rhsVector(index) += localDivergence;

			    ++localCorrectedBubbleCellCount[bubbleRegion];
			}
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_TwoPhasePressureSolver::removeCorrectionDivergenceNullSpace(Vector &rhsVector,
								const SIM_RawIndexField &materialCellLabels,
								const SIM_RawIndexField &activeCellIndices,
								const SIM_RawIndexField &activeRegionIndices,
								const SIM_RawIndexField &bubbleRegionIndices,
								const UT_Array<bool> &areBubblesTracked,
								const UT_Array<SolveReal> &activeRegionAverageDivergence,
								const bool doCorrectVolumeDrift) const
{
    using SIM::FieldUtils::getFieldValue;

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

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::LIQUID_CELL ||
		vit.getValue() == MaterialLabels::BUBBLE_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    UT_Vector3I cell(vit.x(), vit.y(), vit.z());

		    if (vit.getValue() == MaterialLabels::LIQUID_CELL)
		    {
			exint index = getFieldValue(activeCellIndices, cell);
			assert(index >= 0);

			exint activeRegion = getFieldValue(activeRegionIndices, cell);
			assert(activeRegion >= 0);

			rhsVector(index) -= activeRegionAverageDivergence[activeRegion];
		    }
		    else if (vit.getValue() == MaterialLabels::BUBBLE_CELL)
		    {
			exint bubbleRegion = getFieldValue(bubbleRegionIndices, cell);
			assert(bubbleRegion >= 0);

			if (!doCorrectVolumeDrift && areBubblesTracked[bubbleRegion])
			{
			    exint activeRegion = getFieldValue(activeRegionIndices, cell);
			    assert(activeRegion >= 0);

			    exint index = getFieldValue(activeCellIndices, cell);
			    assert(index >= 0);

			    rhsVector(index) -= activeRegionAverageDivergence[activeRegion];
			}
		    }
		}
	    }
	}
    });
}

void
HDK_TwoPhasePressureSolver::applyBubbleCorrection(Vector &rhsVector,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_RawIndexField &activeCellIndices,
						    const SIM_RawIndexField &activeRegionIndices,
						    const SIM_RawIndexField &bubbleRegionIndices,
						    const UT_Array<SolveReal> &newBubbleRestVolumes,
						    const UT_Array<exint> &bubbleRegionSizes,
						    const UT_Array<bool> &areBubblesTracked,
						    const SolveReal dx,
						    const SolveReal dt,
						    const bool doTrackBubbleVolumes,
						    const bool doCorrectBubbleVolumeDrift,
						    const exint activeRegionCount,
						    const exint bubbleRegionCount) const
{
    const int threadCount = UT_Thread::getNumProcessors();

    UT_Array<UT_Array<SolveReal>> parallelCorrectedBubbleDivergence;
    parallelCorrectedBubbleDivergence.setSize(threadCount);

    UT_Array<UT_Array<exint>> parallelUncorrectedCellCount;
    parallelUncorrectedCellCount.setSize(threadCount);

    UT_Array<UT_Array<exint>> parallelCorrectedBubbleCellCount;
    parallelCorrectedBubbleCellCount.setSize(threadCount);

    for (int thread = 0; thread < threadCount; ++thread)
    {
	parallelCorrectedBubbleDivergence[thread].setSize(activeRegionCount);
	parallelCorrectedBubbleDivergence[thread].constant(0);

	parallelUncorrectedCellCount[thread].setSize(activeRegionCount);
	parallelUncorrectedCellCount[thread].constant(0);

	parallelCorrectedBubbleCellCount[thread].setSize(bubbleRegionCount);
	parallelCorrectedBubbleCellCount[thread].constant(0);
    }

    // Apply correction for untracked bubbles and volume correction
    setCorrectedBubbleDivergence(parallelCorrectedBubbleDivergence,
				    parallelUncorrectedCellCount,
				    parallelCorrectedBubbleCellCount,
				    rhsVector,
				    materialCellLabels,
				    activeCellIndices,
				    activeRegionIndices,
				    bubbleRegionIndices,
				    areBubblesTracked,
				    newBubbleRestVolumes,
				    bubbleRegionSizes,
				    dx, dt,
				    doCorrectBubbleVolumeDrift);

    // Accumulate correction divergence and uncorrected cell count for each active region

    UT_Array<SolveReal> activeRegionAverageDivergence;
    activeRegionAverageDivergence.setSize(activeRegionCount);
    activeRegionAverageDivergence.constant(0);

    UT_Array<exint> uncorrectedCellCount;
    uncorrectedCellCount.setSize(activeRegionCount);
    uncorrectedCellCount.constant(0);

    for (int thread = 0; thread < threadCount; ++thread)
	for (exint activeRegion = 0; activeRegion < activeRegionCount; ++activeRegion)
	{
	    activeRegionAverageDivergence[activeRegion] += parallelCorrectedBubbleDivergence[thread][activeRegion];
	    uncorrectedCellCount[activeRegion] += parallelUncorrectedCellCount[thread][activeRegion];
	}
    
    for (exint activeRegion = 0; activeRegion < activeRegionCount; ++activeRegion)
	activeRegionAverageDivergence[activeRegion] /= SolveReal(uncorrectedCellCount[activeRegion]);

    UT_Array<exint> correctedBubbleCellCount;
    correctedBubbleCellCount.setSize(bubbleRegionCount);
    correctedBubbleCellCount.constant(0);

    for (int thread = 0; thread < threadCount; ++thread)
	for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
	    correctedBubbleCellCount[bubbleRegion] += parallelCorrectedBubbleCellCount[thread][bubbleRegion];

    for (exint bubbleRegion = 0; bubbleRegion < bubbleRegionCount; ++bubbleRegion)
    {
	if (!areBubblesTracked[bubbleRegion])
	    std::cout << "    Untracked bubble region " << bubbleRegion << " of size " << correctedBubbleCellCount[bubbleRegion] << std::endl;
	else if (doTrackBubbleVolumes)
	{
	    std::cout << " Bubble: " << bubbleRegion << std::endl;
	    std::cout << " Current bubble size: " << bubbleRegionSizes[bubbleRegion] << std::endl;
	    std::cout << " Rest volume: " << newBubbleRestVolumes[bubbleRegion] << std::endl;

	    if (doCorrectBubbleVolumeDrift)
		std::cout << " Debug check corrected bubble cell count: " << correctedBubbleCellCount[bubbleRegion] << std::endl;
	}
    }

    // Remove divergence in tracked cells
    removeCorrectionDivergenceNullSpace(rhsVector,
					materialCellLabels,
					activeCellIndices,
					activeRegionIndices,
					bubbleRegionIndices,
					areBubblesTracked,
					activeRegionAverageDivergence,
					doCorrectBubbleVolumeDrift);

    std::cout << "    Average divergence in rhs vector after untracked bubble removal: " << rhsVector.sum() / SolveReal(rhsVector.rows()) << std::endl;
}

void
HDK_TwoPhasePressureSolver::buildLinearSystem(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelPoissonElements,
						std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelDiagonalPrecondElements,
						const SIM_RawIndexField &materialCellLabels,
						const SIM_RawIndexField &activeCellIndices,
						const std::array<const SIM_RawField *, 3> &cutCellWeights,
						const SIM_RawField &liquidSurface,
						const SolveReal liquidDensity,
						const SolveReal bubbleDensity) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildLinearSystemAlgorithm;
    buildLinearSystemAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(activeCellIndices.field());
	vit.splitByTile(info);
	
	UT_VoxelTileIteratorI vitt;

	std::vector<Eigen::Triplet<SolveReal>> &localPoissonElements = parallelPoissonElements[info.job()];
	std::vector<Eigen::Triplet<SolveReal>> &localDiagonalPrecondElements = parallelDiagonalPrecondElements[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant())
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    exint index = vitt.getValue();
		    if (index >= 0)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			SolveReal diagonal = 0;

			exint materialLabel = getFieldValue(materialCellLabels, cell);
			assert(materialLabel == MaterialLabels::LIQUID_CELL || materialLabel == MaterialLabels::BUBBLE_CELL);

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I face = cellToFaceMap(cell, axis, direction);
				SolveReal weight = getFieldValue(*cutCellWeights[axis], face);

				if (weight > 0)
				{
				    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= activeCellIndices.getVoxelRes()[axis])
					continue;

				    exint adjacentIndex = getFieldValue(activeCellIndices, adjacentCell);
				    assert(adjacentIndex >= 0);

				    exint adjacentMaterialLabel = getFieldValue(materialCellLabels, adjacentCell);

				    assert(adjacentMaterialLabel == MaterialLabels::LIQUID_CELL || adjacentMaterialLabel == MaterialLabels::BUBBLE_CELL);

				    if (materialLabel != adjacentMaterialLabel)
				    {
					SolveReal phi0 = getFieldValue(liquidSurface, cell);
					SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

					SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
					theta = SYSclamp(theta, SolveReal(0), SolveReal(1));

					weight /= (theta * liquidDensity + (1. - theta) * bubbleDensity);
				    }
				    else
				    {
					if (materialLabel == MaterialLabels::LIQUID_CELL)
					    weight /= liquidDensity;
					else
					    weight /= bubbleDensity;
				    }

				    diagonal += weight;
				    localPoissonElements.push_back(Eigen::Triplet<SolveReal>(index, adjacentIndex, -weight));
				}
			    }

			assert(diagonal > 0);
			localPoissonElements.push_back(Eigen::Triplet<SolveReal>(index, index, diagonal));
			localDiagonalPrecondElements.push_back(Eigen::Triplet<SolveReal>(index, index, 1. / diagonal));
		    }
		}
	    }
	}

	return 0;
    });    
}

void
HDK_TwoPhasePressureSolver::buildValidFaces(SIM_VectorField &validFaces,
						const SIM_RawIndexField &materialCellLabels,
						const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    // Uncompress valid face tiles
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
						[](const exint label){ return label == MaterialLabels::LIQUID_CELL ||
										label == MaterialLabels::BUBBLE_CELL; },
						axis);

	HDK::Utilities::uncompressTiles(*validFaces.getField(axis), isTileOccupiedList);

	HDK::Utilities::classifyValidFaces(*validFaces.getField(axis),
					    materialCellLabels,
					    *cutCellWeights[axis],
					    [](const exint label){ return label == MaterialLabels::LIQUID_CELL ||
									    label == MaterialLabels::BUBBLE_CELL; },
					    axis);
    }
}

void
HDK_TwoPhasePressureSolver::copyPressureVectorToGrid(SIM_RawField &pressure,
							const SIM_RawIndexField &activeCellIndices,
							const Vector &solutionVector) const
{
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(activeCellIndices.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(activeCellIndices.field());
        
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
		    exint activeIndex = vit.getValue();
		    if (activeIndex >= 0)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());
			setFieldValue(pressure, cell, solutionVector[activeIndex]);
		    }
		}
	    }
	}
    });
}

void
HDK_TwoPhasePressureSolver::applySolutionToVelocity(SIM_RawField &velocity,
							const SIM_RawField &validFaces,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &activeCellIndices,
							const SIM_RawField &pressure,
							const SIM_RawField &liquidSurface,
							const SIM_RawField &cutCellWeights,
							const SolveReal liquidDensity,
							const SolveReal bubbleDensity,
							const int axis,
							const SolveReal dx) const
{
    using SIM::FieldUtils::faceToCellMap;
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

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

	    if (!vit.isTileConstant() || vit.getValue() == HDK::Utilities::VALID_FACE)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    UT_Vector3I face(vit.x(), vit.y(), vit.z());

		    if (vit.getValue() == HDK::Utilities::VALID_FACE)
		    {
			UT_Vector3I backwardCell = faceToCellMap(face, axis, 0);
			UT_Vector3I forwardCell = faceToCellMap(face, axis, 1);

			assert(backwardCell[axis] >= 0 && forwardCell[axis] < activeCellIndices.getVoxelRes()[axis]);

			exint backwardMaterial = getFieldValue(materialCellLabels, backwardCell);
			exint forwardMaterial = getFieldValue(materialCellLabels, forwardCell);

			SolveReal localDensity;

			if (backwardMaterial == MaterialLabels::LIQUID_CELL && forwardMaterial == MaterialLabels::LIQUID_CELL)
			    localDensity = liquidDensity;
			else if (backwardMaterial == MaterialLabels::BUBBLE_CELL && forwardMaterial == MaterialLabels::BUBBLE_CELL)
			    localDensity = bubbleDensity;
			else
			{
			    assert((backwardMaterial == MaterialLabels::LIQUID_CELL && forwardMaterial == MaterialLabels::BUBBLE_CELL) ||
				    (backwardMaterial == MaterialLabels::BUBBLE_CELL && forwardMaterial == MaterialLabels::LIQUID_CELL));

			    SolveReal phi0 = getFieldValue(liquidSurface, backwardCell);
			    SolveReal phi1 = getFieldValue(liquidSurface, forwardCell);

			    SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
			    theta = SYSclamp(theta, SolveReal(0), SolveReal(1));
			    localDensity = theta * liquidDensity + (1. - theta) * bubbleDensity;
			}

			SolveReal gradient = getFieldValue(pressure, forwardCell) - getFieldValue(pressure, backwardCell);
			gradient /= localDensity;

			setFieldValue(velocity, face, getFieldValue(velocity, face) - gradient);
		    }
		    else assert(getFieldValue(cutCellWeights, face) == 0);
		}
	    }
	}
    });
}

void
HDK_TwoPhasePressureSolver::computeResultingDivergence(std::array<SolveReal, 3> &divergenceStats,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_VectorField &velocity,
							const std::array<const SIM_RawField *, 3> &cutCellWeights,
							const SIM_VectorField *solidVelocity) const
{
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    const int threadCount = UT_Thread::getNumProcessors();

    UT_Array<SolveReal> parallelAccumulatedDivergence;
    parallelAccumulatedDivergence.setSize(threadCount);
    parallelAccumulatedDivergence.constant(0);

    UT_Array<SolveReal> parallelMaxDivergence;
    parallelMaxDivergence.setSize(threadCount);
    parallelMaxDivergence.constant(0);

    UT_Array<SolveReal> parallelCellCount;
    parallelCellCount.setSize(threadCount);
    parallelCellCount.constant(0);

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm computeResultingDivergenceAlgorithm;
    computeResultingDivergenceAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	SolveReal &localAccumulatedDivergence = parallelAccumulatedDivergence[info.job()];
	SolveReal &localMaxDivergence = parallelMaxDivergence[info.job()];
	SolveReal &localCellCount = parallelCellCount[info.job()];

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
		    if (vitt.getValue() == MaterialLabels::LIQUID_CELL ||
			vitt.getValue() == MaterialLabels::BUBBLE_CELL)
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

			//divergence *= sqr_dx;

			localAccumulatedDivergence += divergence;
			localMaxDivergence = std::max(localMaxDivergence, std::fabs(divergence));
			++localCellCount;
		    }
		}
	    }
	}

	return 0;
    });

    divergenceStats = {0};

    for (int thread = 0; thread < threadCount; ++thread)
    {
	divergenceStats[0] += parallelAccumulatedDivergence[thread];
	divergenceStats[1] = std::max(divergenceStats[1], parallelMaxDivergence[thread]);
	divergenceStats[2] += parallelCellCount[thread];
    }

}