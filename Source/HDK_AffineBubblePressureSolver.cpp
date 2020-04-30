#include "HDK_AffineBubblePressureSolver.h"

#include <Eigen/LU>

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
#include <UT/UT_PerfMonAutoEvent.h>

#include <UT/UT_StopWatch.h>

void initializeSIM(void *)
{
    IMPLEMENT_DATAFACTORY(HDK_AffineBubblePressureSolver);
}

// Standard constructor, note that BaseClass was crated by the
// DECLARE_DATAFACTORY and provides an easy way to chain through
// the class hierarchy.
HDK_AffineBubblePressureSolver::HDK_AffineBubblePressureSolver(const SIM_DataFactory *factory)
    : BaseClass(factory)
{
}

HDK_AffineBubblePressureSolver::~HDK_AffineBubblePressureSolver()
{
}

const SIM_DopDescription* HDK_AffineBubblePressureSolver::getDopDescription()
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

    static PRM_Name 	theExteriorBandSizeName("exteriorBandwidth", "Exteroir Bandwidth");
    static PRM_Default 	theExteriorBandSizeDefault(3);

    static PRM_Name 	theUseTiledInteriorName("useTiledInterior", "Use Tiled Interior");

    static PRM_Name 	theTiledInteriorSizeName("tileSize", "Tiled Interior Size");
    static PRM_Default 	theTiledInteriorSizeDefault(16);

    static PRM_Name	theApplyAffineToLiquidsName("applyAffineToLiquids", "Apply Affine to Liquids");

    static PRM_Name	theUseMGPreconditionerName("useMGPreconditioner", "Use Multigrid Preconditioner");

    static PRM_Name	theTrackBubbleIDsName("trackBubbleIDs", "Track Bubble IDs");

    static PRM_Name	theGeometryName(GAS_NAME_GEOMETRY, "Bubble ID Tracking Geometry");
    static PRM_Default  theGeometryNameDefault(0, "Geometry");

    static PRM_Name	thePrintCellsName("printCells", "Print Cell Activity");
    static PRM_Name	theOnlyPrintCellsName("onlyPrintCells", "Only Print Cells");

    static PRM_Name    thePrintedCellsGeometryName("cellActivityGeometry", "Cell Activity Geometry");
    static PRM_Default thePrintedCellsGeometryDefault(0, "CellActivityGeometry");

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

	PRM_Template(PRM_TOGGLE, 1, &theUseMGPreconditionerName, PRMoneDefaults),

	PRM_Template(PRM_INT, 1, &theExteriorBandSizeName, &theExteriorBandSizeDefault),

	PRM_Template(PRM_TOGGLE, 1, &theUseTiledInteriorName, PRMzeroDefaults),

	PRM_Template(PRM_INT, 1, &theTiledInteriorSizeName, &theTiledInteriorSizeDefault),

	PRM_Template(PRM_TOGGLE, 1, &theApplyAffineToLiquidsName, PRMzeroDefaults),

        PRM_Template(PRM_TOGGLE, 1, &theTrackBubbleIDsName, PRMzeroDefaults),

        PRM_Template(PRM_STRING, 1, &theGeometryName, &theGeometryNameDefault),

	PRM_Template(PRM_TOGGLE, 1, &thePrintCellsName, PRMzeroDefaults),
	PRM_Template(PRM_TOGGLE, 1, &theOnlyPrintCellsName, PRMzeroDefaults),

	PRM_Template(PRM_STRING, 1, &thePrintedCellsGeometryName,
				    &thePrintedCellsGeometryDefault),

    	PRM_Template()
    };

    static SIM_DopDescription theDopDescription(true,
						"HDK_AffineBubblePressureSolver",
						"HDK Affine Bubble Pressure Solver",
						"$OS",
						classname(),
						theTemplates);

    setGasDescription(theDopDescription);

    return &theDopDescription;
}

bool
HDK_AffineBubblePressureSolver::solveGasSubclass(SIM_Engine &engine,
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

	materialCellLabels.match(liquidSurface);
	materialCellLabels.makeConstant(HDK::Utilities::SOLID_CELL);

	HDK::Utilities::buildMaterialCellLabels(materialCellLabels,
						liquidSurface,
						*solidSurface,
						cutCellWeights);

	// Set labels specific to this solver
	setMaterialCellLabels(materialCellLabels);
    }

    ////////////////////////////////////////////
    //
    // Build band of exterior bubble cells
    //
    ////////////////////////////////////////////

    {
	UT_PerfMonAutoSolveEvent event(this, "Build exterior band of active cells at bubble boundary");
	std::cout << "\n// Build exterior band of active cells at bubble boundary" << std::endl;

	buildExteriorCells(materialCellLabels,
			    cutCellWeights,
			    MaterialLabels::EXTERIOR_BUBBLE_CELL,
			    MaterialLabels::INTERIOR_BUBBLE_CELL,
			    MaterialLabels::INTERIOR_LIQUID_CELL);
    }

    ////////////////////////////////////////////
    //
    // Set exterior cells inside the bubble along
    // tile boundaries
    //
    ////////////////////////////////////////////

    if (getUseTiledInterior())
    {
	UT_PerfMonAutoSolveEvent event(this, "Build tile boundary active cells");
    	std::cout << "\n// Build tile boundary active cells" << std::endl;

	int tileSize = getTileSize();
	std::cout << "Tile size: " << tileSize << std::endl;
	buildTiledActiveCells(materialCellLabels,
				tileSize,
				MaterialLabels::EXTERIOR_BUBBLE_CELL,
				MaterialLabels::INTERIOR_BUBBLE_CELL);
    }    

    ////////////////////////////////////////////
    //
    // Build band of exterior liquid cells and
    // cells along tile boundaries
    //
    ////////////////////////////////////////////

    bool doApplyAffineToLiquid = getApplyAffineToLiquids();

    if (doApplyAffineToLiquid)
    {
	{
	    UT_PerfMonAutoSolveEvent event(this, "Build exterior band of active cells at liquid boundary");
	    std::cout << "\n// Build exterior band of active cells at liquid boundary" << std::endl;

	    buildExteriorCells(materialCellLabels,
				cutCellWeights,
				MaterialLabels::EXTERIOR_LIQUID_CELL,
				MaterialLabels::INTERIOR_LIQUID_CELL,
				MaterialLabels::EXTERIOR_BUBBLE_CELL);
	}

	if (getUseTiledInterior())
	{
	    UT_PerfMonAutoSolveEvent event(this, "Build tile boundary active liquid cells");
	    std::cout << "\n// Build tile boundary active liquid cells" << std::endl;

	    int tileSize = getTileSize();
	    std::cout << "Tile size: " << tileSize << std::endl;
	    buildTiledActiveCells(materialCellLabels,
				    tileSize,
				    MaterialLabels::EXTERIOR_LIQUID_CELL,
				    MaterialLabels::INTERIOR_LIQUID_CELL);
	}
    }
    else
    {
	UT_PerfMonAutoSolveEvent event(this, "Set all liquid cells to be active exterior cells");
	std::cout << "\n// Set all liquid cells to be active exterior cells" << std::endl;	

	setAllInteriorCellsExterior(materialCellLabels,
				    MaterialLabels::EXTERIOR_LIQUID_CELL,
				    MaterialLabels::INTERIOR_LIQUID_CELL);
    }

    ////////////////////////////////////////////
    //
    // Any untracked bubble need to be set as fully
    // EXTERIOR so that we can collapse them away
    //
    ////////////////////////////////////////////

    bool doTrackBubbleIDs = getTrackBubbleIDs();
    
    UT_Array<bool> isBubbleTracked;

    SIM_RawIndexField fullBubbleRegionIndices;

    if (doTrackBubbleIDs)
    {
	UT_PerfMonAutoSolveEvent event(this, "Track bubble IDs");
    	std::cout << "\n// Track bubble IDs" << std::endl;

	buildBubbleRegions(fullBubbleRegionIndices, isBubbleTracked, materialCellLabels, obj, liquidSurface, cutCellWeights);
    }

    ////////////////////////////////////////////
    //
    // Build connected components for each interior region
    //
    ////////////////////////////////////////////

    SIM_RawIndexField interiorRegionIndices;
    exint interiorRegionCount = 0;
    {
	UT_PerfMonAutoSolveEvent event(this, "Build interior region labels");
    	std::cout << "\n// Build interior region labels" << std::endl;

	interiorRegionIndices.match(liquidSurface);

	SIM_VolumetricConnectedComponentBuilder interiorRegionBuilder(interiorRegionIndices, materialCellLabels, cutCellWeights.data());

	interiorRegionCount = interiorRegionBuilder.buildConnectedComponents([](const exint label) { return label == MaterialLabels::INTERIOR_BUBBLE_CELL ||
													    label == MaterialLabels::INTERIOR_LIQUID_CELL; });

	HDK::Utilities::overwriteIndices(interiorRegionIndices,
					    SIM_VolumetricConnectedComponentBuilder::INACTIVE_REGION,
					    HDK::Utilities::UNLABELLED_CELL);

	std::cout << "    Interior region count: " << interiorRegionCount << std::endl;
    }

    ////////////////////////////////////////////
    //
    // Build bounding box for each interior region
    // and remove any interior regions with a degenerate
    // bounding box
    //
    ////////////////////////////////////////////

    {
	UT_PerfMonAutoSolveEvent event(this, "Remove single cell regions");
	std::cout << "\n// Remove single cell regions" << std::endl;

	interiorRegionCount = removeDegenerateCellRegions(interiorRegionIndices, materialCellLabels, interiorRegionCount);	
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

	activeCellIndices.match(liquidSurface);
	activeCellIndices.makeConstant(HDK::Utilities::UNLABELLED_CELL);
	activeCellCount = buildActiveCellIndices(activeCellIndices, materialCellLabels);
    }

    ////////////////////////////////////////////
    //
    // Print geometry for display and debugging
    //
    ////////////////////////////////////////////

    if (getPrintCells())
    {
	using SIM::FieldUtils::getFieldValue;
	using SIM::FieldUtils::setFieldValue;

	std::cout << "// Printing cell geometry" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Printing cell geometry");

	printCellGeometry(obj,
			    activeCellIndices,
			    interiorRegionIndices,
			    materialCellLabels,
			    *velocity,
			    dx);

        if (getOnlyPrintCells())
	    return true;
    }

    ////////////////////////////////////////////
    //
    // Build COM for each interior region
    //
    ////////////////////////////////////////////

    const int threadCount = UT_Thread::getNumProcessors();

    UT_Array<UT_Vector3SR> interiorRegionCOM;
    interiorRegionCOM.setSize(interiorRegionCount);
    interiorRegionCOM.constant(UT_Vector3SR(0,0,0));

    exint tbbGrainSize = interiorRegionCount / (4 * threadCount);

    {
	UT_PerfMonAutoSolveEvent event(this, "Build interior region COM");
	std::cout << "\n// Build interior region COM" << std::endl;

	UT_Array<UT_Array<UT_Vector3SR>> parallelInteriorRegionCOM;
	parallelInteriorRegionCOM.setSize(threadCount);

	UT_Array<UT_Array<SolveReal>> parallelInteriorRegionCellCount;
	parallelInteriorRegionCellCount.setSize(threadCount);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    parallelInteriorRegionCOM[thread].setSize(interiorRegionCount);
	    parallelInteriorRegionCOM[thread].constant(UT_Vector3SR(0,0,0));

	    parallelInteriorRegionCellCount[thread].setSize(interiorRegionCount);
	    parallelInteriorRegionCellCount[thread].constant(0);
	}

	buildInteriorRegionCOM(parallelInteriorRegionCOM,
				parallelInteriorRegionCellCount,
				interiorRegionIndices,
				materialCellLabels);

	UT_Array<SolveReal> interiorRegionCellCount;
	interiorRegionCellCount.setSize(interiorRegionCount);
	interiorRegionCellCount.constant(0);

	tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
	{
	    for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
	    {
		for (int thread = 0; thread < threadCount; ++thread)
		{
		    interiorRegionCOM[interiorRegion] += parallelInteriorRegionCOM[thread][interiorRegion];
		    interiorRegionCellCount[interiorRegion] += parallelInteriorRegionCellCount[thread][interiorRegion];
		}
	    }
	});

	tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
	{
	    for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
		interiorRegionCOM[interiorRegion] *= dx / interiorRegionCellCount[interiorRegion];
	});    
    }

    ////////////////////////////////////////////
    //
    // Build best fit affine field for each interior region
    //
    ////////////////////////////////////////////

    UT_Array<UT_Vector3SR> interiorBestFitLinearVelocities;
    interiorBestFitLinearVelocities.setSize(interiorRegionCount);

    UT_Array<GradientMatrix> interiorBestFitVelocityGradients;
    interiorBestFitVelocityGradients.setSize(interiorRegionCount);

    {
	std::cout << "\n// Build best fit affine field for pre-solve velocity" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Build best fit affine field for pre-solve velocity");

	UT_Array<UT_Array<MassMatrix>> parallelBestFitMatrix;
	parallelBestFitMatrix.setSize(threadCount);

	UT_Array<UT_Array<ColumnVector>> parallelBestFitRHS;
	parallelBestFitRHS.setSize(threadCount);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    parallelBestFitMatrix[thread].setSize(interiorRegionCount);
	    parallelBestFitMatrix[thread].constant(MassMatrix::Zero());

	    parallelBestFitRHS[thread].setSize(interiorRegionCount);
	    parallelBestFitRHS[thread].constant(ColumnVector::Zero());
	}

	buildBestFitSystems(parallelBestFitMatrix,
			    parallelBestFitRHS,
			    materialCellLabels,
			    interiorRegionIndices,
			    *velocity,
			    interiorRegionCOM,
			    dx);

	// Compile best fit matrices and rhs vectors
	UT_Array<MassMatrix> interiorBestFitMatrix;
	interiorBestFitMatrix.setSize(interiorRegionCount);
	interiorBestFitMatrix.constant(MassMatrix::Zero());

	UT_Array<ColumnVector> interiorBestFitRHS;
	interiorBestFitRHS.setSize(interiorRegionCount);
	interiorBestFitRHS.constant(ColumnVector::Zero());

	tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
	{
	    for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
	    {
		for (int thread = 0; thread < threadCount; ++thread)
		{
		    interiorBestFitMatrix[interiorRegion] += parallelBestFitMatrix[thread][interiorRegion];
		    interiorBestFitRHS[interiorRegion] += parallelBestFitRHS[thread][interiorRegion];
		}
	    }
	});

	tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
	{
	    for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
	    {
		ColumnVector interiorBestFitSolution = interiorBestFitMatrix[interiorRegion].fullPivLu().solve(interiorBestFitRHS[interiorRegion]);

		// TODO: reduce to simple loops
		interiorBestFitLinearVelocities[interiorRegion][0] = interiorBestFitSolution(0);
		interiorBestFitLinearVelocities[interiorRegion][1] = interiorBestFitSolution(1);
		interiorBestFitLinearVelocities[interiorRegion][2] = interiorBestFitSolution(2);

		interiorBestFitVelocityGradients[interiorRegion](0, 0) = interiorBestFitSolution(3);
		interiorBestFitVelocityGradients[interiorRegion](0, 1) = interiorBestFitSolution(4);
		interiorBestFitVelocityGradients[interiorRegion](0, 2) = interiorBestFitSolution(5);

		interiorBestFitVelocityGradients[interiorRegion](1, 0) = interiorBestFitSolution(6);
		interiorBestFitVelocityGradients[interiorRegion](1, 1) = interiorBestFitSolution(7);
		interiorBestFitVelocityGradients[interiorRegion](1, 2) = interiorBestFitSolution(8);

		interiorBestFitVelocityGradients[interiorRegion](2, 0) = interiorBestFitSolution(9);
		interiorBestFitVelocityGradients[interiorRegion](2, 1) = interiorBestFitSolution(10);
		interiorBestFitVelocityGradients[interiorRegion](2, 2) = -(interiorBestFitSolution(3) + interiorBestFitSolution(7));
	    }
	});
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
		    interiorRegionIndices,
		    activeCellIndices,
		    *velocity,
		    cutCellWeights,
		    solidVelocity,
		    interiorBestFitLinearVelocities,
		    interiorBestFitVelocityGradients,
		    interiorRegionCOM,
		    dx);
    }

    ////////////////////////////////////////////
    //
    // Apply goal divergence for liquid cells
    //
    ////////////////////////////////////////////
    
    // TODO: replace with proper volume controls in bubbles - exterior and affine regions
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
	SIM_RawIndexField activeRegionIndices;
	exint activeRegionCount = 0;

	{
	    std::cout << "\n// Project null space" << std::endl;
	    UT_PerfMonAutoSolveEvent event(this, "Project null space");

	    // Build connected components for active regions
	    {
		std::cout << "\n    Build combined active regions indices" << std::endl;

		activeRegionIndices.match(liquidSurface);

		SIM_VolumetricConnectedComponentBuilder activeRegionBuilder(activeRegionIndices, materialCellLabels, cutCellWeights.data());
		activeRegionCount = activeRegionBuilder.buildConnectedComponents([](const exint label) { return label == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
														label == MaterialLabels::EXTERIOR_LIQUID_CELL; });
	    
		HDK::Utilities::overwriteIndices(activeRegionIndices,
						    SIM_VolumetricConnectedComponentBuilder::INACTIVE_REGION,
						    HDK::Utilities::UNLABELLED_CELL);

		std::cout << "\n    Active region count: " << activeRegionCount << std::endl;
	    }

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

	if (doTrackBubbleIDs)
	{
	    UT_Array<UT_Array<SolveReal>> parallelUntrackedDivergence;
	    parallelUntrackedDivergence.setSize(threadCount);

	    UT_Array<UT_Array<exint>> parallelTrackedCellCount;
	    parallelTrackedCellCount.setSize(threadCount);

	    UT_Array<UT_Array<exint>> parallelUntrackedCellCount;
	    parallelUntrackedCellCount.setSize(threadCount);

	    for (int thread = 0; thread < threadCount; ++thread)
	    {
		parallelUntrackedDivergence[thread].setSize(activeRegionCount);
		parallelUntrackedDivergence[thread].constant(0);

		parallelTrackedCellCount[thread].setSize(activeRegionCount);
		parallelTrackedCellCount[thread].constant(0);

		parallelUntrackedCellCount[thread].setSize(isBubbleTracked.size());
		parallelUntrackedCellCount[thread].constant(0);
	    }

	    setUntrackedDivergence(parallelUntrackedDivergence,
				    parallelTrackedCellCount,
				    parallelUntrackedCellCount,
				    rhsVector,
				    materialCellLabels,
				    activeCellIndices,
				    activeRegionIndices,
				    fullBubbleRegionIndices,
				    isBubbleTracked, dx, dt);


	    UT_Array<SolveReal> untrackedDivergence;
	    untrackedDivergence.setSize(activeRegionCount);
	    untrackedDivergence.constant(0);

	    UT_Array<exint> trackedCellCount;
	    trackedCellCount.setSize(activeRegionCount);
	    trackedCellCount.constant(0);

	    for (int thread = 0; thread < threadCount; ++thread)
		for (exint activeRegion = 0; activeRegion < activeRegionCount; ++activeRegion)
		{
		    untrackedDivergence[activeRegion] += parallelUntrackedDivergence[thread][activeRegion];
		    trackedCellCount[activeRegion] += parallelTrackedCellCount[thread][activeRegion];
		}

	    for (exint activeRegion = 0; activeRegion < activeRegionCount; ++activeRegion)
		untrackedDivergence[activeRegion] /= SolveReal(trackedCellCount[activeRegion]);

	    UT_Array<exint> untrackedCellCount;
	    untrackedCellCount.setSize(isBubbleTracked.size());
	    untrackedCellCount.constant(0);

	    for (int thread = 0; thread < threadCount; ++thread)
		for (exint bubbleRegion = 0; bubbleRegion < isBubbleTracked.size(); ++bubbleRegion)
		{
		    if (!isBubbleTracked[bubbleRegion])
			untrackedCellCount[bubbleRegion] += parallelUntrackedCellCount[thread][bubbleRegion];
		}

	    for (exint bubbleRegion = 0; bubbleRegion < isBubbleTracked.size(); ++bubbleRegion)
	    {
		if (!isBubbleTracked[bubbleRegion])
		    std::cout << "    Untracked bubble region " << bubbleRegion << " of size " << untrackedCellCount[bubbleRegion] << std::endl;
		else assert(untrackedCellCount[bubbleRegion] == 0);
	    }

	    // Remove divergence in tracked cells
	    removeUntrackedDivergence(rhsVector,
					materialCellLabels,
					activeCellIndices,
					activeRegionIndices,
					fullBubbleRegionIndices,
					isBubbleTracked,
					untrackedDivergence);

	    std::cout << "    Average divergence in rhs vector after untracked bubble removal: " << rhsVector.sum() / SolveReal(rhsVector.rows()) << std::endl;
	}
    }

    ////////////////////////////////////////////
    //
    // Build poisson sparse elements for active cells
    //
    ////////////////////////////////////////////

    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> poissonMatrix(activeCellCount, activeCellCount);
    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> diagonalPrecondMatrix(activeCellCount, activeCellCount);

    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> affineCollectionMatrix(11 * interiorRegionCount, activeCellCount);

    {
	std::cout << "\n// Build Poisson system" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Build Poisson system");

	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelPoissonElements(threadCount);
	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelDiagonalPrecondElements(threadCount);
	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelAffineBubbleCollectionElements(threadCount);
	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelAffineLiquidCollectionElements(threadCount);

	buildLinearSystem(parallelPoissonElements,
			    parallelDiagonalPrecondElements,
			    parallelAffineBubbleCollectionElements,
			    parallelAffineLiquidCollectionElements,
			    materialCellLabels,
			    interiorRegionIndices,
			    activeCellIndices,
			    cutCellWeights,
			    liquidSurface,
			    liquidDensity,
			    bubbleDensity,
			    interiorRegionCOM,
			    dx);

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

	std::vector<Eigen::Triplet<SolveReal>> affineCollectionElements;

	// Compile affine elements
	{
	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelAffineBubbleCollectionElements[thread].size();

	    affineCollectionElements.reserve(listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
		affineCollectionElements.insert(affineCollectionElements.end(), parallelAffineBubbleCollectionElements[thread].begin(), parallelAffineBubbleCollectionElements[thread].end());
	}

	{
	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelAffineLiquidCollectionElements[thread].size();

	    affineCollectionElements.reserve(affineCollectionElements.size() + listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
		affineCollectionElements.insert(affineCollectionElements.end(), parallelAffineLiquidCollectionElements[thread].begin(), parallelAffineLiquidCollectionElements[thread].end());
	}

	affineCollectionMatrix.setFromTriplets(affineCollectionElements.begin(), affineCollectionElements.end());
	affineCollectionMatrix.makeCompressed();
    }

    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> affineMassMatrix(11 * interiorRegionCount, 11 * interiorRegionCount);
    
    UT_Array<MassMatrix> invertedAffineMassMatrices;
    {
	std::cout << "\n// Build affine mass matrix" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Build affine mass matrix");

	UT_Array<UT_Array<MassMatrix>> parallelAffineMassMatrix;
	parallelAffineMassMatrix.setSize(threadCount);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    parallelAffineMassMatrix[thread].setSize(interiorRegionCount);
	    parallelAffineMassMatrix[thread].constant(MassMatrix::Zero());
	}

	buildAffineMassMatrixSystems(parallelAffineMassMatrix,
					materialCellLabels,
					interiorRegionIndices,
					interiorRegionCOM,
					dx);

	// Compile affine mass matrix elements
        UT_Array<MassMatrix> interiorRegionMassMatrices;
	interiorRegionMassMatrices.setSize(interiorRegionCount);
	interiorRegionMassMatrices.constant(MassMatrix::Zero());

	tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
	{
	    for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
	    {
		for (int thread = 0; thread < threadCount; ++thread)
		    interiorRegionMassMatrices[interiorRegion] += parallelAffineMassMatrix[thread][interiorRegion];
	    }
	});

	std::vector<Eigen::Triplet<SolveReal>> affineMassMatrixElements;
	affineMassMatrixElements.reserve(11 * 11 * interiorRegionCount);

	invertedAffineMassMatrices.setSize(interiorRegionCount);

	// Build combined mass matrix system
	{
	    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<SolveReal>>> parallelMassMatrixElements;
    	    tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
	    {
		auto &localMassMatrixElements = parallelMassMatrixElements.local();
		for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
		{
		    // Invert mass matrix
		    invertedAffineMassMatrices[interiorRegion] = interiorRegionMassMatrices[interiorRegion].inverse();

		    for (int row = 0; row < 11; ++row)
			for (int col = 0; col < 11; ++col)
			{
			    if (invertedAffineMassMatrices[interiorRegion](row, col) != 0)
				localMassMatrixElements.emplace_back(11 * interiorRegion + row, 11 * interiorRegion + col, invertedAffineMassMatrices[interiorRegion](row, col));
			}
		}
	    });

	    exint listSize = 0;

	    parallelMassMatrixElements.combine_each([&](const std::vector<Eigen::Triplet<SolveReal>> &localMassMatrixElements)
	    {
		listSize += localMassMatrixElements.size();
	    });

	    affineMassMatrixElements.reserve(listSize);

	    parallelMassMatrixElements.combine_each([&](const std::vector<Eigen::Triplet<SolveReal>> &localMassMatrixElements)
	    {
		affineMassMatrixElements.insert(affineMassMatrixElements.end(), localMassMatrixElements.begin(), localMassMatrixElements.end());
	    });
	}

	affineMassMatrix.setFromTriplets(affineMassMatrixElements.begin(), affineMassMatrixElements.end());
	affineMassMatrix.makeCompressed();
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
	// Build list of interior boundary cells to distribute
	// affine vectors with
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

	UT_Array<UT_Vector3I> interiorBoundaryCells;
	{
	    std::cout << "\n// Build interior region boundary cell list" << std::endl;
	    UT_PerfMonAutoSolveEvent event(this, "Build interior region boundary cell list");

	    UT_Array<UT_Array<UT_Vector3I>> parallelInteriorBoundaryCells;
	    parallelInteriorBoundaryCells.setSize(threadCount);

	    buildInteriorBoundaryCells(parallelInteriorBoundaryCells,
					materialCellLabels);

	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelInteriorBoundaryCells[thread].size();

	    interiorBoundaryCells.bumpCapacity(listSize);
	
	    for (int thread = 0; thread < threadCount; ++thread)
		interiorBoundaryCells.concat(parallelInteriorBoundaryCells[thread]);

	    UTparallelSort(interiorBoundaryCells.begin(), interiorBoundaryCells.end(), cellCompare);
	}

	std::cout << "// Solve affine tiled linear system" << std::endl;
    	UT_PerfMonAutoSolveEvent event(this, "Solve affine tiled linear system");

	Vector tempVector0 = Vector::Zero(11 * interiorRegionCount);
	Vector tempVector1 = Vector::Zero(11 * interiorRegionCount);

	auto MatrixVectorMultiply = [&](Vector &destinationVector, const Vector &sourceVector)
	{
	    // UT_StopWatch multiplyTimer;
	    // multiplyTimer.start();

	    destinationVector.noalias() = poissonMatrix * sourceVector;

	    // time = multiplyTimer.stop();
	    // std::cout << "\n    Poisson matrix-vector multiply time: " << time << std::endl;
	    // multiplyTimer.clear();
	    // multiplyTimer.start();

	    if (interiorRegionCount > 0)
	    {
		tempVector0.noalias() = affineCollectionMatrix  * sourceVector;

		// time = multiplyTimer.stop();
		// std::cout << "\n    Collect affine matrix-vector multiply time: " << time << std::endl;
		// multiplyTimer.clear();
		// multiplyTimer.start();

		tempVector1.noalias() = affineMassMatrix * tempVector0;

		// time = multiplyTimer.stop();
		// std::cout << "\n    Apply mass matrix matrix-vector multiply time: " << time << std::endl;
		// multiplyTimer.clear();
		// multiplyTimer.start();

		distributeAffineVectors(destinationVector,
					tempVector1,
					interiorBoundaryCells,
					materialCellLabels,
					interiorRegionIndices,
					activeCellIndices,
					interiorRegionCOM,
					dx);

		// time = multiplyTimer.stop();
		// std::cout << "    Distribute affine matrix-vector multiply time: " << time << std::endl;
	    }
	};

	if (getUseMGPreconditioner() && !doApplyAffineToLiquid)
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

		UT_StopWatch mgBuilderTimer;
		mgBuilderTimer.start();

		UT_VoxelArray<int> baseDomainCellLabels;
		baseDomainCellLabels.size(activeCellIndices.getVoxelRes()[0],
					    activeCellIndices.getVoxelRes()[1],
					    activeCellIndices.getVoxelRes()[2]);

		baseDomainCellLabels.constant(MGCellLabels::EXTERIOR_CELL);

		std::cout << "\n    Build Liquid Multigrid Domain Labels" << std::endl;

		buildMGDomainLabels(baseDomainCellLabels,
				    materialCellLabels,
				    activeCellIndices);

		std::cout << "\n    Build Boundary Weights" << std::endl;

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
					    interiorRegionIndices,
					    baseDomainCellLabels,
					    liquidSurface,
					    *cutCellWeights[axis],
					    liquidDensity,
					    bubbleDensity,
					    axis);
		}
		
		auto time = mgBuilderTimer.stop();
		std::cout << "    Build base MG domain time: " << time << std::endl;
		mgBuilderTimer.clear();
		mgBuilderTimer.start();

		std::cout << "\n    Build Expanded Multigrid Domain Labels" << std::endl;
		
		// Build expanded domain
		auto isExteriorCell = [](const int value) { return value == MGCellLabels::EXTERIOR_CELL; };
		auto isInteriorCell = [](const int value) { return value == MGCellLabels::INTERIOR_CELL; };
		auto isDirichletCell = [](const int value) { return value == MGCellLabels::DIRICHLET_CELL; };

		std::pair<UT_Vector3I, int> mgSettings = HDK::GeometricMultigridOperators::buildExpandedCellLabels(mgDomainCellLabels, baseDomainCellLabels, isExteriorCell, isInteriorCell, isDirichletCell);

		mgExpandedOffset = mgSettings.first;
		mgLevels = mgSettings.second;
		
		std::cout << "\n    Build Expanded Boundary Weights" << std::endl;

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

		time = mgBuilderTimer.stop();
		std::cout << "    Build expaned MG domain time: " << time << std::endl;
		mgBuilderTimer.clear();
		mgBuilderTimer.start();
	    }

	    UT_StopWatch mgSolverConstructionTimer;
	    std::cout << "\n    Build MG Solver" << std::endl;
	    assert(unitTestMGLabels(mgDomainCellLabels, materialCellLabels, activeCellIndices, interiorRegionIndices, mgExpandedOffset));
	    // TODO: include MG build time in profiler
	    HDK::GeometricMultigridPoissonSolver mgPreconditioner(mgDomainCellLabels,
								    mgBoundaryWeights,
								    mgLevels,
								    true /* use Gauss-Seidel smoother */);

	    time = mgSolverConstructionTimer.stop();
	    std::cout << "    Construct MG time: " << time << std::endl;

	    //
	    // Initialize bubble smoother
	    //

	    // Build band of liquid cells at the bubble boundary

	    UT_Array<UT_Vector3I> bubbleSmootherCells;
	    UT_Array<UT_Vector3I> liquidCopyCells;
	    
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

		    exint liquidSmootherBand = 1;
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

			std::cout << "    Liquid smoother layer size: " << newLiquidCells.size() << std::endl;

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

		    UT_Array<UT_Array<UT_Vector3I>> parallelBubbleSmootherCells;
		    parallelBubbleSmootherCells.setSize(threadCount);

		    UT_Array<UT_Array<UT_Vector3I>> parallelLiquidCopyCells;
		    parallelLiquidCopyCells.setSize(threadCount);

		    // Set smoother cells
		    buildBubbleSmootherCells(parallelBubbleSmootherCells,
						parallelLiquidCopyCells,
						materialCellLabels,
						visitedLiquidCells,
						cutCellWeights);

		    exint smootherListSize = 0;
		    exint copyListSize = 0;
		    
		    for (int thread = 0; thread < threadCount; ++thread)
		    {
			smootherListSize += parallelBubbleSmootherCells[thread].size();
			copyListSize += parallelLiquidCopyCells[thread].size();
		    }
		    
		    bubbleSmootherCells.bumpCapacity(smootherListSize);
		    liquidCopyCells.bumpCapacity(copyListSize);

		    for (int thread = 0; thread < threadCount; ++thread)
		    {
			bubbleSmootherCells.concat(parallelBubbleSmootherCells[thread]);
			liquidCopyCells.concat(parallelLiquidCopyCells[thread]);
		    }

		    UTparallelSort(bubbleSmootherCells.begin(), bubbleSmootherCells.end(), cellCompare);
		    UTparallelSort(liquidCopyCells.begin(), liquidCopyCells.end(), cellCompare);
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

	    UT_Array<UT_Array<ColumnVector>> parallelAffineVectors;
	    parallelAffineVectors.setSize(threadCount);

	    for (int thread = 0; thread < threadCount; ++thread)
		parallelAffineVectors[thread].setSize(interiorRegionCount);

	    UT_Array<ColumnVector> affineVectors;
	    affineVectors.setSize(interiorRegionCount);

	    UT_Array<ColumnVector> invertedAffineVectors;
	    invertedAffineVectors.setSize(interiorRegionCount);

	    UT_Array<SolveReal> tempSmootherDestinationValues;
	    tempSmootherDestinationValues.setSize(bubbleSmootherCells.size());

	    // Build tile list
	    UT_Array<bool> isSmootherTileOccupied;
	    isSmootherTileOccupied.setSize(smootherDestinationGrid.numTiles());
	    isSmootherTileOccupied.constant(false);
	    
	    HDK::Utilities::findOccupiedIndexTiles(isSmootherTileOccupied,
						    bubbleSmootherCells,
						    materialCellLabels);

	    HDK::Utilities::findOccupiedIndexTiles(isSmootherTileOccupied,
						    liquidCopyCells,
						    materialCellLabels);

	    auto MultigridPreconditioner = [&](Vector &destinationVector, const Vector &sourceVector)
	    {
		int bubbleSmootherIterations = 20;
		int preconditionerIterations = 1;

		//
		// Apply bubble smoother before applying v-cycle
		//

		UT_StopWatch timer;
		timer.clear();
		timer.start();

		// Reset and uncompress smoother grids
		smootherSourceGrid.constant(0);
		smootherDestinationGrid.constant(0);

		HDK::GeometricMultigridOperators::uncompressTiles(smootherSourceGrid, isSmootherTileOccupied);
		HDK::GeometricMultigridOperators::uncompressTiles(smootherDestinationGrid, isSmootherTileOccupied);

		// Copy source vector to smoother grid
		copySourceToSmootherGrid(smootherSourceGrid,
					    sourceVector,
					    activeCellIndices,
					    bubbleSmootherCells);

		time = timer.stop();
		std::cout << "      Copy source to smoother grid. Compute time: " << time << std::endl;
		timer.clear();
		timer.start();

		for (int smootherIteration = 0; smootherIteration < bubbleSmootherIterations; ++smootherIteration)
		{
		    affineVectors.constant(ColumnVector::Zero());

		    buildAffineVectors(affineVectors,
					materialCellLabels,
					interiorRegionIndices,
					smootherDestinationGrid,
					interiorBoundaryCells,
					interiorRegionCOM,
					interiorRegionCount,
					dx);

		    // Scale by inverse
		    invertedAffineVectors.constant(ColumnVector::Zero());

		    tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
		    {
			for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
			{
			    for (int row = 0; row < 11; ++row)
			    {
				for (int col = 0; col < 11; ++col)
				    invertedAffineVectors[interiorRegion](row) += invertedAffineMassMatrices[interiorRegion](row, col) * affineVectors[interiorRegion](col);

				invertedAffineVectors[interiorRegion](row) /= bubbleDensity;
			    }
			}
		    });

		    tempSmootherDestinationValues.constant(0);

		    applyBubbleSmoother(tempSmootherDestinationValues,
					materialCellLabels,
					interiorRegionIndices,
					liquidSurface,
					cutCellWeights,
					mgDomainCellLabels,
					smootherDestinationGrid,
					smootherSourceGrid,
					invertedAffineVectors,
					interiorRegionCOM,
					bubbleSmootherCells,
					mgExpandedOffset,
					liquidDensity,
					bubbleDensity,
					dx);

		    copyDestinationVectorToGrid(smootherDestinationGrid,
						tempSmootherDestinationValues,
						bubbleSmootherCells);
		}

		time = timer.stop();
		std::cout << "      Jacobi smoothing for exterior and boudary liquid cells. Compute time: " << time << std::endl;
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

		    UTparallelForLightItems(UT_BlockedRange<exint>(0, bubbleSmootherCells.size()), [&](const UT_BlockedRange<exint> &range)
		    {
			if (boss->opInterrupt())
			    return;

			for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
			{
			    UT_Vector3I cell = bubbleSmootherCells[cellIndex];
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
		std::cout << "      Uncompress MG grids. Compute time: " << time << std::endl;
		timer.clear();
		timer.start();

		copySourceToMGGrid(mgSourceGrid,
				    materialCellLabels,
				    activeCellIndices,
				    mgDomainCellLabels,
				    sourceVector,
				    mgExpandedOffset,
				    liquidDensity);

		time = timer.stop();
		std::cout << "      Copy source to MG grid. Compute time: " << time << std::endl;
		timer.clear();
		timer.start();

		for (int precondIter = 0; precondIter < preconditionerIterations; ++precondIter)
		{
		    //
		    // Apply multigrid v-cycle
		    //

		    applyDirichletToMGGrid(mgSourceGrid,
					    mgDestinationGrid,
					    materialCellLabels,
					    activeCellIndices,
					    liquidSurface,
					    cutCellWeights,
					    mgDomainCellLabels,
					    smootherDestinationGrid,
					    sourceVector,
					    bubbleSmootherCells,
					    mgExpandedOffset,
					    liquidDensity,
					    bubbleDensity);

		    time = timer.stop();
		    std::cout << "      Copy to MG and apply Dirichlet boundaries. Compute time: " << time << std::endl;
		    timer.clear();
		    timer.start();

		    // Apply multigrid v-cycle using Dirichlet conditions from the box smoother
		    mgPreconditioner.applyVCycle(mgDestinationGrid, mgSourceGrid, true);

		    time = timer.stop();
		    std::cout << "      Apply v-cycle time: " << time << std::endl;
		    timer.clear();
		    timer.start();

		    copyMGToSmoother(smootherDestinationGrid,
					materialCellLabels,
					mgDomainCellLabels,
					mgDestinationGrid,
					liquidCopyCells,
					mgExpandedOffset);

		    time = timer.stop();
		    std::cout << "      Copy MG grid to smoother grid. Compute time: " << time << std::endl;
		    timer.clear();
		    timer.start();

		    for (int smootherIteration = 0; smootherIteration < bubbleSmootherIterations; ++smootherIteration)
		    {
			UT_StopWatch smootherTimer;
			smootherTimer.start();

			affineVectors.constant(ColumnVector::Zero());

			buildAffineVectors(affineVectors,
					    materialCellLabels,
					    interiorRegionIndices,
					    smootherDestinationGrid,
					    interiorBoundaryCells,
					    interiorRegionCOM,
					    interiorRegionCount,
					    dx);

			time = smootherTimer.stop();
			std::cout << "	Build affine vectors. Compute time: " << time << std::endl;
			smootherTimer.clear();
			smootherTimer.start();

			// Scale by inverse
			invertedAffineVectors.constant(ColumnVector::Zero());

    			tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
			{
			    for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
			    {
				for (int row = 0; row < 11; ++row)
				{
				    for (int col = 0; col < 11; ++col)
					invertedAffineVectors[interiorRegion](row) += invertedAffineMassMatrices[interiorRegion](row, col) * affineVectors[interiorRegion](col);

				    invertedAffineVectors[interiorRegion](row) /= bubbleDensity;
				}
			    }
			});

			time = smootherTimer.stop();
			std::cout << "	Scale by inverted mass matrix. Compute time: " << time << std::endl;
			smootherTimer.clear();
			smootherTimer.start();

			tempSmootherDestinationValues.constant(0);

			applyBubbleSmoother(tempSmootherDestinationValues,
					    materialCellLabels,
					    interiorRegionIndices,
					    liquidSurface,
					    cutCellWeights,
					    mgDomainCellLabels,
					    smootherDestinationGrid,
					    smootherSourceGrid,
					    invertedAffineVectors,
					    interiorRegionCOM,
					    bubbleSmootherCells,
					    mgExpandedOffset,
					    liquidDensity,
					    bubbleDensity,
					    dx);

			time = smootherTimer.stop();
			std::cout << "	Apply bubble smoother. Compute time: " << time << std::endl;
			smootherTimer.clear();
			smootherTimer.start();

			copyDestinationVectorToGrid(smootherDestinationGrid,
						    tempSmootherDestinationValues,
						    bubbleSmootherCells);

			time = smootherTimer.stop();
			std::cout << "	Copy temp values to smoother grid. Compute time: " << time << std::endl;
			smootherTimer.clear();
			smootherTimer.start();
		    }

		    time = timer.stop();
		    std::cout << "      Jacobi smoothing for exterior and boudary liquid cells. Compute time: " << time << std::endl;
		    timer.clear();
		    timer.start();
		}

		// Transfer results back to main vector
		copyMGToDestinationVector(destinationVector,
					    materialCellLabels,
					    activeCellIndices,
					    mgDomainCellLabels,
					    mgDestinationGrid,
					    mgExpandedOffset);

		time = timer.stop();
		std::cout << "      Copy MG grid to output vector. Compute time: " << time << std::endl;
		timer.clear();
		timer.start();

		copySmootherToDestinationVector(destinationVector,
						activeCellIndices,
						smootherDestinationGrid,
						bubbleSmootherCells);

		time = timer.stop();
		std::cout << "      Copy smoother grid to output time: " << time << std::endl;
		timer.clear();
		timer.start();
	    };

	    HDK::Utilities::solveConjugateGradient(solutionVector,
						    rhsVector,
						    MatrixVectorMultiply,
						    MultigridPreconditioner,
						    getTolerance(),
						    getMaxIterations());
	}
	else
	{
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
	}

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
    // Collect affine interior velocity solution
    //
    ////////////////////////////////////////////
    
    UT_Array<UT_Vector3SR> interiorLinearVelocities;
    interiorLinearVelocities.setSize(interiorRegionCount);
    interiorLinearVelocities.constant(UT_Vector3SR(0,0,0));

    UT_Array<GradientMatrix> interiorVelocityGradients(interiorRegionCount);
    interiorVelocityGradients.setSize(interiorRegionCount);
    interiorVelocityGradients.constant(GradientMatrix::Zero());

    {
	Vector affineUpdate = affineMassMatrix * affineCollectionMatrix  * solutionVector;
	tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
	{
	    for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
	    {
		interiorLinearVelocities[interiorRegion] = interiorBestFitLinearVelocities[interiorRegion];

		for (int axis : {0,1,2})
		    interiorLinearVelocities[interiorRegion][axis] += affineUpdate(11 * interiorRegion + axis);

		interiorVelocityGradients[interiorRegion] = interiorBestFitVelocityGradients[interiorRegion];

		interiorVelocityGradients[interiorRegion](0, 0) += affineUpdate(11 * interiorRegion + 3);
		interiorVelocityGradients[interiorRegion](0, 1) += affineUpdate(11 * interiorRegion + 4);
		interiorVelocityGradients[interiorRegion](0, 2) += affineUpdate(11 * interiorRegion + 5);

		interiorVelocityGradients[interiorRegion](1, 0) += affineUpdate(11 * interiorRegion + 6);
		interiorVelocityGradients[interiorRegion](1, 1) += affineUpdate(11 * interiorRegion + 7);
		interiorVelocityGradients[interiorRegion](1, 2) += affineUpdate(11 * interiorRegion + 8);

		interiorVelocityGradients[interiorRegion](2, 0) += affineUpdate(11 * interiorRegion + 9);
		interiorVelocityGradients[interiorRegion](2, 1) += affineUpdate(11 * interiorRegion + 10);
		interiorVelocityGradients[interiorRegion](2, 2) += -(affineUpdate(11 * interiorRegion + 3) + affineUpdate(11 * interiorRegion + 7));
	    }
	});
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
				    interiorRegionIndices,
				    activeCellIndices,
				    *pressure,
				    liquidSurface,
				    *cutCellWeights[axis],
				    interiorRegionCOM,
				    interiorLinearVelocities,
				    interiorVelocityGradients,
				    constantLiquidDensity,
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

	UT_Array<SolveReal> parallelAccumulatedDivergence;
	parallelAccumulatedDivergence.setSize(threadCount);
	parallelAccumulatedDivergence.constant(0);

	UT_Array<SolveReal> parallelMaxDivergence;
	parallelMaxDivergence.setSize(threadCount);
	parallelMaxDivergence.constant(0);

	UT_Array<SolveReal> parallelCellCount;
	parallelCellCount.setSize(threadCount);
	parallelCellCount.constant(0);

	computeResultingDivergence(parallelAccumulatedDivergence,
				    parallelMaxDivergence,
				    parallelCellCount,
				    materialCellLabels,
				    *velocity,
				    cutCellWeights,
				    solidVelocity);

	SolveReal accumulatedDivergence = 0;
	SolveReal maxDivergence = 0;
	SolveReal cellCount = 0;

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    accumulatedDivergence += parallelAccumulatedDivergence[thread];
	    maxDivergence = std::max(maxDivergence, parallelMaxDivergence[thread]);
	    cellCount += parallelCellCount[thread];
	}

	std::cout << "    Max divergence: " << maxDivergence << std::endl;
	std::cout << "    Accumulated divergence: " << accumulatedDivergence << std::endl;
	std::cout << "    Average divergence: " << accumulatedDivergence / cellCount << std::endl;
    }

    pressureField->pubHandleModification();
    velocity->pubHandleModification();
    validFaces->pubHandleModification();

    return true;
}

void
HDK_AffineBubblePressureSolver::setMaterialCellLabels(SIM_RawIndexField &materialCellLabels) const
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
		    cellLabel = MaterialLabels::INTERIOR_LIQUID_CELL;
		else if (label == HDK::Utilities::AIR_CELL)
		    cellLabel = MaterialLabels::INTERIOR_BUBBLE_CELL;
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
			vit.setValue(MaterialLabels::INTERIOR_LIQUID_CELL);
		    else if (label == HDK::Utilities::AIR_CELL)
			vit.setValue(MaterialLabels::INTERIOR_BUBBLE_CELL);
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
HDK_AffineBubblePressureSolver::buildInitialFluidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExteriorCellList,
								const SIM_RawIndexField &materialCellLabels,
								const std::array<const SIM_RawField *, 3> &cutCellWeights,
								const exint interiorMaterialLabel,
								const exint opposingInteriorMaterialLabel) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInitialExteriorLayerAlgorithm;
    buildInitialExteriorLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<UT_Vector3I> &localExteriorCellList = parallelExteriorCellList[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() || vit.getValue() == interiorMaterialLabel)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == interiorMaterialLabel)
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

				if (getFieldValue(*cutCellWeights[axis], face) > 0 && getFieldValue(materialCellLabels, adjacentCell) != interiorMaterialLabel)
				{
				    // Debug check. Remove once verified
				    assert(getFieldValue(materialCellLabels, adjacentCell) == opposingInteriorMaterialLabel);
				    isBoundaryCell = true;
				}
			    }

			if (isBoundaryCell)
			    localExteriorCellList.append(cell);
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::buildNextFluidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExteriorCellList,
							    const UT_Array<UT_Vector3I> &oldExteriorCellList,
							    const SIM_RawIndexField &materialCellLabels,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights,
							    const exint exteriorMaterialLabel,
							    const exint interiorMaterialLabel) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildNextExteriorBandAlgorithm;
    buildNextExteriorBandAlgorithm.run([&](const UT_JobInfo &info)
    {
	if (boss->opInterrupt())
	    return 0;

	UT_Array<UT_Vector3I> &localExteriorCellList = parallelExteriorCellList[info.job()];

	exint startIndex, endIndex;
	const exint elementSize = oldExteriorCellList.entries();
	info.divideWork(elementSize, startIndex, endIndex);

	if (startIndex > 0)
	{
	    while (startIndex != endIndex && oldExteriorCellList[startIndex] == oldExteriorCellList[startIndex - 1])
		++startIndex;
	}
    
	UT_Vector3I oldCell(-1,-1,-1);

	const exint localEndIndex = endIndex;
	for (exint cellIndex = startIndex; cellIndex < localEndIndex; ++cellIndex)
	{
	    UT_Vector3I cell = oldExteriorCellList[cellIndex];

	    if (cell == oldCell)
		continue;

	    oldCell = cell;

	    assert(getFieldValue(materialCellLabels, cell) == exteriorMaterialLabel);

	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
	    
		    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
			continue;

		    UT_Vector3I face = cellToFaceMap(cell, axis, direction);
		    if (getFieldValue(*cutCellWeights[axis], face) > 0 &&
			getFieldValue(materialCellLabels, adjacentCell) == interiorMaterialLabel)
			localExteriorCellList.append(adjacentCell);
		}
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::buildInitialSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExteriorCellList,
								const SIM_RawIndexField &materialCellLabels,
								const std::array<const SIM_RawField *, 3> &cutCellWeights,
								const exint interiorMaterialLabel,
								const exint exteriorMaterialLabel) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInitialExteriorLayerAlgorithm;
    buildInitialExteriorLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<UT_Vector3I> &localExteriorCellList = parallelExteriorCellList[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() ||
		vit.getValue() == interiorMaterialLabel ||
		vit.getValue() == exteriorMaterialLabel)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == interiorMaterialLabel ||
			vitt.getValue() == exteriorMaterialLabel)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			bool isBoundaryCell = false;
			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
				{
				    isBoundaryCell = true;
				    continue;
				}

				UT_Vector3I face = cellToFaceMap(cell, axis, direction);

				if (getFieldValue(*cutCellWeights[axis], face) < 1)
				    isBoundaryCell = true;
			    }

			if (isBoundaryCell)
			    localExteriorCellList.append(cell);
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::setSolidLayerCells(SIM_RawIndexField &materialCellLabels,
						    SIM_RawIndexField &visitedCells,
						    const UT_Array<UT_Vector3I> &newExteriorCellList,
						    const exint interiorMaterialLabel,
						    const exint exteriorMaterialLabel) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    // Assign active cells
    UTparallelForLightItems(UT_BlockedRange<exint>(0, newExteriorCellList.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;

	exint startIndex = range.begin();
    
	if (startIndex > 0)
	{
	    while (startIndex != range.end() && newExteriorCellList[startIndex] == newExteriorCellList[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	for (exint cellIndex = startIndex; cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = newExteriorCellList[cellIndex];

	    if (cell == oldCell)
		continue;
	
	    oldCell = cell;

	    assert(getFieldValue(visitedCells, cell) == UNVISITED_CELL);
	    assert(getFieldValue(materialCellLabels, cell) == interiorMaterialLabel ||
		    getFieldValue(materialCellLabels, cell) == exteriorMaterialLabel);

	    setFieldValue(visitedCells, cell, VISITED_CELL);
	    setFieldValue(materialCellLabels, cell, exteriorMaterialLabel);
	}
    });
}

void
HDK_AffineBubblePressureSolver::buildNextSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExteriorCellList,
							    const UT_Array<UT_Vector3I> &oldExteriorCellList,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &visitedCells,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights,
							    const exint interiorMaterialLabel,
							    const exint exteriorMaterialLabel) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildNextSolidBoundaryLayerAlgorithm;
    buildNextSolidBoundaryLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	if (boss->opInterrupt())
	    return 0;

	UT_Array<UT_Vector3I> &localExteriorCellList = parallelExteriorCellList[info.job()];

	exint startIndex, endIndex;
	const exint elementSize = oldExteriorCellList.entries();
	info.divideWork(elementSize, startIndex, endIndex);

	if (startIndex > 0)
	{
	    while (startIndex != endIndex && oldExteriorCellList[startIndex] == oldExteriorCellList[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	const exint localEndIndex = endIndex;
	for (exint cellIndex = startIndex; cellIndex < localEndIndex; ++cellIndex)
	{
	    UT_Vector3I cell = oldExteriorCellList[cellIndex];

	    if (cell == oldCell)
		continue;

	    oldCell = cell;

	    assert(getFieldValue(materialCellLabels, cell) == exteriorMaterialLabel);
	    assert(getFieldValue(visitedCells, cell) == VISITED_CELL);

	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

		    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
			continue;

		    UT_Vector3I face = cellToFaceMap(cell, axis, direction);

		    if (getFieldValue(*cutCellWeights[axis], face) > 0)
		    {
			if (getFieldValue(visitedCells, adjacentCell) == UNVISITED_CELL &&
			    (getFieldValue(materialCellLabels, adjacentCell) == interiorMaterialLabel ||
			    getFieldValue(materialCellLabels, adjacentCell) == exteriorMaterialLabel))
			    localExteriorCellList.append(adjacentCell);
		    }
		}
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::buildExteriorCells(SIM_RawIndexField &materialCellLabels,
						    const std::array<const SIM_RawField *, 3> &cutCellWeights,
						    const exint exteriorMaterialLabel,
						    const exint interiorMaterialLabel,
						    const exint opposingInteriorMaterialLabel) const
{
    const int threadCount = UT_Thread::getNumProcessors();

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

    UT_Array<UT_Array<UT_Vector3I>> parallelExteriorBandList;
    parallelExteriorBandList.setSize(threadCount);

    UT_Array<UT_Vector3I> newExteriorCellsList;

    UT_Array<bool> isTileOccupiedList;
    isTileOccupiedList.setSize(materialCellLabels.field()->numTiles());

    {
	UT_PerfMonAutoSolveEvent event(this, "Build active cells along fluid-fluid boundary");
	std::cout << "\n// Build active cells along fluid-fluid boundary" << std::endl;

	// Build initial list of exterior cells
	buildInitialFluidBoundaryLayer(parallelExteriorBandList, materialCellLabels, cutCellWeights, interiorMaterialLabel, opposingInteriorMaterialLabel);

	const int exteriorBandwidth = getExteriorBandwidth();
	for (int layer = 0; layer < exteriorBandwidth; ++layer)
	{
	    // Combine parallel exterior band cells lists
	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelExteriorBandList[thread].size();

	    newExteriorCellsList.clear();
	    newExteriorCellsList.bumpCapacity(listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
	    {
		newExteriorCellsList.concat(parallelExteriorBandList[thread]);
		parallelExteriorBandList[thread].clear();
	    }

	    UTparallelSort(newExteriorCellsList.begin(), newExteriorCellsList.end(), cellCompare);

	    std::cout << "Exterior layer: " << layer << ". Exterior cell count: " << newExteriorCellsList.entries() << std::endl;

	    isTileOccupiedList.constant(false);
	    HDK::Utilities::findOccupiedIndexTiles(isTileOccupiedList,
						    newExteriorCellsList,
						    materialCellLabels);

	    HDK::Utilities::uncompressTiles(materialCellLabels, isTileOccupiedList);

	    UTparallelForLightItems(UT_BlockedRange<exint>(0, newExteriorCellsList.size()), [&](const UT_BlockedRange<exint> &range)
	    {
		using SIM::FieldUtils::setFieldValue;
		for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
		{
		    UT_Vector3I cell = newExteriorCellsList[cellIndex];
		    setFieldValue(materialCellLabels, cell, exteriorMaterialLabel);
		}
	    });

	    if (layer < exteriorBandwidth - 1)
	    {
		buildNextFluidBoundaryLayer(parallelExteriorBandList,
					    newExteriorCellsList,
					    materialCellLabels,
					    cutCellWeights,
					    exteriorMaterialLabel,
					    interiorMaterialLabel);
	    }
	}
    }
    {
	UT_PerfMonAutoSolveEvent event(this, "Build active cells along fluid-solid boundary");
    	std::cout << "\n// Build active cells along fluid-solid boundary" << std::endl;

	int solidBoundaryLayerSize = 2;

	for (int thread = 0; thread < threadCount; ++thread)
	    assert(parallelExteriorBandList[thread].isEmpty());

	buildInitialSolidBoundaryLayer(parallelExteriorBandList,
					materialCellLabels,
					cutCellWeights,
					interiorMaterialLabel,
					exteriorMaterialLabel);

    	SIM_RawIndexField visitedCells;
	visitedCells.match(materialCellLabels);
	visitedCells.makeConstant(UNVISITED_CELL);

	for (int layer = 0; layer < solidBoundaryLayerSize; ++layer)
	{
	    // Combine parallel active cell lists
	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelExteriorBandList[thread].size();

	    newExteriorCellsList.clear();
	    newExteriorCellsList.bumpCapacity(listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
	    {
		newExteriorCellsList.concat(parallelExteriorBandList[thread]);
		parallelExteriorBandList[thread].clear();
	    }

	    UTparallelSort(newExteriorCellsList.begin(), newExteriorCellsList.end(), cellCompare);

	    // Build tile list
	    isTileOccupiedList.constant(false);
	    HDK::Utilities::findOccupiedIndexTiles(isTileOccupiedList,
						    newExteriorCellsList,
						    materialCellLabels);

	    // Uncompress tiles of the new active cells
	    HDK::Utilities::uncompressTiles(materialCellLabels,
					    isTileOccupiedList);

	    HDK::Utilities::uncompressTiles(visitedCells,
					    isTileOccupiedList);
	
	    setSolidLayerCells(materialCellLabels, visitedCells, newExteriorCellsList, interiorMaterialLabel, exteriorMaterialLabel);

	    // Build next layer of active cells unless we're at the final layer
	    if (layer < solidBoundaryLayerSize - 1)
	    {
		buildNextSolidBoundaryLayer(parallelExteriorBandList,
					    newExteriorCellsList,
					    materialCellLabels,
					    visitedCells,
					    cutCellWeights,
					    interiorMaterialLabel,
					    exteriorMaterialLabel);
	    }
	}

    }
}

void
HDK_AffineBubblePressureSolver::buildTiledActiveCells(SIM_RawIndexField &materialCellLabels,
							const int tileSize,
							const exint exteriorMaterialLabel,
							const exint interiorMaterialLabel) const
{
    using SIM::FieldUtils::setFieldValue;

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

	    if (!vit.isTileConstant() || vit.getValue() == interiorMaterialLabel)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == interiorMaterialLabel)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());
			
			if (cell[0] % tileSize == 0 ||
			    cell[1] % tileSize == 0 ||
			    cell[2] % tileSize == 0)
			    setFieldValue(materialCellLabels, cell, exteriorMaterialLabel);
		    }
		}
	    }
	}
    });
}

void
HDK_AffineBubblePressureSolver::setAllInteriorCellsExterior(SIM_RawIndexField &materialCellLabels,
							    const exint exteriorMaterialLabel,
							    const exint interiorMaterialLabel) const
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

	    if (vit.isTileConstant() && vit.getValue() == interiorMaterialLabel)
		vit.getTile()->makeConstant(exteriorMaterialLabel);
	    else
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == interiorMaterialLabel)
			vit.setValue(exteriorMaterialLabel);
		}
	    }
	}
    });
}

void
HDK_AffineBubblePressureSolver::trackBubbleIDs(UT_Array<bool> &isBubbleTracked,
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

    UTparallelForLightItems(GA_SplittableRange(particles.getPointRange()), [&](const GA_SplittableRange &range)
    {
    	UT_Array<bool> localIsBubbleTracked;
	localIsBubbleTracked.setSize(isBubbleTracked.size());
	localIsBubbleTracked.constant(false);

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
    		    localIsBubbleTracked[newBubbleRegion] = true;

    		for (int axis : {0,1,2})
    		    for (int direction : {0,1})
    		    {
			UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

			if (adjacentCell[axis] < 0 || adjacentCell[axis] >= voxelRes[axis])
			    continue;

			exint newBubbleRegion = getFieldValue(bubbleRegionIndices, adjacentCell);

			if (newBubbleRegion >= 0)
			    localIsBubbleTracked[newBubbleRegion] = true;
    		    }
    	    }
    	}

    	for (exint newBubbleRegion = 0; newBubbleRegion < isBubbleTracked.size(); ++newBubbleRegion)
    	{
    	    if (localIsBubbleTracked[newBubbleRegion])
                isBubbleTracked[newBubbleRegion] = true;
    	}
    });
}

void
HDK_AffineBubblePressureSolver::updateBubbleIDs(GA_RWHandleI &bubbleIDHandle,
						const UT_Array<bool> &isBubbleTracked,
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
		if (bubbleRegion >= 0 && isBubbleTracked[bubbleRegion])
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

			    if (bubbleRegion >= 0 && isBubbleTracked[bubbleRegion])
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
HDK_AffineBubblePressureSolver::setUntrackedBubbleCells(SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &bubbleRegionIndices,
							const UT_Array<bool> &isBubbleTracked) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForEachNumber(materialCellLabels.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
    {
	UT_VoxelArrayIteratorI vit(materialCellLabels.fieldNC());

    	if (boss->opInterrupt())
	    return;
	
	for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	{
	    vit.myTileStart = tileNumber;
	    vit.myTileEnd = tileNumber + 1;
	    vit.rewind();

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());

			exint bubbleRegion = getFieldValue(bubbleRegionIndices, cell);
			assert(bubbleRegion >= 0);

			if (!isBubbleTracked[bubbleRegion])
			    vit.setValue(MaterialLabels::EXTERIOR_BUBBLE_CELL);
		    }
		}
	    }
	}
    });
}

void
HDK_AffineBubblePressureSolver::buildBubbleRegions(SIM_RawIndexField &fullBubbleRegionIndices,
						    UT_Array<bool> &isBubbleTracked,
						    SIM_RawIndexField &materialCellLabels,
						    SIM_Object *obj,
						    const SIM_RawField &liquidSurface,
						    std::array<const SIM_RawField *, 3> &cutCellWeights)
{
    // First build combined region of bubble cells
    SIM_VolumetricConnectedComponentBuilder fullBubbleRegionBuilder(fullBubbleRegionIndices, materialCellLabels, cutCellWeights.data());

    exint fullBubbleRegionCount = fullBubbleRegionBuilder.buildConnectedComponents([](const exint label)
    { 
	return label == MaterialLabels::INTERIOR_BUBBLE_CELL || label == MaterialLabels::EXTERIOR_BUBBLE_CELL;
    });

    isBubbleTracked.setSize(fullBubbleRegionCount);
    isBubbleTracked.constant(false);

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

    // Track bubble IDs by filling in array. If a particle with a valid bubble ID is adjacent to a
    // newly detected bubble, indicate that the new bubble is indeed tracked.
    {
	trackBubbleIDs(isBubbleTracked,
			bubbleIDHandle,
			*particles,
			fullBubbleRegionIndices,
			liquidSurface);

	// Update bubble ID on particles for particles adjacent to bubbles that were tracked
	// from the previous timestep.
	updateBubbleIDs(bubbleIDHandle, isBubbleTracked, *particles, fullBubbleRegionIndices, liquidSurface);
	bubbleIDHandle.bumpDataId();
    }

    // Now set any untracked bubble cells to exterior
    setUntrackedBubbleCells(materialCellLabels, fullBubbleRegionIndices, isBubbleTracked);

    std::cout << "    Bubble regions: " << fullBubbleRegionCount << std::endl;
}

void
HDK_AffineBubblePressureSolver::buildInteriorBoundingBoxes(UT_Array<UT_Array<UT_Vector3I>> &parallelInteriorRegionBBMin,
							    UT_Array<UT_Array<UT_Vector3I>> &parallelInteriorRegionBBMax,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &interiorRegionIndices) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInteriorBoundingBoxesAlgorithm;
    buildInteriorBoundingBoxesAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_Array<UT_Vector3I> &localInteriorRegionBBMin = parallelInteriorRegionBBMin[info.job()];
	UT_Array<UT_Vector3I> &localInteriorRegionBBMax = parallelInteriorRegionBBMax[info.job()];

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL ||
		vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL ||
			vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
			exint interiorRegion = getFieldValue(interiorRegionIndices, cell);

			assert(interiorRegion >= 0);

			for (int axis : {0,1,2})
			{
			    localInteriorRegionBBMin[interiorRegion][axis] = std::min(localInteriorRegionBBMin[interiorRegion][axis], cell[axis]);
			    localInteriorRegionBBMax[interiorRegion][axis] = std::max(localInteriorRegionBBMax[interiorRegion][axis], cell[axis]);
			}
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::remapInteriorRegions(SIM_RawIndexField &interiorRegionIndices,
						    SIM_RawIndexField &materialCellLabels,
						    const UT_Array<bool> &doRemoveRegion,
						    const UT_Array<exint> &remapRegion) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

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
		vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL ||
		vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL ||
			vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());

			exint interiorRegion = getFieldValue(interiorRegionIndices, cell);
			assert(interiorRegion >= 0);

			if (doRemoveRegion[interiorRegion])
			{
			    assert(remapRegion[interiorRegion] == -1);
			    if (vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
				setFieldValue(materialCellLabels, cell, MaterialLabels::EXTERIOR_BUBBLE_CELL);
			    else
				setFieldValue(materialCellLabels, cell, MaterialLabels::EXTERIOR_LIQUID_CELL);

			    setFieldValue(interiorRegionIndices, cell, HDK::Utilities::UNLABELLED_CELL);
			}
			else if (remapRegion[interiorRegion] < interiorRegion)
			    setFieldValue(interiorRegionIndices, cell, remapRegion[interiorRegion]);
			else assert(remapRegion[interiorRegion] == interiorRegion);
		    }
		}
	    }
	}
    });
}

exint
HDK_AffineBubblePressureSolver::removeDegenerateCellRegions(SIM_RawIndexField &interiorRegionIndices,
							    SIM_RawIndexField &materialCellLabels,
							    const exint interiorRegionCount) const
{
    const int threadCount = UT_Thread::getNumProcessors();

    exint tbbGrainSize = interiorRegionCount / (4 * threadCount);

    UT_Array<UT_Array<UT_Vector3I>> parallelInteriorRegionBBMin;
    UT_Array<UT_Array<UT_Vector3I>> parallelInteriorRegionBBMax;

    parallelInteriorRegionBBMin.setSize(threadCount);
    parallelInteriorRegionBBMax.setSize(threadCount);

    for (int thread = 0; thread < threadCount; ++thread)
    {
	parallelInteriorRegionBBMin[thread].setSize(interiorRegionCount);
	parallelInteriorRegionBBMin[thread].constant(UT_Vector3I(std::numeric_limits<exint>::max()));

	parallelInteriorRegionBBMax[thread].setSize(interiorRegionCount);
	parallelInteriorRegionBBMax[thread].constant(UT_Vector3I(std::numeric_limits<exint>::min()));
    }

    buildInteriorBoundingBoxes(parallelInteriorRegionBBMin, parallelInteriorRegionBBMax, materialCellLabels, interiorRegionIndices);

    UT_Array<UT_Vector3I> interiorRegionBBMin;
    interiorRegionBBMin.setSize(interiorRegionCount);
    interiorRegionBBMin.constant(UT_Vector3I(std::numeric_limits<exint>::max()));

    UT_Array<UT_Vector3I> interiorRegionBBMax;
    interiorRegionBBMax.setSize(interiorRegionCount);
    interiorRegionBBMax.constant(UT_Vector3I(std::numeric_limits<exint>::min()));

    tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
    {
	for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
	{
	    for (int thread = 0; thread < threadCount; ++thread)
	    {
		for (int axis : {0,1,2})
		{
		    interiorRegionBBMin[interiorRegion][axis] = std::min(parallelInteriorRegionBBMin[thread][interiorRegion][axis], interiorRegionBBMin[interiorRegion][axis]);
		    interiorRegionBBMax[interiorRegion][axis] = std::max(parallelInteriorRegionBBMax[thread][interiorRegion][axis], interiorRegionBBMax[interiorRegion][axis]);
		}
	    }
	}
    });

    UT_Array<bool> doRemoveRegion;
    doRemoveRegion.setSize(interiorRegionCount);
    doRemoveRegion.constant(false);

    tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
    {
	for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
	{
	    for (int axis : {0,1,2})
	    {
		if (interiorRegionBBMax[interiorRegion][axis] == interiorRegionBBMin[interiorRegion][axis])
		    doRemoveRegion[interiorRegion] = true;
	    }
	}
    });

    UT_Array<exint> remapRegion;
    remapRegion.setSize(interiorRegionCount);
    remapRegion.constant(-1);

    int regionMap = 0;
    for (exint interiorRegion = 0; interiorRegion < interiorRegionCount; ++interiorRegion)
    {
	if (!doRemoveRegion[interiorRegion])
	    remapRegion[interiorRegion] = regionMap++;
    }

    exint newRegionCount = interiorRegionCount;
    if (newRegionCount > regionMap)
    {
	newRegionCount = regionMap;

	remapInteriorRegions(interiorRegionIndices,
				materialCellLabels,
				doRemoveRegion,
				remapRegion);
    }

    return newRegionCount;
}

exint
HDK_AffineBubblePressureSolver::buildActiveCellIndices(SIM_RawIndexField &activeCellIndices,
							const SIM_RawIndexField &materialCellLabels) const
{
    using SIM::FieldUtils::setFieldValue;

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
	    vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
	    vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
	{
	    vitt.setTile(vit);

	    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
	    {
		if (vitt.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
		    vitt.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
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
HDK_AffineBubblePressureSolver::printCellGeometry(SIM_Object *obj,
						    const SIM_RawIndexField &activeCellIndices,
						    const SIM_RawIndexField &interiorRegionIndices,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_VectorField &velocity,
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
	    vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
	    vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
	    vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
	    vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
	{
	    vitt.setTile(vit);

	    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
	    {
		if (vitt.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
		    vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
		    vitt.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
		    vitt.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

		    GA_Offset pointOffset = cellActivityDetail->appendPoint();
		    cellSizeHandle.set(pointOffset, dx);

		    if (vitt.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
		    {
			assert(getFieldValue(activeCellIndices, cell) >= 0);
			assert(getFieldValue(interiorRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			cellActivityHandle.set(pointOffset, 0);
		    }
		    else if (vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
		    {
			assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			assert(getFieldValue(interiorRegionIndices, cell) >= 0);
			cellActivityHandle.set(pointOffset, 1);
		    }
		    else if (vitt.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
		    {
			assert(getFieldValue(activeCellIndices, cell) >= 0);
			assert(getFieldValue(interiorRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			cellActivityHandle.set(pointOffset, 2);
		    }
		    else
		    {
			assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			assert(getFieldValue(interiorRegionIndices, cell) >= 0);
			cellActivityHandle.set(pointOffset, 3);
		    }

		    UT_Vector3 cellCenter = velocity.getOrig() + UT_Vector3(vitt.x() + .5, vitt.y() + .5, vitt.z() + .5) * dx;
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
HDK_AffineBubblePressureSolver::buildInteriorRegionCOM(UT_Array<UT_Array<UT_Vector3SR>> &parallelInteriorRegionCOM,
							UT_Array<UT_Array<SolveReal>> &parallelInteriorRegionCellCount,
							const SIM_RawIndexField &interiorRegionIndices,
							const SIM_RawIndexField &materialCellLabels) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();
  
    UT_ThreadedAlgorithm buildInteriorRegionCOMAlgorithm;
    buildInteriorRegionCOMAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_Array<UT_Vector3SR> &localInteriorRegionCOM = parallelInteriorRegionCOM[info.job()];
	UT_Array<SolveReal> &localInteriorRegionCellCount = parallelInteriorRegionCellCount[info.job()];

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
			vitt.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			exint interiorRegion = getFieldValue(interiorRegionIndices, cell);
			assert(interiorRegion >= 0);

			++localInteriorRegionCellCount[interiorRegion];
			localInteriorRegionCOM[interiorRegion] += UT_Vector3SR(cell);
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::buildBestFitSystems(UT_Array<UT_Array<MassMatrix>> &parallelBestFitMatrix,
						    UT_Array<UT_Array<ColumnVector>> &parallelBestFitRHS,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_RawIndexField &interiorRegionIndices,
						    const SIM_VectorField &velocity,
						    const UT_Array<UT_Vector3SR> &interiorRegionCOM,
						    const SolveReal dx) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildBestFitSystemAlgorithm;
    buildBestFitSystemAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(interiorRegionIndices.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<MassMatrix> &localBestFitMatrix = parallelBestFitMatrix[info.job()];
	UT_Array<ColumnVector> &localBestFitRHS = parallelBestFitRHS[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() || vit.getValue() >= 0)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

		    exint interiorRegion = vitt.getValue();
		    if (interiorRegion >= 0)
		    {
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_BUBBLE_CELL ||
				getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL);

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
				assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < interiorRegionIndices.getVoxelRes()[axis]);

				if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
				    getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
				{
				    UT_Vector3SR offset(cell);

				    if (direction == 0)
					offset[axis] -= .5;
				    else
					offset[axis] += .5;
				    
				    offset *= dx;
				    offset -= interiorRegionCOM[interiorRegion];

				    ColumnVector columnVector = ColumnVector::Zero();

				    // Account for linear velocity
				    columnVector(axis) = 1;

				    // Account for affine matrix
				    columnVector(3 + 3 * axis) = offset[0];
				    columnVector(3 + 3 * axis + 1) = offset[1];

				    if (axis == 2)
				    {
					// Account for C_33 = -(C_11 + C_22)
					columnVector(3) = -offset[2];
					columnVector(7) = -offset[2];
				    }
				    else
					columnVector(3 + 3 * axis + 2) = offset[2];

				    localBestFitMatrix[interiorRegion] += columnVector * columnVector.transpose();
				    localBestFitRHS[interiorRegion] += getFieldValue(*velocity.getField(axis), cellToFaceMap(cell, axis, direction)) * columnVector;
				}
			    }
		    }
		    else assert(getFieldValue(materialCellLabels, cell) != MaterialLabels::INTERIOR_BUBBLE_CELL &&
				getFieldValue(materialCellLabels, cell) != MaterialLabels::INTERIOR_LIQUID_CELL);
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::buildRHS(Vector &rhsVector,
					    const SIM_RawIndexField &materialCellLabels,
					    const SIM_RawIndexField &interiorRegionIndices,	
					    const SIM_RawIndexField &activeCellIndices,
					    const SIM_VectorField &velocity,
					    const std::array<const SIM_RawField *, 3> &cutCellWeights,
					    const SIM_VectorField *solidVelocity,
					    const UT_Array<UT_Vector3SR> &interiorLinearVelocities,
					    const UT_Array<GradientMatrix> &interiorVelocityGradients,
					    const UT_Array<UT_Vector3SR> &interiorRegionCOM,
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
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_LIQUID_CELL ||
				getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_BUBBLE_CELL);

			assert(getFieldValue(interiorRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL);

			SolveReal divergence = 0;

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I face = cellToFaceMap(cell, axis, direction);
				SolveReal weight = getFieldValue(*cutCellWeights[axis], face);

				SolveReal sign = (direction == 0) ? 1 : -1;

				if (weight > 0)
				{
				    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				    SolveReal localVelocity;

				    if (adjacentCell[axis] >= 0 && adjacentCell[axis] < activeCellIndices.getVoxelRes()[axis] &&
					(getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_BUBBLE_CELL ||
					    getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL))
				    {
					exint interiorRegion = getFieldValue(interiorRegionIndices, adjacentCell);
					assert(interiorRegion >= 0);

					assert(getFieldValue(activeCellIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);

					UT_Vector3T<SolveReal> offset(cell);
					
					if (direction == 0)
					    offset[axis] -= .5;
					else
					    offset[axis] += .5;

					offset *= dx;
					offset -= interiorRegionCOM[interiorRegion];

					localVelocity = interiorLinearVelocities[interiorRegion][axis];

					for (int offsetAxis : {0,1,2})
					    localVelocity += interiorVelocityGradients[interiorRegion](axis, offsetAxis) * offset[offsetAxis];
				    }
				    else localVelocity = getFieldValue(*velocity.getField(axis), face);

				    divergence += sign * weight * localVelocity;
				}

				if (solidVelocity != nullptr && weight < 1)
				{
				    UT_Vector3 point;
				    velocity.getField(axis)->indexToPos(face[0], face[1], face[2], point);

				    divergence += sign * (1. - weight) * solidVelocity->getField(axis)->getValue(point);
				}
			    }

			rhsVector(index) = divergence;
		    }
		    else assert(getFieldValue(materialCellLabels, cell) != MaterialLabels::EXTERIOR_LIQUID_CELL &&
				getFieldValue(materialCellLabels, cell) != MaterialLabels::EXTERIOR_BUBBLE_CELL);
		}
	    }
	}
    });
}

void
HDK_AffineBubblePressureSolver::applyGoalDivergence(Vector &rhsVector,
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

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
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
HDK_AffineBubblePressureSolver::buildActiveRegionDivergence(UT_Array<UT_Array<SolveReal>> &parallelActiveRegionDivergence,
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
		    else assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::removeAverageDivergence(Vector &rhsVector,
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
		    else assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
		}
	    }
	}
    });

}

void
HDK_AffineBubblePressureSolver::setUntrackedDivergence(UT_Array<UT_Array<SolveReal>> &parallelUntrackedDivergence,
							UT_Array<UT_Array<exint>> &parallelTrackedCellCount,
							UT_Array<UT_Array<exint>> &parallelUntrackedCellCount,
							Vector &rhsVector,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &activeCellIndices,
							const SIM_RawIndexField &activeRegionIndices,
							const SIM_RawIndexField &fullBubbleRegionIndices,
							const UT_Array<bool> &isBubbleTracked,
							const SolveReal dx, const SolveReal dt) const
{
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm setUntrackedDivergenceAlgorithm;
    setUntrackedDivergenceAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<SolveReal> &localUntrackedDivergence = parallelUntrackedDivergence[info.job()];
	UT_Array<exint> &localTrackedCellCount = parallelTrackedCellCount[info.job()];
	UT_Array<exint> &localUntrackedCellCount = parallelUntrackedCellCount[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

		    if (vitt.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
		    {
			exint activeRegion = getFieldValue(activeRegionIndices, cell);
			assert(activeRegion >= 0);

			++localTrackedCellCount[activeRegion];
		    }
		    else if (vitt.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
		    {
			exint bubbleRegion = getFieldValue(fullBubbleRegionIndices, cell);

			exint activeRegion = getFieldValue(activeRegionIndices, cell);
			assert(activeRegion >= 0);

			if (isBubbleTracked[bubbleRegion])
			    ++localTrackedCellCount[activeRegion];
			else
			{
			    SolveReal localDivergence = -dx / dt;
			    localUntrackedDivergence[activeRegion] += localDivergence;

			    exint index = getFieldValue(activeCellIndices, cell);
			    assert(index >= 0);

			    rhsVector(index) += localDivergence;

			    ++localUntrackedCellCount[bubbleRegion];
			}
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::removeUntrackedDivergence(Vector &rhsVector,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &activeCellIndices,
							    const SIM_RawIndexField &activeRegionIndices,
							    const SIM_RawIndexField &fullBubbleRegionIndices,
							    const UT_Array<bool> &isBubbleTracked,
							    const UT_Array<SolveReal> &untrackedDivergence) const
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
		vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    UT_Vector3I cell(vit.x(), vit.y(), vit.z());

		    if (vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
		    {
			exint index = getFieldValue(activeCellIndices, cell);
			assert(index >= 0);

			exint activeRegion = getFieldValue(activeRegionIndices, cell);
			assert(activeRegion >= 0);

			rhsVector(index) -= untrackedDivergence[activeRegion];
		    }
		    else if (vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
		    {
			exint bubbleRegion = getFieldValue(fullBubbleRegionIndices, cell);
			assert(bubbleRegion >= 0);

			if (isBubbleTracked[bubbleRegion])
			{
			    exint activeRegion = getFieldValue(activeRegionIndices, cell);
			    assert(activeRegion >= 0);

			    exint index = getFieldValue(activeCellIndices, cell);
			    assert(index >= 0);

			    rhsVector(index) -= untrackedDivergence[activeRegion];
			}
		    }
		}
	    }
	}
    });
}

void
HDK_AffineBubblePressureSolver::buildLinearSystem(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelPoissonElements,
						    std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelDiagonalPrecondElements,
						    std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelAffineBubbleCollectionElements,
						    std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelAffineLiquidCollectionElements,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_RawIndexField &interiorRegionIndices,
						    const SIM_RawIndexField &activeCellIndices,
						    const std::array<const SIM_RawField *, 3> &cutCellWeights,
						    const SIM_RawField &liquidSurface,
						    const SolveReal liquidDensity,
						    const SolveReal bubbleDensity,
						    const UT_Array<UT_Vector3SR> &interiorRegionCOM,
						    const SolveReal dx) const
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
	std::vector<Eigen::Triplet<SolveReal>> &localAffineBubbleCollectionElements = parallelAffineBubbleCollectionElements[info.job()];
    	std::vector<Eigen::Triplet<SolveReal>> &localAffineLiquidCollectionElements = parallelAffineLiquidCollectionElements[info.job()];

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

			assert(getFieldValue(interiorRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL);

			SolveReal diagonal = 0;

			if (getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
			{
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

					if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
					    getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
					{
					    exint adjacentIndex = getFieldValue(activeCellIndices, adjacentCell);
					    assert(adjacentIndex >= 0);

					    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
						weight /= liquidDensity;
					    else
					    {
						SolveReal phi0 = getFieldValue(liquidSurface, cell);
						SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

						SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
						theta = SYSclamp(theta, SolveReal(0), SolveReal(1));

						weight /= (theta * liquidDensity + (1. - theta) * bubbleDensity);
					    }

					    diagonal += weight;
					    localPoissonElements.emplace_back(index, adjacentIndex, -weight);
					}
					else
					{
					    exint interiorRegion = getFieldValue(interiorRegionIndices, adjacentCell);
					    assert(interiorRegion >= 0);
					    assert(weight == 1);

					    assert(getFieldValue(activeCellIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);
					    assert(getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL);

					    UT_Vector3T<SolveReal> offset(cell);

					    if (direction == 0)
						offset[axis] -= .5;
					    else
						offset[axis] += .5;
	    
					    offset *= dx;
					    offset -= interiorRegionCOM[interiorRegion];

					    SolveReal sign = (direction == 0) ? -1 : 1;
					    SolveReal coeff = sign / liquidDensity;

					    // Account for linear velocity
					    localAffineLiquidCollectionElements.emplace_back(11 * interiorRegion + axis, index, coeff);

					    // Account for affine velocity offset
					    localAffineLiquidCollectionElements.emplace_back(11 * interiorRegion + 3 + 3 * axis, index, coeff * offset[0]);
					    localAffineLiquidCollectionElements.emplace_back(11 * interiorRegion + 3 + 3 * axis + 1, index, coeff * offset[1]);

					    // Account for C_33 = -(C_11 + C_22)
					    if (axis == 2)
					    {
						localAffineLiquidCollectionElements.emplace_back(11 * interiorRegion + 3, index, -coeff * offset[2]);
						localAffineLiquidCollectionElements.emplace_back(11 * interiorRegion + 7, index, -coeff * offset[2]);
					    }
					    else
						localAffineLiquidCollectionElements.emplace_back(11 * interiorRegion + 3 + 3 * axis + 2, index, coeff * offset[2]);
					}
				    }
				}
			}
			else
			{
			    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_BUBBLE_CELL);
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

					if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
					    getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
					{
					    exint adjacentIndex = getFieldValue(activeCellIndices, adjacentCell);
					    assert(adjacentIndex >= 0);

					    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL)
						weight /= bubbleDensity;
					    else
					    {
						SolveReal phi0 = getFieldValue(liquidSurface, cell);
						SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

						SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
						theta = SYSclamp(theta, SolveReal(0), SolveReal(1));

						weight /= (theta * liquidDensity + (1. - theta) * bubbleDensity);
					    }

					    diagonal += weight;
					    localPoissonElements.emplace_back(index, adjacentIndex, -weight);
					}
					else
					{
					    exint interiorRegion = getFieldValue(interiorRegionIndices, adjacentCell);
					    assert(interiorRegion >= 0);
					    assert(weight == 1);

					    assert(getFieldValue(activeCellIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);
					    assert(getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_BUBBLE_CELL);

					    UT_Vector3T<SolveReal> offset(cell);

					    if (direction == 0)
						offset[axis] -= .5;
					    else
						offset[axis] += .5;
	    
					    offset *= dx;
					    offset -= interiorRegionCOM[interiorRegion];

					    SolveReal sign = (direction == 0) ? -1 : 1;
					    SolveReal coeff = sign / bubbleDensity;

					    // Account for linear velocity
					    localAffineBubbleCollectionElements.emplace_back(11 * interiorRegion + axis, index, coeff);

					    // Account for affine velocity offset
					    localAffineBubbleCollectionElements.emplace_back(11 * interiorRegion + 3 + 3 * axis, index, coeff * offset[0]);
					    localAffineBubbleCollectionElements.emplace_back(11 * interiorRegion + 3 + 3 * axis + 1, index, coeff * offset[1]);

					    // Account for C_33 = -(C_11 + C_22)
					    if (axis == 2)
					    {
						localAffineBubbleCollectionElements.emplace_back(11 * interiorRegion + 3, index, -coeff * offset[2]);
						localAffineBubbleCollectionElements.emplace_back(11 * interiorRegion + 7, index, -coeff * offset[2]);
					    }
					    else
						localAffineBubbleCollectionElements.emplace_back(11 * interiorRegion + 3 + 3 * axis + 2, index, coeff * offset[2]);
					}
				    }
				}
			}

			assert(diagonal > 0);
			localPoissonElements.emplace_back(index, index, diagonal);
			localDiagonalPrecondElements.emplace_back(index, index, 1. / diagonal);
		    }
		}
	    }
	}

	return 0;
    });    
}

void
HDK_AffineBubblePressureSolver::buildAffineMassMatrixSystems(UT_Array<UT_Array<MassMatrix>> &parallelAffineMassMatrix,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &interiorRegionIndices,
							    const UT_Array<UT_Vector3SR> &interiorRegionCOM,
							    const SolveReal dx) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildBestFitSystemAlgorithm;
    buildBestFitSystemAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(interiorRegionIndices.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	UT_Array<MassMatrix> &localMassMatrix = parallelAffineMassMatrix[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		break;

	    if (!vit.isTileConstant() || vit.getValue() >= 0)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

		    exint interiorRegion = vitt.getValue();
		    if (interiorRegion >= 0)
		    {
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_BUBBLE_CELL ||
				getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL);

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				bool doApplyFace = false;
				if (direction == 0)
				    doApplyFace = true;
				else
				{
				    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
				    assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < interiorRegionIndices.getVoxelRes()[axis]);

				    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
					getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
					doApplyFace = true;
				    else assert(getFieldValue(interiorRegionIndices, adjacentCell) == interiorRegion);
				}

				if (doApplyFace)
				{
				    UT_Vector3SR offset(cell);

				    if (direction == 0)
					offset[axis] -= .5;
				    else
					offset[axis] += .5;
				    
				    offset *= dx;
				    offset -= interiorRegionCOM[interiorRegion];

				    ColumnVector columnVector = ColumnVector::Zero();

				    // Account for linear velocity
				    columnVector(axis) = 1;

				    // Account for affine matrix
				    columnVector(3 + 3 * axis) = offset[0];
				    columnVector(3 + 3 * axis + 1) = offset[1];

				    if (axis == 2)
				    {
					// Account for C_33 = -(C_11 + C_22)
					columnVector(3) = -offset[2];
					columnVector(7) = -offset[2];
				    }
				    else
					columnVector(3 + 3 * axis + 2) = offset[2];

				    localMassMatrix[interiorRegion] += columnVector * columnVector.transpose();
				}
			    }
		    }
		    else assert(getFieldValue(materialCellLabels, cell) != MaterialLabels::INTERIOR_BUBBLE_CELL &&
				getFieldValue(materialCellLabels, cell) != MaterialLabels::INTERIOR_LIQUID_CELL);
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::buildInteriorBoundaryCells(UT_Array<UT_Array<UT_Vector3I>> &parallelInteriorBoundaryCells,
							    const SIM_RawIndexField &materialCellLabels) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::getFieldValue;    

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInteriorBoundaryCellsAlgorithm;
    buildInteriorBoundaryCellsAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_Array<UT_Vector3I> &localInteriorBoundaryCells = parallelInteriorBoundaryCells[info.job()];

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		break;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
			vitt.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			bool isBoundaryCell = false;
			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
				    continue;

				if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_BUBBLE_CELL ||
				    getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL)
				    isBoundaryCell = true;
			    }

			if (isBoundaryCell)
			    localInteriorBoundaryCells.append(cell);
		    }
		}
	    }
	}

    	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::distributeAffineVectors(Vector &destinationVector,
							const Vector &invertedAffineVectors,
							const UT_Array<UT_Vector3I> &interiorBoundaryCells,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &interiorRegionIndices,
							const SIM_RawIndexField &activeCellIndices,
							const UT_Array<UT_Vector3SR> &interiorRegionCOM,
							const SolveReal dx) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();
    
    UTparallelForLightItems(UT_BlockedRange<exint>(0, interiorBoundaryCells.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;
    
	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = interiorBoundaryCells[cellIndex];

	    exint activeIndex = getFieldValue(activeCellIndices, cell);
	    assert(activeIndex >= 0);
	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
		    getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_LIQUID_CELL);

	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

		    assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < materialCellLabels.getVoxelRes()[axis]);

		    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_BUBBLE_CELL ||
			getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL)
		    {
			exint interiorRegion = getFieldValue(interiorRegionIndices, adjacentCell);
			assert(interiorRegion >= 0);

			assert(getFieldValue(activeCellIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);

			UT_Vector3T<SolveReal> offset(cell);

			if (direction == 0)
			    offset[axis] -= .5;
			else
			    offset[axis] += .5;
	    
			offset *= dx;
			offset -= interiorRegionCOM[interiorRegion];

			SolveReal sign = (direction == 0) ? -1 : 1;

			// Account for linear velocity
			destinationVector(activeIndex) += sign * invertedAffineVectors(11 * interiorRegion + axis);

			// Account for affine velocity offset
			destinationVector(activeIndex) += sign * offset[0] * invertedAffineVectors(11 * interiorRegion + 3 + 3 * axis);
			destinationVector(activeIndex) += sign * offset[1] * invertedAffineVectors(11 * interiorRegion + 3 + 3 * axis + 1);

			// Account for C_33 = -(C_11 + C_22)
			if (axis == 2)
			{
			    destinationVector(activeIndex) -= sign * offset[2] * invertedAffineVectors(11 * interiorRegion + 3);
			    destinationVector(activeIndex) -= sign * offset[2] * invertedAffineVectors(11 * interiorRegion + 7);
			}
			else
			    destinationVector(activeIndex) += sign * offset[2] * invertedAffineVectors(11 * interiorRegion + 3 + 3 * axis + 2);
		    }
		}
	}

	return;
    });
}

//
// Multigrid operations
//

void
HDK_AffineBubblePressureSolver::buildMGDomainLabels(UT_VoxelArray<int> &mgDomainCellLabels,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_RawIndexField &activeCellIndices) const
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
		vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
		vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    UT_Vector3I cell(vit.x(), vit.y(), vit.z());

		    if (vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
		    {
			assert(getFieldValue(activeCellIndices, cell) >= 0);
			mgDomainCellLabels.setValue(cell, MGCellLabels::INTERIOR_CELL);
		    }
		    else if (vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
		    {
			assert(getFieldValue(activeCellIndices, cell) >= 0);
			mgDomainCellLabels.setValue(cell, MGCellLabels::DIRICHLET_CELL);
		    }
		    else if (vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
		    {
			assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
			mgDomainCellLabels.setValue(cell, MGCellLabels::DIRICHLET_CELL);
		    }
		    else
		    {
			assert(mgDomainCellLabels(cell) == MGCellLabels::EXTERIOR_CELL);
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::SOLID_CELL);
			assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
		    }
		}
	    }
	    else assert(vit.getValue() != MaterialLabels::INTERIOR_LIQUID_CELL);
	}
    });
}

void
HDK_AffineBubblePressureSolver::buildMGBoundaryWeights(UT_VoxelArray<SolveReal> &boundaryWeights,
							const SIM_RawField &validFaces,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &interiorRegionIndices,
							const UT_VoxelArray<int> &domainCellLabels,
							const SIM_RawField &liquidSurface,
							const SIM_RawField &cutCellWeights,
							const SolveReal liquidDensity,
							const SolveReal bubbleDensity,
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
			    assert(getFieldValue(materialCellLabels, backwardCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
				    getFieldValue(materialCellLabels, forwardCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL);

			    SolveReal phi0 = getFieldValue(liquidSurface, backwardCell);
			    SolveReal phi1 = getFieldValue(liquidSurface, forwardCell);

			    SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
			    theta = SYSclamp(theta, .01, 1.);

			    weight *= liquidDensity / (theta * liquidDensity + (1. - theta) * bubbleDensity);
			}

			boundaryWeights.setValue(face, weight);
		    }
		}
	    }
	}
    });
}

bool
HDK_AffineBubblePressureSolver::unitTestMGLabels(const UT_VoxelArray<int> &mgDomainCellLabels,
						const SIM_RawIndexField &materialCellLabels,
						const SIM_RawIndexField &activeCellIndices,
						const SIM_RawIndexField &interiorRegionIndices,
						const UT_Vector3I &mgExpandedOffset) const
{
    using SIM::FieldUtils::getFieldValue;

    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

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

		if (vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
		{
		    if (mgDomainCellLabels(expandedCell) != MGCellLabels::INTERIOR_CELL &&
			mgDomainCellLabels(expandedCell) != MGCellLabels::BOUNDARY_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(interiorRegionIndices, cell) != HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		}
		else if (vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
		{
		    if (mgDomainCellLabels(expandedCell) != MGCellLabels::DIRICHLET_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(activeCellIndices, cell) != HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(interiorRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		}
		else if (vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
		{
		    if (mgDomainCellLabels(expandedCell) != MGCellLabels::DIRICHLET_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(interiorRegionIndices, cell) != HDK::Utilities::UNLABELLED_CELL)
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
		    if (getFieldValue(activeCellIndices, cell) != HDK::Utilities::UNLABELLED_CELL)
		    {
			assert(false);
			passedTest = false;
		    }
		    if (getFieldValue(interiorRegionIndices, cell) != HDK::Utilities::UNLABELLED_CELL)
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
HDK_AffineBubblePressureSolver::buildInitialLiquidSmootherLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidBandCells,
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

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
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
				    getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL)
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
HDK_AffineBubblePressureSolver::setVisitedLiquidCells(SIM_RawIndexField &visitedLiquidCells,
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
HDK_AffineBubblePressureSolver::buildNextLiquidSmootherLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidLayerCells,
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

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_LIQUID_CELL);
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
			if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_LIQUID_CELL &&
			    getFieldValue(visitedLiquidCells, adjacentCell) == UNVISITED_CELL)
			    localLiquidLayerCells.append(adjacentCell);
		    }
		}
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::buildBubbleSmootherCells(UT_Array<UT_Array<UT_Vector3I>> &parallelBubbleSmootherCells,
							    UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidCopyCells,
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

	UT_Array<UT_Vector3I> &localBubbleSmootherCells = parallelBubbleSmootherCells[info.job()];
	UT_Array<UT_Vector3I> &localLiquidCopyCells = parallelLiquidCopyCells[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
		    if (vitt.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
		    {
			if (getFieldValue(visitedLiquidCells, cell) == VISITED_CELL)
			{
			    localBubbleSmootherCells.append(cell);
			    localLiquidCopyCells.append(cell);
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

				    if (getFieldValue(*cutCellWeights[axis], face) > 0 && 
					getFieldValue(visitedLiquidCells, adjacentCell) == VISITED_CELL)
					isSmootherBoundaryCell = true;
				}

			    if (isSmootherBoundaryCell)
				localLiquidCopyCells.append(cell);
			}			
		    }
		    else if (vitt.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL)
			localBubbleSmootherCells.append(cell);
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineBubblePressureSolver::findOccupiedTiles(UT_Array<bool> &isTileOccupied,
						const UT_VoxelArray<SolveReal> &voxelArrayGrid,
						const UT_Array<UT_Vector3I> &activeCells) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UTparallelForLightItems(UT_BlockedRange<exint>(0, activeCells.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;

	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = activeCells[cellIndex];

	    int tileNumber = voxelArrayGrid.indexToLinearTile(cell[0], cell[1], cell[2]);

	    if (!isTileOccupied[tileNumber])
		isTileOccupied[tileNumber] = true;
	}
    });
}

void
HDK_AffineBubblePressureSolver::copySourceToSmootherGrid(UT_VoxelArray<SolveReal> &smootherSourceGrid,
							const Vector &sourceVector,
							const SIM_RawIndexField &activeCellIndices,
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
	    exint index = getFieldValue(activeCellIndices, cell);
	    assert(index >= 0);
	    smootherSourceGrid.setValue(cell, sourceVector(index));
	}
    });
}

void
HDK_AffineBubblePressureSolver::buildAffineVectors(UT_Array<ColumnVector> &affineVectors,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_RawIndexField &interiorRegionIndices,
						    const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
						    const UT_Array<UT_Vector3I> &interiorBoundaryCells,
						    const UT_Array<UT_Vector3SR> &interiorRegionCOM,
						    const exint interiorRegionCount,
						    const SolveReal dx) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    tbb::enumerable_thread_specific<std::vector<ColumnVector>> parallelAffineVectors(interiorRegionCount, ColumnVector::Zero());

    tbb::parallel_for(tbb::blocked_range<exint>(0, interiorBoundaryCells.size(), 1000), [&](const tbb::blocked_range<exint> &range)
    {
	std::vector<ColumnVector> &localAffineVectors = parallelAffineVectors.local();

	if (boss->opInterrupt())
	    return;
	
	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = interiorBoundaryCells[cellIndex];

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_BUBBLE_CELL);

#if !defined(NDEBUG)
	    bool hasInteriorNeighbour = false;
#endif

	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

		    assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < materialCellLabels.getVoxelRes()[axis]);

		    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_BUBBLE_CELL)
		    {

#if !defined(NDEBUG)
			hasInteriorNeighbour = true;
#endif

			exint interiorRegion = getFieldValue(interiorRegionIndices, adjacentCell);
			assert(interiorRegion >= 0);

			UT_Vector3T<SolveReal> offset(cell);

			if (direction == 0)
			    offset[axis] -= .5;
			else
			    offset[axis] += .5;
	    
			offset *= dx;
			offset -= interiorRegionCOM[interiorRegion];

			SolveReal sign = (direction == 0) ? -1 : 1;

			SolveReal signedValue = sign * smootherDestinationGrid(cell);

			// Account for linear velocity
			localAffineVectors[interiorRegion](axis) += signedValue;

			// Account for affine velocity offset
			localAffineVectors[interiorRegion](3 + 3 * axis) += offset[0] * signedValue;
			localAffineVectors[interiorRegion](3 + 3 * axis + 1) += offset[1] * signedValue;

			// Account for C_33 = -(C_11 + C_22)
			if (axis == 2)
			{
			    localAffineVectors[interiorRegion](3) -= offset[2] * signedValue;
			    localAffineVectors[interiorRegion](7) -= offset[2] * signedValue;
			}
			else
			    localAffineVectors[interiorRegion](3 + 3 * axis + 2) += offset[2] * signedValue;
		    }
		}

	    assert(hasInteriorNeighbour);
	}

	return;
    });

    parallelAffineVectors.combine_each([&](const std::vector<ColumnVector> &localAffineVectors)
    {
	for (exint interiorRegion = 0; interiorRegion < interiorRegionCount; ++interiorRegion)
	    affineVectors[interiorRegion] += localAffineVectors[interiorRegion];
    });
}

void
HDK_AffineBubblePressureSolver::applyBubbleSmoother(UT_Array<SolveReal> &tempSmootherDestinationValues,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_RawIndexField &interiorRegionIndices,
						    const SIM_RawField &liquidSurface,
						    const std::array<const SIM_RawField *, 3> &cutCellWeights,
						    const UT_VoxelArray<int> &mgDomainCellLabels,
						    const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
						    const UT_VoxelArray<SolveReal> &smootherSourceGrid,
						    const UT_Array<ColumnVector> &invertedAffineVectors,
						    const UT_Array<UT_Vector3SR> &interiorRegionCOM,
						    const UT_Array<UT_Vector3I> &bubbleSmootherCells,
						    const UT_Vector3I &mgExpandedOffset,
						    const SolveReal liquidDensity,
						    const SolveReal bubbleDensity,
						    const SolveReal dx) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    using MGCellLabels = HDK::GeometricMultigridOperators::CellLabels;

    UT_Interrupt *boss = UTgetInterrupt();
    
    constexpr SolveReal dampedWeight = 2. / 3.;

    tbb::parallel_for(tbb::blocked_range<exint>(0, bubbleSmootherCells.size(), 1000), [&](const tbb::blocked_range<exint> &range)
    {
	if (boss->opInterrupt())
	    return;
	
	for (exint cellIndex = range.begin(); cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = bubbleSmootherCells[cellIndex];

	    SolveReal diagonal = 0;
	    SolveReal laplacian = 0;
	    if (getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
	    {
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

				if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
				    weight /= liquidDensity;
				else
				{
				    assert(getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL);

				    SolveReal phi0 = getFieldValue(liquidSurface, cell);
				    SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

				    assert(phi1 > 0);

				    SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
				    theta = SYSclamp(theta, SolveReal(0), SolveReal(1));

				    weight /= (theta * liquidDensity + (1. - theta) * bubbleDensity);
				}

				diagonal += weight;
				laplacian -= weight * smootherDestinationGrid(adjacentCell);
			    }
			}
		}
	    }
	    else
	    {
		assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_BUBBLE_CELL);
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

			    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
				getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
			    {
				if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL)
				    weight /= bubbleDensity;
				else
				{
				    SolveReal phi0 = getFieldValue(liquidSurface, cell);
				    SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

				    SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
				    theta = SYSclamp(theta, SolveReal(0), SolveReal(1));

				    weight /= (theta * liquidDensity + (1. - theta) * bubbleDensity);
				}

				diagonal += weight;
				laplacian -= weight * smootherDestinationGrid(adjacentCell);
			    }
			    else
			    {
				exint interiorRegion = getFieldValue(interiorRegionIndices, adjacentCell);
				assert(interiorRegion >= 0);
				assert(weight == 1);

				assert(getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_BUBBLE_CELL);

				UT_Vector3T<SolveReal> offset(cell);

				if (direction == 0)
				    offset[axis] -= .5;
				else
				    offset[axis] += .5;
	    
				offset *= dx;
				offset -= interiorRegionCOM[interiorRegion];

				SolveReal sign = (direction == 0) ? -1 : 1;

				// Account for linear velocity
				laplacian += sign * invertedAffineVectors[interiorRegion](axis);

				// Account for affine velocity offset
				laplacian += sign * offset[0] * invertedAffineVectors[interiorRegion](3 + 3 * axis);
				laplacian += sign * offset[1] * invertedAffineVectors[interiorRegion](3 + 3 * axis + 1);

				// Account for C_33 = -(C_11 + C_22)
				if (axis == 2)
				{
				    laplacian -= sign * offset[2] * invertedAffineVectors[interiorRegion](3);
				    laplacian -= sign * offset[2] * invertedAffineVectors[interiorRegion](7);
				}
				else
				    laplacian += sign * offset[2] * invertedAffineVectors[interiorRegion](3 + 3 * axis + 2);
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
HDK_AffineBubblePressureSolver::copyDestinationVectorToGrid(UT_VoxelArray<SolveReal> &smootherDestinationGrid,
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
HDK_AffineBubblePressureSolver::copySourceToMGGrid(UT_VoxelArray<SolveReal> &mgSourceGrid,
						const SIM_RawIndexField &materialCellLabels,
						const SIM_RawIndexField &activeCellIndices,
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

			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_LIQUID_CELL);

			exint activeIndex = getFieldValue(activeCellIndices, cell);
			assert(activeIndex >= 0);

			mgSourceGrid.setValue(expandedCell, liquidDensity * sourceVector(activeIndex));
		    }
		}
	    }
	}
    });
}

void
HDK_AffineBubblePressureSolver::applyDirichletToMGGrid(UT_VoxelArray<SolveReal> &mgSourceGrid,
							UT_VoxelArray<SolveReal> &mgDestinationGrid,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &activeCellIndices,
							const SIM_RawField &liquidSurface,
							const std::array<const SIM_RawField *, 3> &cutCellWeights,
							const UT_VoxelArray<int> &mgDomainCellLabels,
							const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
							const Vector &sourceVector,
							const UT_Array<UT_Vector3I> &bubbleSmootherCells,
							const UT_Vector3I &mgExpandedOffset,
							const SolveReal liquidDensity,
							const SolveReal bubbleDensity) const
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
	    
	    if (getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_LIQUID_CELL)
	    {
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

				if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::EXTERIOR_BUBBLE_CELL)
				{
				    SolveReal phi0 = getFieldValue(liquidSurface, cell);
				    SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

				    assert(phi1 > 0);

				    SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
				    theta = SYSclamp(theta, SolveReal(0), SolveReal(1));

				    weight *= liquidDensity / (theta * liquidDensity + (1. - theta) * bubbleDensity);

				    laplacian -= weight * smootherDestinationGrid(adjacentCell);
				}
			    }
			}

		    if (laplacian != 0)
		    {
			exint activeIndex = getFieldValue(activeCellIndices, cell);
			assert(activeIndex >= 0);
			mgSourceGrid.setValue(expandedCell, liquidDensity * sourceVector(activeIndex) - laplacian);
		    }
		}
	    }
	    else
	    {
		assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_BUBBLE_CELL);
		assert(mgDomainCellLabels(expandedCell) == MGCellLabels::DIRICHLET_CELL);
	    }
	}
    });
}

void
HDK_AffineBubblePressureSolver::copyMGToSmoother(UT_VoxelArray<SolveReal> &smootherDestinationGrid,
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

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::EXTERIOR_LIQUID_CELL);
	    assert(mgDomainCellLabels(expandedCell) == MGCellLabels::INTERIOR_CELL ||
		    mgDomainCellLabels(expandedCell) == MGCellLabels::BOUNDARY_CELL);

	    smootherDestinationGrid.setValue(cell, mgDestinationGrid(expandedCell));
	}
    });
}

void
HDK_AffineBubblePressureSolver::copyMGToDestinationVector(Vector &destinationVector,
    							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &activeCellIndices,
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

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());
			UT_Vector3I expandedCell = cell + mgExpandedOffset;

			exint activeIndex = getFieldValue(activeCellIndices, cell);
			assert(activeIndex >= 0);

			assert(mgDomainCellLabels(expandedCell) == MGCellLabels::BOUNDARY_CELL ||
				mgDomainCellLabels(expandedCell) == MGCellLabels::INTERIOR_CELL);

			destinationVector(activeIndex) = mgDestinationGrid(expandedCell);
		    }
		}
	    }
	}
    });
}

void
HDK_AffineBubblePressureSolver::copySmootherToDestinationVector(Vector &destinationVector,
								const SIM_RawIndexField &activeCellIndices,
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

	    exint activeIndex = getFieldValue(activeCellIndices, cell);
	    assert(activeIndex >= 0);

	    destinationVector(activeIndex) = smootherDestinationGrid(cell);
	}
    });
}


void
HDK_AffineBubblePressureSolver::buildValidFaces(SIM_VectorField &validFaces,
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
						[](const exint label){ return label == MaterialLabels::EXTERIOR_LIQUID_CELL ||
										label == MaterialLabels::INTERIOR_LIQUID_CELL ||
										label == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
										label == MaterialLabels::INTERIOR_BUBBLE_CELL; },
						axis);

	HDK::Utilities::uncompressTiles(*validFaces.getField(axis), isTileOccupiedList);

	HDK::Utilities::classifyValidFaces(*validFaces.getField(axis),
					    materialCellLabels,
					    *cutCellWeights[axis],
					    [](const exint label){ return label == MaterialLabels::EXTERIOR_LIQUID_CELL ||
									    label == MaterialLabels::INTERIOR_LIQUID_CELL ||
									    label == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
									    label == MaterialLabels::INTERIOR_BUBBLE_CELL; },
					    axis);
    }
}

void
HDK_AffineBubblePressureSolver::copyPressureVectorToGrid(SIM_RawField &pressure,
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
HDK_AffineBubblePressureSolver::applySolutionToVelocity(SIM_RawField &velocity,
							const SIM_RawField &validFaces,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &interiorRegionIndices,
							const SIM_RawIndexField &activeCellIndices,
							const SIM_RawField &pressure,
							const SIM_RawField &liquidSurface,
							const SIM_RawField &cutCellWeights,
							const UT_Array<UT_Vector3SR> &interiorRegionCOM,
							const UT_Array<UT_Vector3SR> &interiorLinearVelocities,
							const UT_Array<GradientMatrix> &interiorVelocityGradients,
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
		    if (vit.getValue() == HDK::Utilities::VALID_FACE)
		    {
			UT_Vector3I face(vit.x(), vit.y(), vit.z());

			UT_Vector3I backwardCell = faceToCellMap(face, axis, 0);
			UT_Vector3I forwardCell = faceToCellMap(face, axis, 1);

			assert(backwardCell[axis] >= 0 && forwardCell[axis] < activeCellIndices.getVoxelRes()[axis]);

			exint backwardMaterial = getFieldValue(materialCellLabels, backwardCell);
			exint forwardMaterial = getFieldValue(materialCellLabels, forwardCell);

			if (backwardMaterial == MaterialLabels::INTERIOR_BUBBLE_CELL ||
			    forwardMaterial == MaterialLabels::INTERIOR_BUBBLE_CELL ||
			    backwardMaterial == MaterialLabels::INTERIOR_LIQUID_CELL ||
			    forwardMaterial == MaterialLabels::INTERIOR_LIQUID_CELL)
			{
			    exint interiorRegion;

			    if (backwardMaterial == MaterialLabels::INTERIOR_BUBBLE_CELL ||
				backwardMaterial == MaterialLabels::INTERIOR_LIQUID_CELL)
			    {
				interiorRegion = getFieldValue(interiorRegionIndices, backwardCell);

				assert(interiorRegion >= 0);
				assert(getFieldValue(activeCellIndices, backwardCell) == HDK::Utilities::UNLABELLED_CELL);

				if (forwardMaterial == MaterialLabels::INTERIOR_BUBBLE_CELL ||
				    forwardMaterial == MaterialLabels::INTERIOR_LIQUID_CELL)
				{
				    if (backwardMaterial == MaterialLabels::INTERIOR_BUBBLE_CELL)
					assert(forwardMaterial == MaterialLabels::INTERIOR_BUBBLE_CELL);
				    else
					assert(forwardMaterial == MaterialLabels::INTERIOR_LIQUID_CELL);

				    assert(interiorRegion == getFieldValue(interiorRegionIndices, forwardCell));
				    assert(getFieldValue(activeCellIndices, forwardCell) == HDK::Utilities::UNLABELLED_CELL);
				}
				else
				{
				    if (backwardMaterial == MaterialLabels::INTERIOR_BUBBLE_CELL)
					assert(forwardMaterial == MaterialLabels::EXTERIOR_BUBBLE_CELL);
				    else
					assert(forwardMaterial == MaterialLabels::EXTERIOR_LIQUID_CELL);

				    assert(getFieldValue(interiorRegionIndices, forwardCell) == HDK::Utilities::UNLABELLED_CELL);
				    assert(getFieldValue(activeCellIndices, forwardCell) >= 0);
				}
			    }
			    else
			    {
				interiorRegion = getFieldValue(interiorRegionIndices, forwardCell);

				assert(interiorRegion >= 0);
				assert(getFieldValue(activeCellIndices, forwardCell) == HDK::Utilities::UNLABELLED_CELL);

				if (forwardMaterial == MaterialLabels::INTERIOR_BUBBLE_CELL)
				    assert(backwardMaterial == MaterialLabels::EXTERIOR_BUBBLE_CELL);
				else
				    assert(backwardMaterial == MaterialLabels::EXTERIOR_LIQUID_CELL);

				assert(getFieldValue(interiorRegionIndices, backwardCell) == HDK::Utilities::UNLABELLED_CELL);
				assert(getFieldValue(activeCellIndices, backwardCell) >= 0);
			    }

			    UT_Vector3SR offset(backwardCell);

			    offset[axis] += .5;
	    
			    offset *= dx;
			    offset -= interiorRegionCOM[interiorRegion];

			    SolveReal localVelocity = interiorLinearVelocities[interiorRegion][axis] +
							offset[0] * interiorVelocityGradients[interiorRegion](axis, 0) +
							offset[1] * interiorVelocityGradients[interiorRegion](axis, 1) +
							offset[2] * interiorVelocityGradients[interiorRegion](axis, 2);

			    setFieldValue(velocity, face, localVelocity);
			}
			else
			{
			    assert(backwardMaterial != MaterialLabels::SOLID_CELL && forwardMaterial != MaterialLabels::SOLID_CELL);

			    SolveReal localDensity;

			    if (backwardMaterial == MaterialLabels::EXTERIOR_LIQUID_CELL && forwardMaterial == MaterialLabels::EXTERIOR_LIQUID_CELL)
				localDensity = liquidDensity;
			    else if (backwardMaterial == MaterialLabels::EXTERIOR_BUBBLE_CELL && forwardMaterial == MaterialLabels::EXTERIOR_BUBBLE_CELL)
				localDensity = bubbleDensity;
			    else
			    {
				assert((backwardMaterial == MaterialLabels::EXTERIOR_LIQUID_CELL && forwardMaterial == MaterialLabels::EXTERIOR_BUBBLE_CELL) ||
					(backwardMaterial == MaterialLabels::EXTERIOR_BUBBLE_CELL && forwardMaterial == MaterialLabels::EXTERIOR_LIQUID_CELL));

				SolveReal phi0 = getFieldValue(liquidSurface, backwardCell);
				SolveReal phi1 = getFieldValue(liquidSurface, forwardCell);

				SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
				theta = SYSclamp(theta, SolveReal(0), SolveReal(1));
				localDensity = (theta * liquidDensity + (1. - theta) * bubbleDensity);
			    }

			    SolveReal gradient = getFieldValue(pressure, forwardCell) - getFieldValue(pressure, backwardCell);
			    gradient /= localDensity;

			    setFieldValue(velocity, face, getFieldValue(velocity, face) - gradient);
			}
		    }
		}
	    }
	}
    });
}

void
HDK_AffineBubblePressureSolver::computeResultingDivergence(UT_Array<SolveReal> &parallelAccumulatedDivergence,
							UT_Array<SolveReal> &parallelMaxDivergence,
							UT_Array<SolveReal> &parallelCellCount,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_VectorField &velocity,
							const std::array<const SIM_RawField *, 3> &cutCellWeights,
							const SIM_VectorField *solidVelocity) const
{
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

	SolveReal &localAccumulatedDivergence = parallelAccumulatedDivergence[info.job()];
	SolveReal &localMaxDivergence = parallelMaxDivergence[info.job()];
	SolveReal &localCellCount = parallelCellCount[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		return 0;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
		vit.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::EXTERIOR_LIQUID_CELL ||
			vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
			vitt.getValue() == MaterialLabels::EXTERIOR_BUBBLE_CELL ||
			vitt.getValue() == MaterialLabels::INTERIOR_BUBBLE_CELL)
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

			localAccumulatedDivergence += divergence;
			localMaxDivergence = std::max(localMaxDivergence, divergence);
			++localCellCount;
		    }
		}
	    }
	}

	return 0;
    });
}