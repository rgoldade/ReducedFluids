#include "HDK_AffineFreeSurfacePressureSolver.h"

#include <Eigen/LU>

#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"

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
   IMPLEMENT_DATAFACTORY(HDK_AffineFreeSurfacePressureSolver);
}

// Standard constructor, note that BaseClass was crated by the
// DECLARE_DATAFACTORY and provides an easy way to chain through
// the class hierarchy.
HDK_AffineFreeSurfacePressureSolver::HDK_AffineFreeSurfacePressureSolver(const SIM_DataFactory *factory)
    : BaseClass(factory)
{
}

HDK_AffineFreeSurfacePressureSolver::~HDK_AffineFreeSurfacePressureSolver()
{
}

const SIM_DopDescription* HDK_AffineFreeSurfacePressureSolver::getDopDescription()
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

    static PRM_Name    theUseOldPressureName("useOldPressure", "Use old pressure as an initial guess");
    static PRM_Name    theExtrapolatePressureName("extrapolatePressure", "Extrapolate pressure into interior");

    static PRM_Name 	theDensityName(GAS_NAME_DENSITY, "Liquid Density Field");
    static PRM_Default 	theDensityNameDefault(0, "massdensity");

    static PRM_Name 	theValidFacesName("validFaces", "Valid Faces Field");

    static PRM_Name 	theToleranceName(SIM_NAME_TOLERANCE, "Solver Tolerance");
    static PRM_Default 	theToleranceDefault(1e-5);

    static PRM_Name 	theMaxIterations("maxSolverIterations", "Max Solver Iterations");
    static PRM_Default 	theMaxIterationsDefault(2500);

    static PRM_Name	theActiveCellWidthName("boundaryLayerSize", "Boundary Layer Size");
    static PRM_Default  theActiveCellWidthDefault(2);

    static PRM_Name 	theUseTiledInteriorName("useTiledInterior", "Use Tiled Interior");

    static PRM_Name	theTileSizeName("tileSize", "Interior Tile Size");
    static PRM_Default  theTileSizeDefault(16);

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

	PRM_Template(PRM_TOGGLE, 1, &theExtrapolatePressureName, PRMoneDefaults),

	PRM_Template(PRM_STRING, 1, &theDensityName, &theDensityNameDefault),
        
        PRM_Template(PRM_STRING, 1, &theValidFacesName),

    	PRM_Template(PRM_FLT, 1, &theToleranceName, &theToleranceDefault),
    	PRM_Template(PRM_FLT, 1, &theMaxIterations, &theMaxIterationsDefault),

	PRM_Template(PRM_INT, 1, &theActiveCellWidthName, &theActiveCellWidthDefault),
	
	PRM_Template(PRM_TOGGLE, 1, &theUseTiledInteriorName, PRMoneDefaults),
	PRM_Template(PRM_INT, 1, &theTileSizeName, &theTileSizeDefault),

	PRM_Template(PRM_TOGGLE, 1, &thePrintCellsName, PRMzeroDefaults),
	PRM_Template(PRM_TOGGLE, 1, &theOnlyPrintCellsName, PRMzeroDefaults),

	PRM_Template(PRM_STRING, 1, &thePrintedCellsGeometryName,
				    &thePrintedCellsGeometryDefault),

    	PRM_Template()
    };

    static SIM_DopDescription theDopDescription(true,
						"HDK_AffineFreeSurfacePressureSolver",
						"HDK Affine Free Surface Pressure Solver",
						"$OS",
						classname(),
						theTemplates);

    setGasDescription(theDopDescription);

    return &theDopDescription;
}

bool HDK_AffineFreeSurfacePressureSolver::solveGasSubclass(SIM_Engine &engine,
						    SIM_Object *obj,
						    SIM_Time time,
						    SIM_Time timestep)
{
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
    
    std::array<const SIM_RawField *, 3> cutCellWeights;
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

    // Load pressure

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
    const fpreal dx = velocity->getVoxelSize().maxComponent();

    const SIM_ScalarField *solidSurfaceField = getConstScalarField(obj, GAS_NAME_COLLISION);
    const SIM_RawField *solidSurface;
    SIM_RawField localSolidSurface;

    if (solidSurfaceField == nullptr)
    {
        // Treat as all fluid.
        localSolidSurface.init(SIM_SAMPLE_CENTER,
				velocity->getOrig(),
				velocity->getSize(),
				velocity->getTotalVoxelRes()[0],
				velocity->getTotalVoxelRes()[1],
				velocity->getTotalVoxelRes()[2]);

        // Solid level set convention in Houdini is negative outside and positive inside.
        localSolidSurface.makeConstant(-10. * dx);
        solidSurface = &localSolidSurface;
    }
    else
	solidSurface = solidSurfaceField->getField();

    // Load liquid density

    const SIM_ScalarField *liquidDensityField = getConstScalarField(obj, GAS_NAME_DENSITY);

    if (liquidDensityField == nullptr)
    {
        addError(obj, SIM_MESSAGE, "There is no liquid density to simulate with", UT_ERROR_WARNING);
        return false;
    }

    const SIM_RawField &liquidDensity = *liquidDensityField->getField();

    if (!liquidDensity.isAligned(&liquidSurface))
    {
        addError(obj, SIM_MESSAGE, "Density must align with the surface volume", UT_ERROR_WARNING);
        return false;
    }

    fpreal32 constantLiquidDensity = 0.;
    if (!liquidDensity.field()->isConstant(&constantLiquidDensity))
    {
	addError(obj, SIM_MESSAGE, "Variable density is not currently supported", UT_ERROR_WARNING);
        return false;
    }

    ////////////////////////////////////////////
    //
    // Build material cell labels
    //
    ////////////////////////////////////////////

    std::cout << "//\n//\n// Starting tiled free surface pressure solver\n//\n//" << std::endl;

    SIM_RawIndexField materialCellLabels;
    
    {
    	UT_PerfMonAutoSolveEvent event(this, "Build liquid cell labels");
    	std::cout << "\n// Build liquid cell labels" << std::endl;

	materialCellLabels.match(liquidSurface);
	materialCellLabels.makeConstant(HDK::Utilities::SOLID_CELL);

	HDK::Utilities::buildMaterialCellLabels(materialCellLabels,
						liquidSurface, 
						*solidSurface,
						cutCellWeights);

	// Set labels specific to this solver
	setMaterialCellLabels(materialCellLabels);

	materialCellLabels.fieldNC()->collapseAllTiles();
    }

    ////////////////////////////////////////////
    //
    // Build layer of active cells along liquid-air boundary
    //
    ////////////////////////////////////////////

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

    {
	UT_PerfMonAutoSolveEvent event(this, "Build active cells along liquid boundary");
    	std::cout << "\n// Build active cells along liquid boundary" << std::endl;

	UT_Array<UT_Array<UT_Vector3I>> parallelActiveCellLayer;
	parallelActiveCellLayer.setSize(threadCount);

	// Build list of cells along the liquid boundary
	buildInitialAirBoundaryLayer(parallelActiveCellLayer,
					materialCellLabels,
					cutCellWeights);

	UT_Array<UT_Vector3I> activeCellLayer;

	UT_Array<bool> isTileOccupiedList;
	isTileOccupiedList.setSize(materialCellLabels.field()->numTiles());

	// Now flood inwards one layer at a time, repeating the process
	const int boundaryLayerSize = getBoundaryLayerSize();
	for (int layer = 0; layer < boundaryLayerSize; ++layer)
	{
	    // Combine parallel active cell lists
	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelActiveCellLayer[thread].size();

	    activeCellLayer.clear();
	    activeCellLayer.bumpCapacity(listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
	    {
		activeCellLayer.concat(parallelActiveCellLayer[thread]);
		parallelActiveCellLayer[thread].clear();
	    }

	    UTparallelSort(activeCellLayer.begin(), activeCellLayer.end(), cellCompare);

	    // Build tile list
	    isTileOccupiedList.constant(false);
	    HDK::Utilities::findOccupiedIndexTiles(isTileOccupiedList,
						    activeCellLayer,
						    materialCellLabels);

	    // Uncompress tiles of the new active cells
	    HDK::Utilities::uncompressTiles(materialCellLabels,
					    isTileOccupiedList);
	
	    setActiveLiquidLayerCells(materialCellLabels,
					activeCellLayer);
	    
	    // Build next layer of active cells unless we're at the final layer
	    if (layer < boundaryLayerSize - 1)
	    {
		buildNextLiquidBoundaryLayer(parallelActiveCellLayer,
						activeCellLayer,
						materialCellLabels,
						cutCellWeights);
	    }
	}
    }

    {
	UT_PerfMonAutoSolveEvent event(this, "Build active cells along liquid-solid boundary");
    	std::cout << "\n// Build active cells along liquid-solid boundary" << std::endl;

	UT_Array<UT_Array<UT_Vector3I>> parallelActiveCellLayer;
	parallelActiveCellLayer.setSize(threadCount);

	// Build list of cells along the liquid boundary
	buildInitialSolidBoundaryLayer(parallelActiveCellLayer,
					materialCellLabels,
					cutCellWeights);

	UT_Array<UT_Vector3I> activeCellLayer;

	UT_Array<bool> isTileOccupiedList;
	isTileOccupiedList.setSize(materialCellLabels.field()->numTiles());

	SIM_RawIndexField visitedCells;
	visitedCells.match(materialCellLabels);
	visitedCells.makeConstant(UNVISITED_CELL);

	int solidBoundaryLayerSize = 2;
	for (int layer = 0; layer < solidBoundaryLayerSize; ++layer)
	{
	    // Combine parallel active cell lists
	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelActiveCellLayer[thread].size();

	    activeCellLayer.clear();
	    activeCellLayer.bumpCapacity(listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
	    {
		activeCellLayer.concat(parallelActiveCellLayer[thread]);
		parallelActiveCellLayer[thread].clear();
	    }

	    UTparallelSort(activeCellLayer.begin(), activeCellLayer.end(), cellCompare);

	    // Build tile list
	    isTileOccupiedList.constant(false);
	    HDK::Utilities::findOccupiedIndexTiles(isTileOccupiedList,
						    activeCellLayer,
						    materialCellLabels);

	    // Uncompress tiles of the new active cells
	    HDK::Utilities::uncompressTiles(materialCellLabels,
					    isTileOccupiedList);

	    HDK::Utilities::uncompressTiles(visitedCells,
					    isTileOccupiedList);
	
	    setActiveSolidLayerCells(materialCellLabels, visitedCells, activeCellLayer);

	    // Build next layer of active cells unless we're at the final layer
	    if (layer < solidBoundaryLayerSize - 1)
	    {
		buildNextSolidBoundaryLayer(parallelActiveCellLayer,
					    activeCellLayer,
					    materialCellLabels,
					    visitedCells,
					    cutCellWeights);
	    }
	}
    }

    ////////////////////////////////////////////
    //
    // Build layer of active cells along liquid-solid boundary
    //
    ////////////////////////////////////////////

    if (getUseTiledInterior())
    {
	UT_PerfMonAutoSolveEvent event(this, "Build tile boundary active cells");
    	std::cout << "\n// Build tile boundary active cells" << std::endl;

	int tileSize = getTileSize();
	std::cout << "Tile size: " << tileSize << std::endl;
	buildTiledActiveCells(materialCellLabels,
				tileSize);
    }

    ////////////////////////////////////////////
    //
    // Build unique indices for each interior region
    //
    ////////////////////////////////////////////

    SIM_RawIndexField interiorRegionIndices;
    exint interiorRegionCount = 0;
    {
	UT_PerfMonAutoSolveEvent event(this, "Build interior regions");
    	std::cout << "\n// Build interior regions" << std::endl;

	interiorRegionIndices.match(liquidSurface);

	SIM_VolumetricConnectedComponentBuilder<> interiorRegionBuilder(interiorRegionIndices, materialCellLabels, cutCellWeights.data());
	interiorRegionCount = interiorRegionBuilder.buildConnectedComponents([](const exint label) { return label == MaterialLabels::INTERIOR_LIQUID_CELL; });

	HDK::Utilities::overwriteIndices(interiorRegionIndices,
					    SIM_VolumetricConnectedComponentBuilder<>::INACTIVE_REGION,
					    HDK::Utilities::UNLABELLED_CELL);
    }

    ////////////////////////////////////////////
    //
    // Build bounding box for each interior region
    //
    ////////////////////////////////////////////

    exint tbbGrainSize = interiorRegionCount / (4 * threadCount);

    {
	UT_PerfMonAutoSolveEvent event(this, "Build interior region bounding boxes");
	std::cout << "\n// Build interior region bounding boxes" << std::endl;

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

	if (interiorRegionCount > regionMap)
	{
	    interiorRegionCount = regionMap;

	    remapInteriorRegions(interiorRegionIndices,
				    materialCellLabels,
				    doRemoveRegion,
				    remapRegion);
	}
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

    std::cout << "    Active cell count: " << activeCellCount << std::endl;
    std::cout << "    Interior region count: " << interiorRegionCount << std::endl;

    ////////////////////////////////////////////
    //
    // Print tiled geometry for display and debugging
    //
    ////////////////////////////////////////////

    if (getPrintCells())
    {
	using SIM::FieldUtils::getFieldValue;
	using SIM::FieldUtils::setFieldValue;

	std::cout << "// Build cell geometry" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Printing cell geometry");

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

	using PointPair = std::pair<UT_Vector3, int>;
	std::vector<PointPair> pointPairs;

	{
	    tbb::enumerable_thread_specific<std::vector<PointPair>> parallelPointPairs;

	    tbb::parallel_for(tbb::blocked_range<int>(0, materialCellLabels.field()->numTiles()), [&](const tbb::blocked_range<int> &range)
	    {
		auto &localPointPairs = parallelPointPairs.local();

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
			vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
			vit.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
		    {
			for (; !vit.atEnd(); vit.advance())
			{
			    if (vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
				vit.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
			    {
				UT_Vector3I cell(vit.x(), vit.y(), vit.z());

				int activityLabel;

				if (vit.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
				{
				    assert(getFieldValue(activeCellIndices, cell) >= 0);
				    assert(getFieldValue(interiorRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
				    activityLabel = 0;
				}
				else
				{
				    assert(getFieldValue(activeCellIndices, cell) == HDK::Utilities::UNLABELLED_CELL);
				    assert(getFieldValue(interiorRegionIndices, cell) >= 0);
				    activityLabel = 1;
				}

				UT_Vector3 cellCenter = velocity->getOrig() + UT_Vector3(SolveReal(cell[0]) + .5, SolveReal(cell[1]) + .5, SolveReal(cell[2]) + .5) * dx;

				localPointPairs.emplace_back(cellCenter, activityLabel);
			    }
			}
		    }
		}
	    });

	    exint listSize = 0;

	    parallelPointPairs.combine_each([&](const std::vector<PointPair> &localPointPairs)
	    {
		listSize += localPointPairs.size();
	    });

	    pointPairs.reserve(listSize);

	    parallelPointPairs.combine_each([&](const std::vector<PointPair> &localPointPairs)
	    {
		pointPairs.insert(pointPairs.end(), localPointPairs.begin(), localPointPairs.end());
	    });	    
	}

	// Add points to details
	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());

	for (const auto &pointPair : pointPairs)
	{
	    GA_Offset pointOffset = cellActivityDetail->appendPoint();
	    cellSizeHandle.set(pointOffset, dx);

	    cellActivityHandle.set(pointOffset, pointPair.second);
	    cellActivityDetail->setPos3(pointOffset, pointPair.first);
	}

	cellSizeHandle.bumpDataId();
	cellActivityHandle.bumpDataId();
	cellActivityDetail->getAttributes().bumpAllDataIds(GA_ATTRIB_POINT);

        if (getOnlyPrintCells())
            return true;
    }

    ////////////////////////////////////////////
    //
    // Build COM for each interior region
    //
    ////////////////////////////////////////////

    UT_Array<UT_Vector3T<SolveReal>> interiorRegionCOM;
    interiorRegionCOM.setSize(interiorRegionCount);
    interiorRegionCOM.constant(UT_Vector3T<SolveReal>(0,0,0));

    {
	UT_PerfMonAutoSolveEvent event(this, "Build interior region COM");
	std::cout << "\n// Build interior region COM" << std::endl;

	UT_Array<UT_Array<UT_Vector3T<SolveReal>>> parallelInteriorRegionCOM;
	parallelInteriorRegionCOM.setSize(threadCount);

	UT_Array<UT_Array<SolveReal>> parallelInteriorRegionCellCount;
	parallelInteriorRegionCellCount.setSize(threadCount);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    parallelInteriorRegionCOM[thread].setSize(interiorRegionCount);
	    parallelInteriorRegionCOM[thread].constant(UT_Vector3T<SolveReal>(0,0,0));

	    parallelInteriorRegionCellCount[thread].setSize(interiorRegionCount);
	    parallelInteriorRegionCellCount[thread].constant(0);
	}

	buildInteriorRegionCOM(parallelInteriorRegionCOM,
				parallelInteriorRegionCellCount,
				materialCellLabels,
				interiorRegionIndices,
				activeCellIndices);

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

    UT_Array<UT_Vector3T<SolveReal>> interiorBestFitLinearVelocities;
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
    // Build rhs using fluid velocity and affine interior field
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
    // Build sparse matrix consisting of a standard poisson system and affine interior regions
    //
    ////////////////////////////////////////////

    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> poissonMatrix(activeCellCount, activeCellCount);
    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> diagonalPrecondMatrix(activeCellCount, activeCellCount);

    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> affineCollectionMatrix(11 * interiorRegionCount, activeCellCount);

    {
	std::cout << "\n// Build poisson and affine transfer systems" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Build poisson and affine transfer systems");

	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelPoissonElements(threadCount);
	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelDiagonalPrecondElements(threadCount);
	std::vector<std::vector<Eigen::Triplet<SolveReal>>> parallelAffineCollectionElements(threadCount);

	buildLinearSystem(parallelPoissonElements,
			    parallelDiagonalPrecondElements,
			    parallelAffineCollectionElements,
			    materialCellLabels,
			    interiorRegionIndices,
			    activeCellIndices,
			    cutCellWeights,
			    liquidSurface,
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

	// Compile affine elements
	{
	    exint listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelAffineCollectionElements[thread].size();

	    std::vector<Eigen::Triplet<SolveReal>> affineCollectionElements;
	    affineCollectionElements.reserve(listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
		affineCollectionElements.insert(affineCollectionElements.end(), parallelAffineCollectionElements[thread].begin(), parallelAffineCollectionElements[thread].end());

	    affineCollectionMatrix.setFromTriplets(affineCollectionElements.begin(), affineCollectionElements.end());
	    affineCollectionMatrix.makeCompressed();
	}
    }
    
    Eigen::SparseMatrix<SolveReal, Eigen::RowMajor> affineMassMatrix(11 * interiorRegionCount, 11 * interiorRegionCount);

    {
	std::cout << "\n// Build affine mass matrix" << std::endl;
	UT_PerfMonAutoSolveEvent event(this, "Build affine mass matrix");

	UT_Array<UT_Array<MassMatrix>> parallelAffineMassMatrices;
	parallelAffineMassMatrices.setSize(threadCount);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    parallelAffineMassMatrices[thread].setSize(interiorRegionCount);
	    parallelAffineMassMatrices[thread].constant(MassMatrix::Zero());
	}

	buildAffineMassMatrixSystems(parallelAffineMassMatrices,
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
		    interiorRegionMassMatrices[interiorRegion] += parallelAffineMassMatrices[thread][interiorRegion];
	    }
	});

	std::vector<Eigen::Triplet<SolveReal>> affineMassMatrixElements;
	affineMassMatrixElements.reserve(11 * 11 * interiorRegionCount);

	// Build combined mass matrix system
	{
	    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<SolveReal>>> parallelMassMatrixElements;
    	    tbb::parallel_for(tbb::blocked_range<exint>(0, interiorRegionCount, tbbGrainSize), [&](const tbb::blocked_range<exint> &range)
	    {
		auto &localMassMatrixElements = parallelMassMatrixElements.local();
		for (exint interiorRegion = range.begin(); interiorRegion != range.end(); ++interiorRegion)
		{
		    // Invert mass matrix
		    MassMatrix invertedMassMatrix = interiorRegionMassMatrices[interiorRegion].inverse();

		    for (exint i = 0; i < 11; ++i)
			for (exint j = 0; j < 11; ++j)
			{
			    if (invertedMassMatrix(i,j) != 0)
				localMassMatrixElements.emplace_back(11 * interiorRegion + i, 11 * interiorRegion + j, invertedMassMatrix(i,j));
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

	{
	    std::cout << "// Solve affine tiled linear system" << std::endl;
		UT_PerfMonAutoSolveEvent event(this, "Solve affine tiled linear system");

	    Vector tempVector0 = Vector::Zero(11 * interiorRegionCount);
	    Vector tempVector1 = Vector::Zero(11 * interiorRegionCount);

	    auto MatrixVectorMultiply = [&](Vector &destinationVector, const Vector &sourceVector)
	    {
		// UT_StopWatch multiplyTimer;
		// multiplyTimer.start();

		destinationVector.noalias() = poissonMatrix * sourceVector;

		// auto time = multiplyTimer.stop();
		// std::cout << "\n    Poisson matrix-vector multiply time: " << time << std::endl;
		// multiplyTimer.clear();
		// multiplyTimer.start();

		if (interiorRegionCount > 0)
		{
		    tempVector0.noalias() = affineCollectionMatrix  * sourceVector;

		    tempVector1.noalias() = affineMassMatrix * tempVector0;

		    distributeAffineVectors(destinationVector,
					    tempVector1,
					    interiorBoundaryCells,
					    materialCellLabels,
					    interiorRegionIndices,
					    activeCellIndices,
					    interiorRegionCOM,
					    dx);

		    //destinationVector.noalias() += affineDistributeMatrix * tempVector1;
		}		

		// time = multiplyTimer.stop();
		// std::cout << "    Affine matrix-vector multiply time: " << time << std::endl;

	    };

	    auto DiagonalPreconditioner = [&](Vector &destinationVector, const Vector &sourceVector)
	    {
		destinationVector.noalias() = diagonalPrecondMatrix * sourceVector;
	    };


	    HDK::Utilities::solveConjugateGradient(solutionVector,
						    rhsVector,
						    MatrixVectorMultiply,
						    DiagonalPreconditioner,
						    getSolverTolerance(),
						    getMaxSolverIterations(),
						    false,
						    false);

	    Vector residualVector = Vector::Zero(activeCellCount);
	    MatrixVectorMultiply(residualVector, solutionVector);
	    residualVector = rhsVector - residualVector;

	    std::cout << "    Re-computed residual: " << std::sqrt(residualVector.squaredNorm() / rhsVector.squaredNorm()) << std::endl;
	    std::cout << "    L-infinity residual: " << residualVector.lpNorm<Eigen::Infinity>() << std::endl;
	}
    }

    ////////////////////////////////////////////
    //
    // Copy pressure to grid
    //
    ////////////////////////////////////////////

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
    // Extrapolate pressure into interior regions
    //
    ////////////////////////////////////////////

    if (getDoExtrapolatePressure())
    {
	UT_PerfMonAutoSolveEvent event(this, "Extrapolate pressure into interior regions");
	std::cout << "// Extrapolate pressure into interior regions" << std::endl;

	UT_Array<bool> isTileOccupiedList;
	isTileOccupiedList.setSize(materialCellLabels.field()->numTiles());

	SIM_RawIndexField visitedCells;
	visitedCells.match(materialCellLabels);
	visitedCells.makeConstant(UNVISITED_CELL);

	// Build initial list of interior cells with affine neighbours
	UT_Array<UT_Array<UT_Vector3I>> parallelExtrapolationCellLayer;

	parallelExtrapolationCellLayer.setSize(threadCount);

	buildInitialExtrapolationLayer(parallelExtrapolationCellLayer,
					visitedCells,
					materialCellLabels,
					cutCellWeights);

	// Compile list and sort to handle duplicates
	UT_Array<UT_Vector3I> extrapolationCellLayer;

	exint listSize = 0;
	for (int thread = 0; thread < threadCount; ++thread)
	    listSize += parallelExtrapolationCellLayer[thread].size();

	extrapolationCellLayer.clear();
	extrapolationCellLayer.bumpCapacity(listSize);

	for (int thread = 0; thread < threadCount; ++thread)
	{
	    extrapolationCellLayer.concat(parallelExtrapolationCellLayer[thread]);
	    parallelExtrapolationCellLayer[thread].clear();
	}

	while (extrapolationCellLayer.size() > 0)
	{
	    // Sort list to handle duplicates
	    UTparallelSort(extrapolationCellLayer.begin(), extrapolationCellLayer.end(), cellCompare);

	    // Extrapolate from adjacent visited cells
	    extrapolatePressure(*pressure,
				materialCellLabels,
				visitedCells,
				cutCellWeights,
				extrapolationCellLayer);

	    // Build tile list
	    isTileOccupiedList.constant(false);
	    HDK::Utilities::findOccupiedIndexTiles(isTileOccupiedList,
						    extrapolationCellLayer,
						    materialCellLabels);

	    // Uncompress tiles of the new active cells
	    HDK::Utilities::uncompressTiles(*pressure,
					    isTileOccupiedList);

	    // Uncompress tiles of the new active cells
	    HDK::Utilities::uncompressTiles(visitedCells,
					    isTileOccupiedList);

	    setVisitedExtrapolationCells(visitedCells, materialCellLabels, extrapolationCellLayer);

	    buildNextExtrapolationLayer(parallelExtrapolationCellLayer,
					extrapolationCellLayer,
					materialCellLabels,
					visitedCells,
					cutCellWeights);

	    listSize = 0;
	    for (int thread = 0; thread < threadCount; ++thread)
		listSize += parallelExtrapolationCellLayer[thread].size();

	    extrapolationCellLayer.clear();
	    extrapolationCellLayer.bumpCapacity(listSize);

	    for (int thread = 0; thread < threadCount; ++thread)
	    {
		extrapolationCellLayer.concat(parallelExtrapolationCellLayer[thread]);
		parallelExtrapolationCellLayer[thread].clear();
	    }
	}
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
HDK_AffineFreeSurfacePressureSolver::setMaterialCellLabels(SIM_RawIndexField &materialCellLabels) const
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
		    cellLabel = MaterialLabels::AIR_CELL;
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
			vit.setValue(MaterialLabels::AIR_CELL);
		    else
		    {
			assert(label == HDK::Utilities::SOLID_CELL);
			vit.setValue(MaterialLabels::SOLID_CELL);
		    }
		}
	    }
	}
    });
}

void
HDK_AffineFreeSurfacePressureSolver::buildInitialAirBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelActiveCellLayer,
							    const SIM_RawIndexField &materialCellLabels,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInitialActiveCellLayerAlgorithm;
    buildInitialActiveCellLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_Array<UT_Vector3I> &localActiveCellLayer = parallelActiveCellLayer[info.job()];

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		    break;

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
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
				    getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::AIR_CELL)
				    isBoundaryCell = true;
			    }

			if (isBoundaryCell)
			    localActiveCellLayer.append(cell);
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineFreeSurfacePressureSolver::setActiveLiquidLayerCells(SIM_RawIndexField &materialCellLabels,
								const UT_Array<UT_Vector3I> &activeCellLayer) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    // Assign active cells
    UTparallelForLightItems(UT_BlockedRange<exint>(0, activeCellLayer.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;

	exint startIndex = range.begin();
    
	if (startIndex > 0)
	{
	    while (startIndex != range.end() && activeCellLayer[startIndex] == activeCellLayer[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	for (exint cellIndex = startIndex; cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = activeCellLayer[cellIndex];

	    if (cell == oldCell)
		continue;
	
	    oldCell = cell;

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL);
	    setFieldValue(materialCellLabels, cell, MaterialLabels::ACTIVE_LIQUID_CELL);
	}
    });

}

void
HDK_AffineFreeSurfacePressureSolver::buildNextLiquidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelNewActiveCellLayer,
								    const UT_Array<UT_Vector3I> &oldActiveCellLayer,
								    const SIM_RawIndexField &materialCellLabels,
								    const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildNextActiveCellLayerAlgorithm;
    buildNextActiveCellLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	if (boss->opInterrupt())
	    return 0;

	UT_Array<UT_Vector3I> &localNewActiveCellLayer = parallelNewActiveCellLayer[info.job()];

	exint startIndex, endIndex;
	const exint elementSize = oldActiveCellLayer.entries();
	info.divideWork(elementSize, startIndex, endIndex);

	if (startIndex > 0)
	{
	    while (startIndex != endIndex && oldActiveCellLayer[startIndex] == oldActiveCellLayer[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	const exint localEndIndex = endIndex;
	for (exint cellIndex = startIndex; cellIndex < localEndIndex; ++cellIndex)
	{
	    UT_Vector3I cell = oldActiveCellLayer[cellIndex];

	    if (cell == oldCell)
		continue;

	    oldCell = cell;

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::ACTIVE_LIQUID_CELL);

	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

		    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
			continue;

		    UT_Vector3I face = cellToFaceMap(cell, axis, direction);

		    if (getFieldValue(*cutCellWeights[axis], face) > 0 &&
			getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL)
			localNewActiveCellLayer.append(adjacentCell);
		}
	}

	return 0;
    });
}

void
HDK_AffineFreeSurfacePressureSolver::buildInitialSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelActiveCellLayer,
								const SIM_RawIndexField &materialCellLabels,
								const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInitialActiveCellLayerAlgorithm;
    buildInitialActiveCellLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_Array<UT_Vector3I> &localActiveCellLayer = parallelActiveCellLayer[info.job()];

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		    break;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
			vitt.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
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
			    localActiveCellLayer.append(cell);
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineFreeSurfacePressureSolver::setActiveSolidLayerCells(SIM_RawIndexField &materialCellLabels,
								SIM_RawIndexField &visitedCells,
								const UT_Array<UT_Vector3I> &activeCellLayer) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    // Assign active cells
    UTparallelForLightItems(UT_BlockedRange<exint>(0, activeCellLayer.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;

	exint startIndex = range.begin();
    
	if (startIndex > 0)
	{
	    while (startIndex != range.end() && activeCellLayer[startIndex] == activeCellLayer[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	for (exint cellIndex = startIndex; cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = activeCellLayer[cellIndex];

	    if (cell == oldCell)
		continue;
	
	    oldCell = cell;

	    assert(getFieldValue(visitedCells, cell) == UNVISITED_CELL);
	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL ||
		    getFieldValue(materialCellLabels, cell) == MaterialLabels::ACTIVE_LIQUID_CELL);

	    setFieldValue(visitedCells, cell, VISITED_CELL);
	    setFieldValue(materialCellLabels, cell, MaterialLabels::ACTIVE_LIQUID_CELL);
	}
    });
}

void
HDK_AffineFreeSurfacePressureSolver::buildNextSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelNewActiveCellLayer,
							    const UT_Array<UT_Vector3I> &oldActiveCellLayer,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &visitedCells,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildNextActiveCellLayerAlgorithm;
    buildNextActiveCellLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	if (boss->opInterrupt())
	    return 0;

	UT_Array<UT_Vector3I> &localNewActiveCellLayer = parallelNewActiveCellLayer[info.job()];

	exint startIndex, endIndex;
	const exint elementSize = oldActiveCellLayer.entries();
	info.divideWork(elementSize, startIndex, endIndex);

	if (startIndex > 0)
	{
	    while (startIndex != endIndex && oldActiveCellLayer[startIndex] == oldActiveCellLayer[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	const exint localEndIndex = endIndex;
	for (exint cellIndex = startIndex; cellIndex < localEndIndex; ++cellIndex)
	{
	    UT_Vector3I cell = oldActiveCellLayer[cellIndex];

	    if (cell == oldCell)
		continue;

	    oldCell = cell;

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::ACTIVE_LIQUID_CELL);
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
			    (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::ACTIVE_LIQUID_CELL ||
			    getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL))
			    localNewActiveCellLayer.append(adjacentCell);
		    }
		}
	}

	return 0;
    });
}

void
HDK_AffineFreeSurfacePressureSolver::buildTiledActiveCells(SIM_RawIndexField &materialCellLabels,
							const int tileSize) const
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

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());

			if (cell[0] % tileSize == 0 ||
			    cell[1] % tileSize == 0 ||
			    cell[2] % tileSize == 0)
				setFieldValue(materialCellLabels, cell, MaterialLabels::ACTIVE_LIQUID_CELL);
		    }
		}
	    }
	}
    });
}

void
HDK_AffineFreeSurfacePressureSolver::buildInteriorBoundingBoxes(UT_Array<UT_Array<UT_Vector3I>> &parallelInteriorRegionBBMin,
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
		break;

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
			exint interiorRegion = getFieldValue(interiorRegionIndices, cell);

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
HDK_AffineFreeSurfacePressureSolver::remapInteriorRegions(SIM_RawIndexField &interiorRegionIndices,
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

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
	    {
		for (; !vit.atEnd(); vit.advance())
		{
		    if (vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
		    {
			UT_Vector3I cell(vit.x(), vit.y(), vit.z());

			exint interiorRegion = getFieldValue(interiorRegionIndices, cell);
			assert(interiorRegion >= 0);

			if (doRemoveRegion[interiorRegion])
			{
			    assert(remapRegion[interiorRegion] == -1);
			    setFieldValue(materialCellLabels, cell, MaterialLabels::ACTIVE_LIQUID_CELL);
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
HDK_AffineFreeSurfacePressureSolver::buildActiveCellIndices(SIM_RawIndexField &activeCellIndices,
							const SIM_RawIndexField &materialCellLabels) const
{
    using SIM::FieldUtils::setFieldValue;

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(materialCellLabels.field());

    UT_VoxelTileIteratorI vitt;

    exint activeCellCount = 0;

    // Build active cell indices
    UT_Interrupt *boss = UTgetInterrupt();
    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
	if (boss->opInterrupt())
	    break;

	if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
	{
	    vitt.setTile(vit);

	    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
	    {
		if (vitt.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
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
HDK_AffineFreeSurfacePressureSolver::buildInteriorRegionCOM(UT_Array<UT_Array<UT_Vector3T<SolveReal>>> &parallelInteriorRegionCOM,
							UT_Array<UT_Array<SolveReal>> &parallelInteriorRegionCellCount,
							const SIM_RawIndexField &materialCellLabels,
							const SIM_RawIndexField &interiorRegionIndices,
							const SIM_RawIndexField &activeCellIndices) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();
  
    UT_ThreadedAlgorithm buildInteriorRegionCOMAlgorithm;
    buildInteriorRegionCOMAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_Array<UT_Vector3T<SolveReal>> &localInteriorRegionCOM = parallelInteriorRegionCOM[info.job()];
	UT_Array<SolveReal> &localInteriorRegionCellCount = parallelInteriorRegionCellCount[info.job()];

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		break;

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
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
HDK_AffineFreeSurfacePressureSolver::buildBestFitSystems(UT_Array<UT_Array<MassMatrix>> &parallelBestFitMatrix,
						    UT_Array<UT_Array<ColumnVector>> &parallelBestFitRHS,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_RawIndexField &interiorRegionIndices,
						    const SIM_VectorField &velocity,
						    const UT_Array<UT_Vector3T<SolveReal>> &interiorRegionCOM,
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
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL);

			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);
				assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < interiorRegionIndices.getVoxelRes()[axis]);

				if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::ACTIVE_LIQUID_CELL)
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
		    else assert(getFieldValue(materialCellLabels, cell) != MaterialLabels::INTERIOR_LIQUID_CELL);
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineFreeSurfacePressureSolver::buildRHS(Vector &rhsVector,
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
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::ACTIVE_LIQUID_CELL);
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
					getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL)
				    {
					exint interiorRegion = getFieldValue(interiorRegionIndices, adjacentCell);
					assert(interiorRegion >= 0);

					assert(getFieldValue(activeCellIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);

					UT_Vector3SR offset(cell);
					
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
		    else assert(getFieldValue(materialCellLabels, cell) != MaterialLabels::ACTIVE_LIQUID_CELL);
		}
	    }
	}
    });
}

void
HDK_AffineFreeSurfacePressureSolver::buildLinearSystem(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelPoissonElements,
						    std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelDiagonalPrecondElements,
						    std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelAffineCollectionElements,
						    const SIM_RawIndexField &materialCellLabels,
						    const SIM_RawIndexField &interiorRegionIndices,
						    const SIM_RawIndexField &activeCellIndices,
						    const std::array<const SIM_RawField *, 3> &cutCellWeights,
						    const SIM_RawField &liquidSurface,
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
	std::vector<Eigen::Triplet<SolveReal>> &localAffineCollectionElements = parallelAffineCollectionElements[info.job()];

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		break;

	    if (!vit.isTileConstant())
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    exint index = vitt.getValue();
		    if (index >= 0)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::ACTIVE_LIQUID_CELL);
			assert(getFieldValue(interiorRegionIndices, cell) == HDK::Utilities::UNLABELLED_CELL);

			SolveReal diagonal = 0;

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

				    auto adjacentCellLabel = getFieldValue(materialCellLabels, adjacentCell);
				    if (adjacentCellLabel == MaterialLabels::AIR_CELL)
				    {
					assert(getFieldValue(interiorRegionIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);
					assert(getFieldValue(activeCellIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);

					SolveReal phi0 = getFieldValue(liquidSurface, cell);
					SolveReal phi1 = getFieldValue(liquidSurface, adjacentCell);

					assert(phi1 > 0);

					SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);

					theta = SYSclamp(theta, SolveReal(.01), SolveReal(1.));
					diagonal += weight / theta;
				    }
				    else if (adjacentCellLabel == MaterialLabels::ACTIVE_LIQUID_CELL)
				    {
					assert(getFieldValue(interiorRegionIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);

					exint adjacentIndex = getFieldValue(activeCellIndices, adjacentCell);
					assert(adjacentIndex >= 0);

					localPoissonElements.push_back(Eigen::Triplet<SolveReal>(index, adjacentIndex, -weight));
					diagonal += weight;
				    }
				    else
				    {
					exint interiorRegion = getFieldValue(interiorRegionIndices, adjacentCell);
					assert(interiorRegion >= 0);
					assert(weight == 1);

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
					localAffineCollectionElements.push_back(Eigen::Triplet<SolveReal>(11 * interiorRegion + axis, index, sign));

					// Account for affine velocity offset
					localAffineCollectionElements.push_back(Eigen::Triplet<SolveReal>(11 * interiorRegion + 3 + 3 * axis, index, sign * offset[0]));
					localAffineCollectionElements.push_back(Eigen::Triplet<SolveReal>(11 * interiorRegion + 3 + 3 * axis + 1, index, sign * offset[1]));

					// Account for C_33 = -(C_11 + C_22)
					if (axis == 2)
					{
					    localAffineCollectionElements.push_back(Eigen::Triplet<SolveReal>(11 * interiorRegion + 3, index, -sign * offset[2]));
					    localAffineCollectionElements.push_back(Eigen::Triplet<SolveReal>(11 * interiorRegion + 7, index, -sign * offset[2]));
					}
					else
					    localAffineCollectionElements.push_back(Eigen::Triplet<SolveReal>(11 * interiorRegion + 3 + 3 * axis + 2, index, sign * offset[2]));
				    }
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
HDK_AffineFreeSurfacePressureSolver::buildAffineMassMatrixSystems(UT_Array<UT_Array<MassMatrix>> &parallelAffineMassMatrices,
								    const SIM_RawIndexField &materialCellLabels,
								    const SIM_RawIndexField &interiorRegionIndices,
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

	UT_Array<MassMatrix> &localMassMatrices = parallelAffineMassMatrices[info.job()];

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
			assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL);

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

				    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::ACTIVE_LIQUID_CELL)
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

				    localMassMatrices[interiorRegion] += columnVector * columnVector.transpose();
				}
			    }
		    }
		    else assert(getFieldValue(materialCellLabels, cell) != MaterialLabels::INTERIOR_LIQUID_CELL);
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineFreeSurfacePressureSolver::buildInteriorBoundaryCells(UT_Array<UT_Array<UT_Vector3I>> &parallelInteriorBoundaryCells,
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

	    if (!vit.isTileConstant() || vit.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
		    {
			UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

			bool isBoundaryCell = false;
			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				if (adjacentCell[axis] < 0 || adjacentCell[axis] >= materialCellLabels.getVoxelRes()[axis])
				    continue;

				if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL)
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
HDK_AffineFreeSurfacePressureSolver::distributeAffineVectors(Vector &destinationVector,
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
	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::ACTIVE_LIQUID_CELL);

	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

		    assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < materialCellLabels.getVoxelRes()[axis]);

		    if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL)
		    {
			exint interiorRegion = getFieldValue(interiorRegionIndices, adjacentCell);
			assert(interiorRegion >= 0);

			assert(getFieldValue(activeCellIndices, adjacentCell) == HDK::Utilities::UNLABELLED_CELL);

			UT_Vector3SR offset(cell);

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

void
HDK_AffineFreeSurfacePressureSolver::buildValidFaces(SIM_VectorField &validFaces,
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
						[](const exint label){ return label == MaterialLabels::ACTIVE_LIQUID_CELL || label == MaterialLabels::INTERIOR_LIQUID_CELL; },
						axis);

	HDK::Utilities::uncompressTiles(*validFaces.getField(axis), isTileOccupiedList);

	HDK::Utilities::classifyValidFaces(*validFaces.getField(axis),
					    materialCellLabels,
					    *cutCellWeights[axis],
					    [](const exint label){ return label == MaterialLabels::ACTIVE_LIQUID_CELL || label == MaterialLabels::INTERIOR_LIQUID_CELL; },
					    axis);
    }
}

void
HDK_AffineFreeSurfacePressureSolver::copyPressureVectorToGrid(SIM_RawField &pressure,
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
HDK_AffineFreeSurfacePressureSolver::buildInitialExtrapolationLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExtrapolationCellLayer,
								    SIM_RawIndexField &visitedCells,
								    const SIM_RawIndexField &materialCellLabels,
								    const std::array<const SIM_RawField *, 3> &cutCellWeights) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildInitialExtrapolationCellLayerAlgorithm;
    buildInitialExtrapolationCellLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	UT_Array<UT_Vector3I> &localExtrapolationCellLayer = parallelExtrapolationCellLayer[info.job()];

	UT_VoxelArrayIteratorI vit;
	vit.setConstArray(materialCellLabels.field());
	vit.splitByTile(info);

	UT_VoxelTileIteratorI vitt;

	for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
	{
	    if (boss->opInterrupt())
		    break;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

		    if (vitt.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
		    {
			setFieldValue(visitedCells, cell, VISITED_CELL);
		    }
		    else if (vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL)
		    {
			bool isBoundaryCell = false;
			for (int axis : {0,1,2})
			    for (int direction : {0,1})
			    {
				UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

				assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < materialCellLabels.getVoxelRes()[axis]);
				assert(getFieldValue(*cutCellWeights[axis], cellToFaceMap(cell, axis, direction)) == 1);

				if (getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::ACTIVE_LIQUID_CELL)
				    isBoundaryCell = true;
			    }

			if (isBoundaryCell)
			    localExtrapolationCellLayer.append(cell);
		    }
		}
	    }
	}

	return 0;
    });
}

void
HDK_AffineFreeSurfacePressureSolver::extrapolatePressure(SIM_RawField &pressure,
							    const SIM_RawIndexField &materialCellLabels,
							    const SIM_RawIndexField &visitedCells,
							    const std::array<const SIM_RawField *, 3> &cutCellWeights,
							    const UT_Array<UT_Vector3I> &extrapolationCellLayer) const
{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    // Assign active cells
    UTparallelForLightItems(UT_BlockedRange<exint>(0, extrapolationCellLayer.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;

	exint startIndex = range.begin();
    
	if (startIndex > 0)
	{
	    while (startIndex != range.end() && extrapolationCellLayer[startIndex] == extrapolationCellLayer[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	for (exint cellIndex = startIndex; cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = extrapolationCellLayer[cellIndex];

	    if (cell == oldCell)
		continue;
	
	    oldCell = cell;

	    assert(getFieldValue(visitedCells, cell) == UNVISITED_CELL);
	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL);

	    fpreal accumulatedPressure = 0;
	    fpreal cellCount = 0;

	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

		    assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < materialCellLabels.getVoxelRes()[axis]);
		    assert(getFieldValue(*cutCellWeights[axis], cellToFaceMap(cell, axis, direction)) == 1);
		    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL ||
			    getFieldValue(materialCellLabels, cell) == MaterialLabels::ACTIVE_LIQUID_CELL);

		    if (getFieldValue(visitedCells, adjacentCell) == VISITED_CELL)
		    {
			accumulatedPressure += getFieldValue(pressure, adjacentCell);
			++cellCount;
		    }
		}

	    assert(cellCount > 0);

	    setFieldValue(pressure, cell, accumulatedPressure / cellCount);
	}
    });
}

void
HDK_AffineFreeSurfacePressureSolver::setVisitedExtrapolationCells(SIM_RawIndexField &visitedCells,
								    const SIM_RawIndexField &materialCellLabels,
								    const UT_Array<UT_Vector3I> &extrapolationCellLayer) const
{
    using SIM::FieldUtils::getFieldValue;
    using SIM::FieldUtils::setFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    // Assign active cells
    UTparallelForLightItems(UT_BlockedRange<exint>(0, extrapolationCellLayer.size()), [&](const UT_BlockedRange<exint> &range)
    {
	if (boss->opInterrupt())
	    return;

	exint startIndex = range.begin();
    
	if (startIndex > 0)
	{
	    while (startIndex != range.end() && extrapolationCellLayer[startIndex] == extrapolationCellLayer[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	for (exint cellIndex = startIndex; cellIndex != range.end(); ++cellIndex)
	{
	    UT_Vector3I cell = extrapolationCellLayer[cellIndex];

	    if (cell == oldCell)
		continue;
	
	    oldCell = cell;

	    assert(getFieldValue(visitedCells, cell) == UNVISITED_CELL);
	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL);

	    setFieldValue(visitedCells, cell, VISITED_CELL);
	}
    });
}

void
HDK_AffineFreeSurfacePressureSolver::buildNextExtrapolationLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExtrapolationCellLayer,
								    const UT_Array<UT_Vector3I> &oldExtrapolationCellLayer,
								    const SIM_RawIndexField &materialCellLabels,
								    const SIM_RawIndexField &visitedCells,
								    const std::array<const SIM_RawField *, 3> &cutCellWeights) const

{
    using SIM::FieldUtils::cellToCellMap;
    using SIM::FieldUtils::cellToFaceMap;
    using SIM::FieldUtils::getFieldValue;

    UT_Interrupt *boss = UTgetInterrupt();

    UT_ThreadedAlgorithm buildNextExtrapolationLayerAlgorithm;
    buildNextExtrapolationLayerAlgorithm.run([&](const UT_JobInfo &info)
    {
	if (boss->opInterrupt())
	    return 0;

	UT_Array<UT_Vector3I> &localNewExtrapolationLayer = parallelExtrapolationCellLayer[info.job()];

	exint startIndex, endIndex;
	const exint elementSize = oldExtrapolationCellLayer.entries();
	info.divideWork(elementSize, startIndex, endIndex);

	if (startIndex > 0)
	{
	    while (startIndex != endIndex && oldExtrapolationCellLayer[startIndex] == oldExtrapolationCellLayer[startIndex - 1])
		++startIndex;
	}

	UT_Vector3I oldCell(-1,-1,-1);

	const exint localEndIndex = endIndex;
	for (exint cellIndex = startIndex; cellIndex < localEndIndex; ++cellIndex)
	{
	    UT_Vector3I cell = oldExtrapolationCellLayer[cellIndex];

	    if (cell == oldCell)
		continue;

	    oldCell = cell;

	    assert(getFieldValue(materialCellLabels, cell) == MaterialLabels::INTERIOR_LIQUID_CELL);
	    assert(getFieldValue(visitedCells, cell) == VISITED_CELL);

	    for (int axis : {0,1,2})
		for (int direction : {0,1})
		{
		    UT_Vector3I adjacentCell = cellToCellMap(cell, axis, direction);

		    assert(adjacentCell[axis] >= 0 && adjacentCell[axis] < materialCellLabels.getVoxelRes()[axis]);
		    assert(getFieldValue(*cutCellWeights[axis], cellToFaceMap(cell, axis, direction)) == 1);

		    if (getFieldValue(visitedCells, adjacentCell) == UNVISITED_CELL)
		    {
			assert(getFieldValue(materialCellLabels, adjacentCell) == MaterialLabels::INTERIOR_LIQUID_CELL);
			localNewExtrapolationLayer.append(adjacentCell);
		    }
		}
	}

	return 0;
    });
}

void
HDK_AffineFreeSurfacePressureSolver::applySolutionToVelocity(SIM_RawField &velocity,
								const SIM_RawField &validFaces,
								const SIM_RawIndexField &materialCellLabels,
								const SIM_RawIndexField &interiorRegionIndices,
								const SIM_RawIndexField &activeCellIndices,
								const SIM_RawField &pressure,
								const SIM_RawField &liquidSurface,
								const SIM_RawField &cutCellWeights,
								const UT_Array<UT_Vector3T<SolveReal>> &interiorRegionCOM,
								const UT_Array<UT_Vector3T<SolveReal>> &interiorLinearVelocities,
								const UT_Array<GradientMatrix> &interiorVelocityGradients,
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

			assert(getFieldValue(cutCellWeights, face) > 0);

			UT_Vector3I backwardCell = faceToCellMap(face, axis, 0);
			UT_Vector3I forwardCell = faceToCellMap(face, axis, 1);

			assert(backwardCell[axis] >= 0 && forwardCell[axis] < materialCellLabels.getVoxelRes()[axis]);

			exint backwardMaterial = getFieldValue(materialCellLabels, backwardCell);
			exint forwardMaterial = getFieldValue(materialCellLabels, forwardCell);

			if (backwardMaterial == MaterialLabels::INTERIOR_LIQUID_CELL || 
			    forwardMaterial == MaterialLabels::INTERIOR_LIQUID_CELL)
			{
			    exint interiorRegion;

			    if (backwardMaterial == MaterialLabels::INTERIOR_LIQUID_CELL)
			    {
				interiorRegion = getFieldValue(interiorRegionIndices, backwardCell);

				assert(interiorRegion >= 0);
				assert(getFieldValue(activeCellIndices, backwardCell) == HDK::Utilities::UNLABELLED_CELL);

				if (forwardMaterial == MaterialLabels::INTERIOR_LIQUID_CELL)
				{
				    assert(interiorRegion == getFieldValue(interiorRegionIndices, forwardCell));
				    assert(getFieldValue(activeCellIndices, forwardCell) == HDK::Utilities::UNLABELLED_CELL);
				}
				else
				{
				    assert(forwardMaterial == MaterialLabels::ACTIVE_LIQUID_CELL);
				    assert(getFieldValue(interiorRegionIndices, forwardCell) == HDK::Utilities::UNLABELLED_CELL);
				    assert(getFieldValue(activeCellIndices, forwardCell) >= 0);
				}
			    }
			    else
			    {
				interiorRegion = getFieldValue(interiorRegionIndices, forwardCell);

				assert(interiorRegion >= 0);
				assert(getFieldValue(activeCellIndices, forwardCell) == HDK::Utilities::UNLABELLED_CELL);

				assert(backwardMaterial == MaterialLabels::ACTIVE_LIQUID_CELL);
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

#if !defined(NDEBUG)
			    if (backwardMaterial == MaterialLabels::ACTIVE_LIQUID_CELL)
			    {
				assert(getFieldValue(activeCellIndices, backwardCell) >= 0);
				assert(getFieldValue(interiorRegionIndices, backwardCell) == HDK::Utilities::UNLABELLED_CELL);
			    }
			    else
			    {
				assert(backwardMaterial == MaterialLabels::AIR_CELL);
				assert(forwardMaterial == MaterialLabels::ACTIVE_LIQUID_CELL);				
			    }

			    if (forwardMaterial == MaterialLabels::ACTIVE_LIQUID_CELL)
			    {
				assert(getFieldValue(activeCellIndices, forwardCell) >= 0);
				assert(getFieldValue(interiorRegionIndices, forwardCell) == HDK::Utilities::UNLABELLED_CELL);
			    }
			    else
			    {
				assert(backwardMaterial == MaterialLabels::ACTIVE_LIQUID_CELL);
				assert(forwardMaterial == MaterialLabels::AIR_CELL);				
			    }

#endif

			    SolveReal gradient = getFieldValue(pressure, forwardCell) - getFieldValue(pressure, backwardCell);

			    if (backwardMaterial == MaterialLabels::AIR_CELL || forwardMaterial == MaterialLabels::AIR_CELL)
			    {
				SolveReal phi0 = getFieldValue(liquidSurface, backwardCell);
				SolveReal phi1 = getFieldValue(liquidSurface, forwardCell);

				SolveReal theta = HDK::Utilities::computeGhostFluidWeight(phi0, phi1);
				theta = SYSclamp(theta, .01, 1.);
				gradient /= theta;
			    }

			    setFieldValue(velocity, face, getFieldValue(velocity, face) - gradient);
			}
		    }
		}
	    }
	}
    });
}

void
HDK_AffineFreeSurfacePressureSolver::computeResultingDivergence(UT_Array<SolveReal> &parallelAccumulatedDivergence,
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
		break;

	    if (!vit.isTileConstant() ||
		vit.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
		vit.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
	    {
		vitt.setTile(vit);

		for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
		{
		    if (vitt.getValue() == MaterialLabels::INTERIOR_LIQUID_CELL ||
			vitt.getValue() == MaterialLabels::ACTIVE_LIQUID_CELL)
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

			localAccumulatedDivergence += fabs(divergence);
			localMaxDivergence = std::max(localMaxDivergence, fabs(divergence));
			++localCellCount;
		    }
		}
	    }
	}

	return 0;
    });
}