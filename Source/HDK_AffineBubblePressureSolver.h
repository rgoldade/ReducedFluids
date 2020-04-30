#ifndef HDK_AFFINE_BUBBLE_PRESSURE_SOLVER_H
#define HDK_AFFINE_BUBBLE_PRESSURE_SOLVER_H

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>

#include <SIM/SIM_RawIndexField.h>

#include <UT/UT_ThreadedAlgorithm.h>

class SIM_VectorField;
class SIM_ScalarField;

#include "../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridPoissonSolver.h"

#include "HDK_ReducedFluidUtilities.h"

class GAS_API HDK_AffineBubblePressureSolver : public GAS_SubSolver
{
    using SolveReal = double;
    using Vector = Eigen::VectorXd;
    using UT_Vector3SR = UT_Vector3T<SolveReal>;

    static constexpr exint VISITED_CELL = -1;
    static constexpr exint UNVISITED_CELL = -2;

    static constexpr int UNINITIALIZED_PARTICLE = -2;
    static constexpr int UNLABELLED_PARTICLE = -1;

    enum MaterialLabels { SOLID_CELL, EXTERIOR_LIQUID_CELL, INTERIOR_LIQUID_CELL, EXTERIOR_BUBBLE_CELL, INTERIOR_BUBBLE_CELL };
    
    using GradientMatrix = Eigen::Matrix<SolveReal, 3, 3>;
    using MassMatrix = Eigen::Matrix<SolveReal, 11, 11>;
    using ColumnVector = Eigen::Matrix<SolveReal, 11, 1>;

public:

    GET_DATA_FUNC_F(SIM_NAME_TOLERANCE, Tolerance);
    GET_DATA_FUNC_I("maxIterations", MaxIterations);

    GET_DATA_FUNC_B("useOldPressure", UseOldPressure);

    GET_DATA_FUNC_I("exteriorBandwidth", ExteriorBandwidth);

    GET_DATA_FUNC_F("bubbleDensity", BubbleDensity);

    GET_DATA_FUNC_B("useMGPreconditioner", UseMGPreconditioner);

    GET_DATA_FUNC_B("trackBubbleIDs", TrackBubbleIDs);

    GET_DATA_FUNC_B("useTiledInterior", UseTiledInterior);
    GET_DATA_FUNC_I("tileSize", TileSize);

    GET_DATA_FUNC_B("applyAffineToLiquids", ApplyAffineToLiquids);

    GET_DATA_FUNC_B("printCells", PrintCells);
    GET_DATA_FUNC_B("onlyPrintCells", OnlyPrintCells);

    GET_DATA_FUNC_S(GAS_NAME_GEOMETRY, GeometryName);

protected:

    explicit HDK_AffineBubblePressureSolver(const SIM_DataFactory *factory);
    virtual ~HDK_AffineBubblePressureSolver();

    // The overloaded callback that GAS_SubSolver will invoke to
    // perform our actual computation.  We are giving a single object
    // at a time to work on.
    virtual bool solveGasSubclass(SIM_Engine &engine,
				    SIM_Object *obj,
				    SIM_Time time,
				    SIM_Time timestep);

private:

    /// These macros are necessary to bind our node to the factory and
    /// ensure useful constants like BaseClass are defined.
    static const SIM_DopDescription	*getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(HDK_AffineBubblePressureSolver,
                        GAS_SubSolver,
                        "HDK Affine Bubble Pressure Solver",
                        getDopDescription());

    ////////////////////////////////////////////
    //
    // Solver methods
    //
    ////////////////////////////////////////////

    void
    setMaterialCellLabels(SIM_RawIndexField &materialCellLabels) const;

    //
    // Helper methods to build affine regions for liquid and bubble volumes
    //

    void
    buildInitialFluidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExteriorCellList,
				    const SIM_RawIndexField &materialCellLabels,
				    const std::array<const SIM_RawField *, 3> &cutCellWeights,
				    const exint interiorMaterialLabel,
				    const exint opposingInteriorMaterialLabel) const;

    void
    buildNextFluidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExteriorCellList,
				const UT_Array<UT_Vector3I> &oldExteriorCellList,
				const SIM_RawIndexField &materialCellLabels,
				const std::array<const SIM_RawField *, 3> &cutCellWeights,
				const exint exteriorMaterialLabel,
				const exint interiorMaterialLabel) const;

    void
    buildInitialSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExteriorCellList,
				    const SIM_RawIndexField &materialCellLabels,
				    const std::array<const SIM_RawField *, 3> &cutCellWeights,
				    const exint interiorMaterialLabel,
				    const exint exteriorMaterialLabel) const;

    void
    setSolidLayerCells(SIM_RawIndexField &materialCellLabels,
			SIM_RawIndexField &visitedCells,
			const UT_Array<UT_Vector3I> &newExteriorCellList,
			const exint interiorMaterialLabel,
			const exint exteriorMaterialLabel) const;

    void
    buildNextSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExteriorCellList,
				const UT_Array<UT_Vector3I> &oldExteriorCellList,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_RawIndexField &visitedCells,
				const std::array<const SIM_RawField *, 3> &cutCellWeights,
				const exint interiorMaterialLabel,
				const exint exteriorMaterialLabel) const;

    void
    buildExteriorCells(SIM_RawIndexField &materialCellLabels,
			const std::array<const SIM_RawField *, 3> &cutCellWeights,
			const exint exteriorMaterialLabel,
			const exint interiorMaterialLabel,
			const exint opposingInteriorMaterialLabel) const;

    void
    buildTiledActiveCells(SIM_RawIndexField &materialCellLabels,
			    const int tileSize,
			    const exint exteriorMaterialLabel,
			    const exint interiorMaterialLabel) const;

    void
    setAllInteriorCellsExterior(SIM_RawIndexField &materialCellLabels,
				const exint exteriorMaterialLabel,
				const exint interiorMaterialLabel) const;

    //
    // Helper methods to track bubble IDs
    //

    void
    trackBubbleIDs(UT_Array<bool> &isBubbleTracked,
		    const GA_RWHandleI &bubbleIDHandle,
		    const GU_Detail &particles,
		    const SIM_RawIndexField &bubbleRegionIndices,
		    const SIM_RawField &liquidSurface) const;

    void
    updateBubbleIDs(GA_RWHandleI &bubbleIDHandle,
		    const UT_Array<bool> &isBubbleTracked,
		    const GU_Detail &particles,
		    const SIM_RawIndexField &bubbleRegionIndices,
		    const SIM_RawField &liquidSurface) const;

    void
    setUntrackedBubbleCells(SIM_RawIndexField &materialCellLabels,
			    const SIM_RawIndexField &bubbleRegionIndices,
			    const UT_Array<bool> &isBubbleTracked) const;

    void
    buildBubbleRegions(SIM_RawIndexField &fullBubbleRegionIndices,
			UT_Array<bool> &isBubbleTracked,
			SIM_RawIndexField &materialCellLabels,
			SIM_Object *obj,
			const SIM_RawField &liquidSurface,
			std::array<const SIM_RawField *, 3> &cutCellWeights);

    //
    // Helper methods to remove single cell regions
    //

    void
    buildInteriorBoundingBoxes(UT_Array<UT_Array<UT_Vector3I>> &parallelInteriorRegionBBMin,
				UT_Array<UT_Array<UT_Vector3I>> &parallelInteriorRegionBBMax,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_RawIndexField &interiorRegionIndices) const;

    void
    remapInteriorRegions(SIM_RawIndexField &interiorRegionIndices,
			    SIM_RawIndexField &materialCellLabels,
			    const UT_Array<bool> &doRemoveRegion,
			    const UT_Array<exint> &remapRegion) const;
   
    exint
    removeDegenerateCellRegions(SIM_RawIndexField &interiorRegionIndices,
				SIM_RawIndexField &materialCellLabels,
				const exint interiorRegionCount) const;

    exint
    buildActiveCellIndices(SIM_RawIndexField &activeCellLabels,
			    const SIM_RawIndexField &materialCellLabels) const;

    void
    printCellGeometry(SIM_Object *obj,
			const SIM_RawIndexField &activeCellIndices,
			const SIM_RawIndexField &interiorRegionIndices,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_VectorField &velocity,
			const fpreal dx);

    void
    buildInteriorRegionCOM(UT_Array<UT_Array<UT_Vector3SR>> &parallelInteriorRegionCOM,
			    UT_Array<UT_Array<SolveReal>> &parallelInteriorRegionCellCount,
			    const SIM_RawIndexField &interiorRegionIndices,
			    const SIM_RawIndexField &materialCellLabels) const;

    void
    buildBestFitSystems(UT_Array<UT_Array<MassMatrix>> &parallelBestFitMatrix,
			UT_Array<UT_Array<ColumnVector>> &parallelBestFitRHS,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &interiorRegionIndices,
			const SIM_VectorField &velocity,
			const UT_Array<UT_Vector3SR> &interiorRegionCOM,
			const SolveReal dx) const;

    void
    buildRHS(Vector &rhsVector,
		const SIM_RawIndexField &materialCellLabels,
		const SIM_RawIndexField &interiorRegionIndices,	
		const SIM_RawIndexField &activeCellIndices,
		const SIM_VectorField &velocity,
		const std::array<const SIM_RawField *, 3> &cutCellWeights,
		const SIM_VectorField *solidVelocity,
		const UT_Array<UT_Vector3SR> &interiorLinearVelocities,
		const UT_Array<GradientMatrix> &interiorVelocityGradients,
		const UT_Array<UT_Vector3SR> &interiorRegionCOM,
		const SolveReal dx) const;

    void
    applyGoalDivergence(Vector &rhsVector,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &activeCellIndices,
			const SIM_RawField &goalDivergence) const;

    void
    buildActiveRegionDivergence(UT_Array<UT_Array<SolveReal>> &parallelActiveRegionDivergence,
				UT_Array<UT_Array<SolveReal>> &parallelActiveRegionCellCount,
				const SIM_RawIndexField &activeCellIndices,
				const SIM_RawIndexField &activeRegionIndices,
				const Vector &rhsVector) const;

    void
    removeAverageDivergence(Vector &rhsVector,
			    const SIM_RawIndexField &activeCellIndices,
			    const SIM_RawIndexField &activeRegionIndices,
			    const UT_Array<SolveReal> &averageActiveRegionDivergence) const;

    void
    setUntrackedDivergence(UT_Array<UT_Array<SolveReal>> &parallelUntrackedDivergence,
			    UT_Array<UT_Array<exint>> &parallelTrackedCellCount,
			    UT_Array<UT_Array<exint>> &parallelUntrackedCellCount,
			    Vector &rhsVector,
			    const SIM_RawIndexField &materialCellLabels,
			    const SIM_RawIndexField &activeCellIndices,
			    const SIM_RawIndexField &activeRegionIndices,
			    const SIM_RawIndexField &fullBubbleRegionIndices,
			    const UT_Array<bool> &isBubbleTracked,
			    const SolveReal dx, const SolveReal dt) const;

    void
    removeUntrackedDivergence(Vector &rhsVector,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_RawIndexField &activeCellIndices,
				const SIM_RawIndexField &activeRegionIndices,
				const SIM_RawIndexField &fullBubbleRegionIndices,
				const UT_Array<bool> &isBubbleTracked,
				const UT_Array<SolveReal> &untrackedDivergence) const;

    void
    buildLinearSystem(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelPoissonElements,
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
			const SolveReal dx) const;

    void
    buildAffineMassMatrixSystems(UT_Array<UT_Array<MassMatrix>> &parallelAffineMassMatrix,
				    const SIM_RawIndexField &materialCellLabels,
				    const SIM_RawIndexField &interiorRegionIndices,
				    const UT_Array<UT_Vector3SR> &interiorRegionCOM,
				    const SolveReal dx) const;

    void
    buildInteriorBoundaryCells(UT_Array<UT_Array<UT_Vector3I>> &parallelInteriorBoundaryCells,
				const SIM_RawIndexField &materialCellLabels) const;

    void
    distributeAffineVectors(Vector &destinationVector,
			    const Vector &invertedAffineVectors,
			    const UT_Array<UT_Vector3I> &interiorBoundaryCells,
			    const SIM_RawIndexField &materialCellLabels,
			    const SIM_RawIndexField &interiorRegionIndices,
			    const SIM_RawIndexField &activeCellIndices,
			    const UT_Array<UT_Vector3SR> &interiorRegionCOM,
			    const SolveReal dx) const;

    //
    // Multigrid operations
    //

    void
    buildMGDomainLabels(UT_VoxelArray<int> &mgDomainCellLabels,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &activeCellIndices) const;

    void
    buildMGBoundaryWeights(UT_VoxelArray<SolveReal> &boundaryWeights,
			    const SIM_RawField &validFaces,
			    const SIM_RawIndexField &materialCellLabels,
			    const SIM_RawIndexField &interiorRegionIndices,
			    const UT_VoxelArray<int> &domainCellLabels,
			    const SIM_RawField &liquidSurface,
			    const SIM_RawField &cutCellWeights,
			    const SolveReal liquidDensity,
			    const SolveReal bubbleDensity,
			    const int axis) const;

    bool
    unitTestMGLabels(const UT_VoxelArray<int> &mgDomainCellLabels,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &activeCellIndices,
			const SIM_RawIndexField &interiorRegionIndices,
			const UT_Vector3I &mgExpandedOffset) const;

    void
    buildInitialLiquidSmootherLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidBandCells,
				    const SIM_RawIndexField &materialCellLabels,
				    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;
    void
    setVisitedLiquidCells(SIM_RawIndexField &visitedLiquidCells,
			    const UT_Array<UT_Vector3I> &newLiquidCells) const;

    void
    buildNextLiquidSmootherLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidLayerCells,
				    const UT_Array<UT_Vector3I> &oldLiquidLayerCells,
				    const SIM_RawIndexField &materialCellLabels,
				    const SIM_RawIndexField &visitedLiquidCells,
				    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    buildBubbleSmootherCells(UT_Array<UT_Array<UT_Vector3I>> &parallelBubbleSmootherCells,
				UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidCopyCells,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_RawIndexField &visitedLiquidCells,
				const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    findOccupiedTiles(UT_Array<bool> &isSmootherTileOccupied,
			const UT_VoxelArray<SolveReal> &voxelArrayGrid,
			const UT_Array<UT_Vector3I> &activeCells) const;

    void
    copySourceToSmootherGrid(UT_VoxelArray<SolveReal> &smootherSourceGrid,
				const Vector &sourceVector,
				const SIM_RawIndexField &activeCellIndices,
				const UT_Array<UT_Vector3I> &bubbleSmootherCells) const;

    void
    buildAffineVectors(UT_Array<ColumnVector> &affineVectors,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &interiorRegionIndices,
			const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
			const UT_Array<UT_Vector3I> &interiorBoundaryCells,
			const UT_Array<UT_Vector3SR> &interiorRegionCOM,
			const exint interioRegionCount,
			const SolveReal dx) const;

    void
    applyBubbleSmoother(UT_Array<SolveReal> &tempDestinationValues,
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
			const SolveReal dx) const;

    void
    copyDestinationVectorToGrid(UT_VoxelArray<SolveReal> &smootherDestinationGrid,
				const UT_Array<SolveReal> &tempDestinationValues,
				const UT_Array<UT_Vector3I> &bubbleSmootherCells) const;

    void
    copySourceToMGGrid(UT_VoxelArray<SolveReal> &mgSourceGrid,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &activeCellIndices,
			const UT_VoxelArray<int> &mgDomainCellLabels,
			const Vector &sourceVector,
			const UT_Vector3I &mgExpandedOffset,
			const SolveReal liquidDensity) const;

    void
    applyDirichletToMGGrid(UT_VoxelArray<SolveReal> &mgSourceGrid,
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
			    const SolveReal bubbleDensity) const;

    void
    copyMGToSmoother(UT_VoxelArray<SolveReal> &smootherDestinationGrid,
			const SIM_RawIndexField &materialCellLabels,
			const UT_VoxelArray<int> &mgDomainCellLabels,
			const UT_VoxelArray<SolveReal> &mgDestinationGrid,
			const UT_Array<UT_Vector3I> &liquidCopyCells,
			const UT_Vector3I &mgExpandedOffset) const;

    void
    copyMGToDestinationVector(Vector &destinationVector,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_RawIndexField &activeCellIndices,
				const UT_VoxelArray<int> &mgDomainCellLabels,
				const UT_VoxelArray<SolveReal> &mgDestinationGrid,
				const UT_Vector3I &mgExpandedOffset) const;

    void
    copySmootherToDestinationVector(Vector &destinationVector,
				    const SIM_RawIndexField &activeCellIndices,
				    const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
				    const UT_Array<UT_Vector3I> &bubbleSmootherCells) const;

    void
    buildValidFaces(SIM_VectorField &validFaces,
		    const SIM_RawIndexField &materialCellLabels,
		    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    copyPressureVectorToGrid(SIM_RawField &pressure,
				const SIM_RawIndexField &activeCellIndices,
				const Vector &solutionVector) const;

    void
    applySolutionToVelocity(SIM_RawField &velocity,
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
			    const SolveReal dx) const;

    void
    computeResultingDivergence(UT_Array<SolveReal> &parallelAccumulatedDivergence,
				UT_Array<SolveReal> &parallelMaxDivergence,
				UT_Array<SolveReal> &parallelCellCount,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_VectorField &velocity,
				const std::array<const SIM_RawField *, 3> &cutCellWeights,
				const SIM_VectorField *solidVelocity) const;
};

#endif