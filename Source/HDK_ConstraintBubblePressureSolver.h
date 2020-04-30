#ifndef HDK_CONSTRAINT_BUBBLE_PRESSURE_SOLVER_H
#define HDK_CONSTRAINT_BUBBLE_PRESSURE_SOLVER_H

#include "Eigen/Sparse"

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>

#include <GA/GA_Handle.h>

#include <SIM/SIM_RawIndexField.h>

#include <UT/UT_ThreadedAlgorithm.h>

class GU_Detail;

class SIM_ScalarField;
class SIM_VectorField;

#include "../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridPoissonSolver.h"

#include "HDK_ReducedFluidUtilities.h"

class GAS_API HDK_ConstraintBubblePressureSolver : public GAS_SubSolver
{
    using SolveReal = double;
    using Vector = Eigen::VectorXd;

    static constexpr exint VISITED_CELL = -1;
    static constexpr exint UNVISITED_CELL = -2;

    static constexpr int UNINITIALIZED_PARTICLE = -2;
    static constexpr int UNLABELLED_PARTICLE = -1;

    enum MaterialLabels { SOLID_CELL, LIQUID_CELL, BUBBLE_CELL };

public:

    GET_DATA_FUNC_F(SIM_NAME_TOLERANCE, Tolerance);
    GET_DATA_FUNC_I("maxIterations", MaxIterations);

    GET_DATA_FUNC_B("useOldPressure", UseOldPressure);

    GET_DATA_FUNC_B("useOpenBoundaries", UseOpenBoundaries);

    GET_DATA_FUNC_B("trackBubbleIDs", TrackBubbleIDs);
    GET_DATA_FUNC_B("trackBubbleVolumes", TrackBubbleVolumes);
    GET_DATA_FUNC_B("correctVolumeDrift", CorrectVolumeDrift);

    GET_DATA_FUNC_B("useMGPreconditioner", UseMGPreconditioner);

    GET_DATA_FUNC_S(GAS_NAME_GEOMETRY, GeometryName);

protected:

    explicit HDK_ConstraintBubblePressureSolver(const SIM_DataFactory *factory);
    virtual ~HDK_ConstraintBubblePressureSolver();

    // Used to determine if the field is complicated enough to justify
    // the overhead of multithreading.
    bool shouldMultiThread(const SIM_RawField *field) const
    {
	return field->field()->numTiles() > 1;
    }

    // The overloaded callback that GAS_SubSolver will invoke to
    // perform our actual computation.  We are giving a single object
    // at a time to work on.
    virtual bool solveGasSubclass(SIM_Engine &engine,
				    SIM_Object *obj,
				    SIM_Time time,
				    SIM_Time timestep);
private:

    // We define this to be a DOP_Auto node which means we do not
    // need to implement a DOP_Node derivative for this data.  Instead,
    // this description is used to define the interface.
    static const SIM_DopDescription *getDopDescription();
    
    /// These macros are necessary to bind our node to the factory and
    /// ensure useful constants like BaseClass are defined.
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(HDK_ConstraintBubblePressureSolver,
                        GAS_SubSolver,
                        "HDK Constraint Bubble Pressure Solver",
                        getDopDescription());

    ////////////////////////////////////////////
    //
    // Solver methods
    //
    ////////////////////////////////////////////

    void
    setMaterialCellLabels(SIM_RawIndexField &materialCellLabels) const;

    void
    buildValidFaces(SIM_VectorField &validFaces,
		    const SIM_RawIndexField &liquidCellIndices,
		    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    exint
    buildLiquidCellIndices(SIM_RawIndexField &liquidCellIndices,
			    const SIM_RawIndexField &materialCellLabels) const;

    void
    integrateSurfacePressure(SIM_RawField &velocity,
				const SIM_RawField &validFaces,
				const SIM_RawField &surfacePressure,
				const SIM_RawField &liquidSurface,
				const SIM_RawIndexField &liquidCellIndices,
				const SolveReal liquidDensity,
				const SolveReal dt,
				const int axis) const;

    void
    computeBubbleRegionSize(UT_Array<UT_Array<exint>> &parallelBubbleRegionSize,
			    const SIM_RawIndexField &bubbleRegionIndices) const;

    exint
    buildBubbleRegionIndices(SIM_RawIndexField &bubbleRegionIndices,
				UT_Array<exint> &bubbleRegionSize,
				const SIM_RawIndexField &materialCellLabels,
				std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    exint
    buildCombinedRegionIndices(SIM_RawIndexField &combinedRegionIndices,
				const SIM_RawIndexField &materialCellLabels,
				std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    buildRegionBooleans(UT_Array<UT_Array<bool>> &isBubbleInRegion,
			UT_Array<bool> &isLiquidInRegion,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &liquidCellIndices,
			const SIM_RawIndexField &bubbleRegionIndices,
			const SIM_RawIndexField &combinedRegionIndices,
			const exint combinedRegionCount,
			const exint bubbleRegionCount) const;

    void
    trackBubbleIDs(UT_Array<bool> &isBubbleTracked,
		    UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
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
    distributeOldBubbleVolumes(UT_Array<SolveReal> &newBubbleRegionVolumes,
				UT_Array<bool> &isBubbleUninitialized,
				const GA_RWHandleF &bubbleVolumeHandle,
				const UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
				const UT_Array<bool> &isBubbleTracked,
				const UT_Array<exint> &bubbleRegionSizes) const;

    void
    removeBubblesAtOpenBoundaries(UT_Array<bool> &isBubbleRemoved,
				    const SIM_RawIndexField &bubbleRegionIndices,
				    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    buildLiquidRHS(Vector &rhs,
		    const SIM_RawIndexField &materialCellLabels,
		    const SIM_RawIndexField &liquidCellIndices,
		    const SIM_RawIndexField &bubbleRegionIndices,
		    const SIM_RawIndexField &combinedRegionIndices,
		    const SIM_VectorField &velocity,
		    const SIM_VectorField *solidVelocity,
		    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    buildBubbleRHS(UT_Array<UT_Array<SolveReal>> &parallelBubbleRegionDivergence,
		    const SIM_RawIndexField &materialCellLabels,
		    const SIM_RawIndexField &liquidCellIndices,
		    const SIM_RawIndexField &bubbleRegionIndices,
		    const SIM_RawIndexField &combinedRegionIndices,
		    const std::array<const SIM_RawField *, 3> &cutCellWeights,
		    const UT_Array<bool> &isBubbleRemoved,
		    const SIM_VectorField &velocity,
		    const SIM_VectorField *solidVelocity,
		    const exint bubbleRegionCount) const;

    void
    buildLiquidRegionDivergence(UT_Array<UT_Array<SolveReal>> &parallelLiquidRegionDivergence,
				UT_Array<UT_Array<exint>> &parallelLiquidRegionCellCount,
				const SIM_RawIndexField &combinedRegionIndices,
				const SIM_RawIndexField &liquidCellIndices,
				const UT_Array<bool> &doesRegionHaveRemovedBubble,
				const Vector &rhs) const;

    void
    removeDivergenceFromLiquidCells(Vector &rhs,
				    const UT_Array<SolveReal> &regionDivergence,
				    const SIM_RawIndexField &combinedRegionIndices,
				    const SIM_RawIndexField &liquidCellIndices,
				    const UT_Array<bool> &doesRegionHaveRemovedBubble) const;

    void
    buildLiquidRows(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelSparseElements,
		    const SIM_RawIndexField &liquidCellIndices,
		    const SIM_RawIndexField &bubbleRegionIndices,
		    const SIM_RawField &liquidSurface,
		    const std::array<const SIM_RawField *, 3> &cutCellWeights,
		    const UT_Array<bool> &isBubbleRemoved,
		    const SolveReal liquidDensity,
		    const exint liquidCellCount) const;

    void
    buildBubbleRows(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelSparseElements,
		    UT_Array<UT_Array<SolveReal>> &parallelBubbleRegionDiagonals,
		    const SIM_RawIndexField &liquidCellIndices,
		    const SIM_RawIndexField &bubbleRegionIndices,
		    const SIM_RawField &liquidSurface,
		    const std::array<const SIM_RawField *, 3> &cutCellWeights,
		    const UT_Array<bool> &isBubbleRemoved,
		    const SolveReal liquidDensity,
		    const exint liquidCellCount) const;

    void
    computeBubblePressure(UT_Array<UT_Array<SolveReal>> &parallelAccumulatedPressure,
			    const SIM_RawField &pressure,
			    const SIM_RawIndexField &bubbleRegionIndices,
			    const UT_Array<bool> &isBubbleRemoved,
			    const exint bubbleRegionCount) const;
    //
    // Multigrid operations
    //

    void
    buildMGDomainLabels(UT_VoxelArray<int> &mgDomainCellLabels,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &liquidCellIndices) const;

    void
    buildMGBoundaryWeights(UT_VoxelArray<SolveReal> &boundaryWeights,
			    const SIM_RawField &validFaces,
			    const SIM_RawIndexField &materialCellLabels,
			    const SIM_RawIndexField &bubbleRegionIndices,
			    const UT_VoxelArray<int> &domainCellLabels,
			    const SIM_RawField &liquidSurface,
			    const SIM_RawField &cutCellWeights,
			    const SolveReal liquidDensity,
			    const int axis) const;

    bool
    unitTestMGLabels(const UT_VoxelArray<int> &mgDomainCellLabels,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &liquidCellIndices,
			const SIM_RawIndexField &bubbleRegionIndices,
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
    buildSmootherCells(UT_Array<UT_Array<UT_Vector3I>> &parallelBubbleSmootherCells,
			UT_Array<UT_Array<UT_Vector3I>> &parallelLiquidCopyCells,
			UT_Array<UT_Array<UT_Vector3I>> &parallelBubbleBoundaryCells,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &visitedLiquidCells,
			const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    copySourceToSmootherGrid(UT_VoxelArray<SolveReal> &smootherSourceGrid,
				const Vector &sourceVector,
				const SIM_RawIndexField &liquidCellIndices,
				const UT_Array<UT_Vector3I> &bubbleSmootherCells) const;

    void
    buildBubbleLaplacian(UT_Array<UT_Array<SolveReal>> &parallelBubbleLaplacian,
			    const SIM_RawIndexField &materialCellLabels,
			    const SIM_RawIndexField &bubbleRegionIndices,
			    const SIM_RawField &liquidSurface,
			    const std::array<const SIM_RawField *, 3> &cutCellWeights,
			    const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
			    const UT_Array<UT_Vector3I> &bubbleBoundaryCells,
			    const UT_Array<bool> &isBubbleRemoved,
			    const SolveReal liquidDensity) const;

    void
    applyBubbleSmoother(UT_Array<SolveReal> &bubbleRegionDestinations,
			const UT_Array<SolveReal> &bubbleRegionSources,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &bubbleRegionIndices,
			const SIM_RawField &liquidSurface,
			const std::array<const SIM_RawField *, 3> &cutCellWeights,
			const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
			const UT_Array<SolveReal> &bubbleRegionDiagonals,
			const UT_Array<UT_Vector3I> &bubbleBoundaryCells,
			const UT_Array<bool> &isBubbleRemoved,
			const SolveReal liquidDensity) const;

    void
    applyLiquidSmoother(UT_Array<SolveReal> &tempSmootherDestinationValues,
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
			const SolveReal liquidDensity) const;

    void
    copyTempDestinationToGrid(UT_VoxelArray<SolveReal> &smootherDestinationGrid,
				const UT_Array<SolveReal> &tempSmootherDestinationValues,
				const UT_Array<UT_Vector3I> &bubbleSmootherCells) const;

    void
    copySourceToMGGrid(UT_VoxelArray<SolveReal> &mgSourceGrid,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &liquidCellIndices,
			const UT_VoxelArray<int> &mgDomainCellLabels,
			const Vector &sourceVector,
			const UT_Vector3I &mgExpandedOffset,
			const SolveReal liquidDensity) const;

    void
    applyDirichletToMGGrid(UT_VoxelArray<SolveReal> &mgSourceGrid,
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
			    const SolveReal liquidDensity) const;

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
				const SIM_RawIndexField &liquidCellIndices,
				const UT_VoxelArray<int> &mgDomainCellLabels,
				const UT_VoxelArray<SolveReal> &mgDestinationGrid,
				const UT_Vector3I &mgExpandedOffset) const;

    void
    copySmootherToDestinationVector(Vector &destinationVector,
				    const SIM_RawIndexField &liquidCellIndices,
				    const UT_VoxelArray<SolveReal> &smootherDestinationGrid,
				    const UT_Array<UT_Vector3I> &bubbleSmootherCells) const;

    void
    applySolutionToPressure(SIM_RawField &pressure,
			    const SIM_RawIndexField &indices,
			    const Vector &solution,
			    const exint offset) const;

    void
    applyPressureGradient(SIM_RawField &velocity,
			    const SIM_RawField &validFaces,
			    const SIM_RawField &pressure,
			    const SIM_RawField &liquidSurface,
			    const SIM_RawIndexField &materialCellLabels,
			    const SIM_RawIndexField &bubbleRegionIndices,
			    const SolveReal &liquidDensity,
			    const int axis) const;

    void
    computeResultingDivergence(UT_Array<SolveReal> &parallelAccumulatedLiquidDivergence,
				UT_Array<SolveReal> &parallelMaxLiquidDivergence,
				UT_Array<SolveReal> &parallelLiquidCellCount,
				UT_Array<UT_Array<SolveReal>> &parallelBubbleRegionDivergence,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_RawIndexField &bubbleRegionIndices,
				const SIM_VectorField &velocity,
				const std::array<const SIM_RawField *, 3> &cutCellWeights,
				const SIM_VectorField *solidVelocity,
				const UT_Array<bool> &isBubbleTracked,
				const UT_Array<bool> &isBubbleRemoved) const;
};
#endif