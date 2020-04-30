#ifndef HDK_TWO_PHASE_PRESSURE_SOLVER_H
#define HDK_TWO_PHASE_PRESSURE_SOLVER_H

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>

#include <SIM/SIM_RawIndexField.h>

#include <UT/UT_ThreadedAlgorithm.h>

class SIM_VectorField;
class SIM_ScalarField;

#include "../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridPoissonSolver.h"

#include "HDK_ReducedFluidUtilities.h"

class GAS_API HDK_TwoPhasePressureSolver : public GAS_SubSolver
{
    using SolveReal = double;
    using Vector = Eigen::VectorXd;
    using UT_Vector3SR = UT_Vector3T<SolveReal>;

    static constexpr exint VISITED_CELL = -1;
    static constexpr exint UNVISITED_CELL = -2;

    static constexpr int UNINITIALIZED_PARTICLE = -2;
    static constexpr int UNLABELLED_PARTICLE = -1;

    enum MaterialLabels { SOLID_CELL, LIQUID_CELL, BUBBLE_CELL};
    
public:

    GET_DATA_FUNC_F(SIM_NAME_TOLERANCE, Tolerance);
    GET_DATA_FUNC_I("maxIterations", MaxIterations);

    GET_DATA_FUNC_B("useOldPressure", UseOldPressure);

    GET_DATA_FUNC_F("bubbleDensity", BubbleDensity);

    GET_DATA_FUNC_B("trackBubbleIDs", TrackBubbleIDs);
    GET_DATA_FUNC_B("trackBubbleVolumes", TrackBubbleVolumes);
    GET_DATA_FUNC_B("correctVolumeDrift", CorrectVolumeDrift);

    GET_DATA_FUNC_B("printCells", PrintCells);
    GET_DATA_FUNC_B("onlyPrintCells", OnlyPrintCells);

    GET_DATA_FUNC_S(GAS_NAME_GEOMETRY, GeometryName);

protected:

    explicit HDK_TwoPhasePressureSolver(const SIM_DataFactory *factory);
    virtual ~HDK_TwoPhasePressureSolver();

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
    DECLARE_DATAFACTORY(HDK_TwoPhasePressureSolver,
                        GAS_SubSolver,
                        "HDK Two Phase Pressure Solver",
                        getDopDescription());

    ////////////////////////////////////////////
    //
    // Solver methods
    //
    ////////////////////////////////////////////

    void
    setMaterialCellLabels(SIM_RawIndexField &materialCellLabels) const;

    exint
    buildActiveCellIndices(SIM_RawIndexField &activeCellLabels,
			    const SIM_RawIndexField &materialCellLabels) const;

    void
    computeBubbleRegionSize(UT_Array<UT_Array<exint>> &parallelBubbleRegionSizes,
			    const SIM_RawIndexField &bubbleRegionIndices) const;
    exint
    buildBubbleRegionIndices(SIM_RawIndexField &bubbleRegionIndices,
				UT_Array<exint> &bubbleRegionSizes,
				const SIM_RawIndexField &materialCellLabels,
				std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    exint
    buildActiveRegionIndices(SIM_RawIndexField &activeRegionIndices,
				const SIM_RawIndexField &materialCellLabels,
				std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    //
    // Helper methods to track bubble IDs
    //

    void
    trackBubbleIDs(UT_Array<bool> &areBubblesTracked,
		    UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
		    const GA_RWHandleI &bubbleIDHandle,
		    const GU_Detail &particles,
		    const SIM_RawIndexField &bubbleRegionIndices,
		    const SIM_RawField &liquidSurface) const;

    void
    updateBubbleIDs(GA_RWHandleI &bubbleIDHandle,
		    const UT_Array<bool> &areBubblesTracked,
		    const GU_Detail &particles,
		    const SIM_RawIndexField &bubbleRegionIndices,
		    const SIM_RawField &liquidSurface) const;

    void
    buildBubbleIDs(UT_Array<bool> &areBubblesTracked,
		    UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
		    SIM_RawIndexField &materialCellLabels,
		    SIM_Object *obj,
		    const SIM_RawField &liquidSurface,
		    const SIM_RawIndexField &bubbleRegionIndices,
		    std::array<const SIM_RawField *, 3> &cutCellWeights,
		    const exint bubbleRegionCount,
		    const bool doTrackBubbleVolumes);

    void
    distributeOldBubbleVolumes(UT_Array<SolveReal> &newBubbleRestVolumes,
				UT_Array<bool> &areBubblesUninitialized,
				const GA_RWHandleF &bubbleVolumeHandle,
				const UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
				const UT_Array<bool> &areBubblesTracked,
				const UT_Array<exint> &newBubbleRegionSizes) const;

    void
    buildBubbleRestVolumes(UT_Array<SolveReal> &newBubbleRestVolumes,
			    SIM_Object *obj,
			    const UT_Array<UT_Array<bool>> &oldToNewBubbleMap,
			    const UT_Array<bool> &areBubblesTracked,
			    const UT_Array<exint> &newBubbleRegionSizes,
			    const exint bubbleRegionCount);

    void
    printCellGeometry(SIM_Object *obj,
			const SIM_RawIndexField &activeCellIndices,
			const SIM_RawIndexField &materialCellLabels,
			const UT_Vector3 origin,
			const fpreal dx);

    void
    buildRHS(Vector &rhsVector,
		const SIM_RawIndexField &materialCellLabels,
		const SIM_RawIndexField &activeCellIndices,
		const SIM_VectorField &velocity,
		const std::array<const SIM_RawField *, 3> &cutCellWeights,
		const SIM_VectorField *solidVelocity,
		const SolveReal dx) const;

    void
    applyGoalDivergence(Vector &rhsVector,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &activeCellIndices,
			const SIM_RawField &goalDivergence) const;

    //
    // Helper methods to remove null space
    //

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
    removeBackgroundNullSpace(Vector &rhsVector,
				const SIM_RawIndexField &activeCellIndices,
				const SIM_RawIndexField &activeRegionIndices,
				const exint activeRegionCount) const;

    void
    setCorrectedBubbleDivergence(UT_Array<UT_Array<SolveReal>> &parallelCorrectedBubbleDivergence,
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
				    const bool doCorrectVolumeDrift) const;

    void
    removeCorrectionDivergenceNullSpace(Vector &rhsVector,
					const SIM_RawIndexField &materialCellLabels,
					const SIM_RawIndexField &activeCellIndices,
					const SIM_RawIndexField &activeRegionIndices,
					const SIM_RawIndexField &bubbleRegionIndices,
					const UT_Array<bool> &areBubblesTracked,
					const UT_Array<SolveReal> &activeRegionAverageDivergence,
					const bool doCorrectVolumeDrift) const;

    void
    applyBubbleCorrection(Vector &rhsVector,
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
			    const exint bubbleRegionCount) const;

    void
    buildLinearSystem(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelPoissonElements,
			std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelDiagonalPrecondElements,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &activeCellIndices,
			const std::array<const SIM_RawField *, 3> &cutCellWeights,
			const SIM_RawField &liquidSurface,
			const SolveReal liquidDensity,
			const SolveReal bubbleDensity) const;

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
			    const SIM_RawIndexField &activeCellIndices,
			    const SIM_RawField &pressure,
			    const SIM_RawField &liquidSurface,
			    const SIM_RawField &cutCellWeights,
			    const SolveReal liquidDensity,
			    const SolveReal bubbleDensity,
			    const int axis,
			    const SolveReal dx) const;

    void
    computeResultingDivergence(std::array<SolveReal, 3> &divergenceStats,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_VectorField &velocity,
				const std::array<const SIM_RawField *, 3> &cutCellWeights,
				const SIM_VectorField *solidVelocity) const;
};

#endif