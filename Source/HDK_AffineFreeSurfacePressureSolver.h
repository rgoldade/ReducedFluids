#ifndef HDK_AFFINE_FREE_SURFACE_PRESSURE_SOLVER_H
#define HDK_AFFINE_FREE_SURFACE_PRESSURE_SOLVER_H

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>

#include <SIM/SIM_RawIndexField.h>

#include <UT/UT_ThreadedAlgorithm.h>

#include "HDK_ReducedFluidUtilities.h"

class SIM_VectorField;
class SIM_ScalarField;

class GAS_API HDK_AffineFreeSurfacePressureSolver : public GAS_SubSolver
{
    using SolveReal = double;
    using Vector = Eigen::VectorXd;
    using UT_Vector3SR = UT_Vector3T<SolveReal>;

    enum MaterialLabels { SOLID_CELL, AIR_CELL, INTERIOR_LIQUID_CELL, ACTIVE_LIQUID_CELL };
    
    using GradientMatrix = Eigen::Matrix<SolveReal, 3, 3>;
    using MassMatrix = Eigen::Matrix<SolveReal, 11, 11>;
    using ColumnVector = Eigen::Matrix<SolveReal, 11, 1>;

    static constexpr exint VISITED_CELL = -1;
    static constexpr exint UNVISITED_CELL = -2;

public:

    GET_DATA_FUNC_F(SIM_NAME_TOLERANCE, SolverTolerance);
    GET_DATA_FUNC_I("maxSolverIterations", MaxSolverIterations);

    GET_DATA_FUNC_B("useOldPressure", UseOldPressure);
    GET_DATA_FUNC_B("extrapolatePressure", DoExtrapolatePressure);

    GET_DATA_FUNC_I("boundaryLayerSize", BoundaryLayerSize);

    GET_DATA_FUNC_I("tileSize", TileSize);
    GET_DATA_FUNC_B("useTiledInterior", UseTiledInterior);

    GET_DATA_FUNC_B("printCells", PrintCells);
    GET_DATA_FUNC_B("onlyPrintCells", OnlyPrintCells);

protected:

    explicit HDK_AffineFreeSurfacePressureSolver(const SIM_DataFactory *factory);
    virtual ~HDK_AffineFreeSurfacePressureSolver();

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
    DECLARE_DATAFACTORY(HDK_AffineFreeSurfacePressureSolver,
                        GAS_SubSolver,
                        "HDK Affine Free Surface Pressure Solver",
                        getDopDescription());

    ////////////////////////////////////////////
    //
    // Solver methods
    //
    ////////////////////////////////////////////

    void
    setMaterialCellLabels(SIM_RawIndexField &materialCellLabels) const;

    void
    buildInitialAirBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelActiveCellList,
				const SIM_RawIndexField &materialCellLabels,
				const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    setActiveLiquidLayerCells(SIM_RawIndexField &materialCellLabels,
				const UT_Array<UT_Vector3I> &activeCellLayer) const;

    void
    buildNextLiquidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelActiveCellList,
				    const UT_Array<UT_Vector3I> &activeCellList,
				    const SIM_RawIndexField &materialCellLabels,
				    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;
    void
    buildInitialSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelActiveCellList,
				    const SIM_RawIndexField &materialCellLabels,
				    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    setActiveSolidLayerCells(SIM_RawIndexField &materialCellLabels,
				SIM_RawIndexField &visitedCells,
				const UT_Array<UT_Vector3I> &activeCellLayer) const;

    void
    buildNextSolidBoundaryLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelActiveCellList,
				const UT_Array<UT_Vector3I> &oldActiveCellList,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_RawIndexField &visitedCells,
				const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    buildTiledActiveCells(SIM_RawIndexField &materialCellLabels,
			    const int tileSize) const;

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
    buildActiveCellIndices(SIM_RawIndexField &activeCellLabels,
			    const SIM_RawIndexField &materialCellLabels) const;

    void
    buildInteriorRegionCOM(UT_Array<UT_Array<UT_Vector3T<SolveReal>>> &parallelInteriorRegionCOM,
			    UT_Array<UT_Array<SolveReal>> &parallelInteriorRegionCellCount,
			    const SIM_RawIndexField &materialCellLabels,
			    const SIM_RawIndexField &interiorRegionIndices,
			    const SIM_RawIndexField &activeRegionIndices) const;

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
		const UT_Array<UT_Vector3T<SolveReal>> &interiorLinearVelocities,
		const UT_Array<GradientMatrix> &interiorVelocityGradients,
		const UT_Array<UT_Vector3T<SolveReal>> &interiorRegionCOM,
		const SolveReal dx) const;

    void
    buildLinearSystem(std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelPoissonElements,
			std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelDiagonalPrecondElements,
			std::vector<std::vector<Eigen::Triplet<SolveReal>>> &parallelAffineCollectionElements,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &interiorRegionIndices,
			const SIM_RawIndexField &activeCellIndices,
			const std::array<const SIM_RawField *, 3> &cutCellWeights,
			const SIM_RawField &liquidSurface,
			const UT_Array<UT_Vector3T<SolveReal>> &interiorRegionCOM,
			const SolveReal dx) const;

    void
    buildAffineMassMatrixSystems(UT_Array<UT_Array<MassMatrix>> &parallelAffineMassMatrixElements,
				    const SIM_RawIndexField &materialCellLabels,
				    const SIM_RawIndexField &interiorRegionIndices,
				    const UT_Array<UT_Vector3T<SolveReal>> &interiorRegionCOM,
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

    void
    buildValidFaces(SIM_VectorField &validFacesField,
		    const SIM_RawIndexField &liquidCellLabels,
		    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    copyPressureVectorToGrid(SIM_RawField &pressure,
				const SIM_RawIndexField &activeCellIndices,
				const Vector &solutionVector) const;

    //
    // Helper methods to extrapolate pressure into interior regions
    //

    void
    buildInitialExtrapolationLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExtrapolationCellLayer,
				    SIM_RawIndexField &visitedCells,
				    const SIM_RawIndexField &materialCellLabels,
				    const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    extrapolatePressure(SIM_RawField &pressure,
			const SIM_RawIndexField &materialCellLabels,
			const SIM_RawIndexField &visitedCells,
			const std::array<const SIM_RawField *, 3> &cutCellWeights,
			const UT_Array<UT_Vector3I> &extrapolationCellLayer) const;

    void
    setVisitedExtrapolationCells(SIM_RawIndexField &visitedCells,
				    const SIM_RawIndexField &materialCellLabels,
				    const UT_Array<UT_Vector3I> &extrapolationCellLayer) const;

    void
    buildNextExtrapolationLayer(UT_Array<UT_Array<UT_Vector3I>> &parallelExtrapolationCellLayer,
				const UT_Array<UT_Vector3I> &oldExtrapolationCellLayer,
				const SIM_RawIndexField &materialCellLabels,
				const SIM_RawIndexField &visitedCells,
				const std::array<const SIM_RawField *, 3> &cutCellWeights) const;

    void
    applySolutionToVelocity(SIM_RawField &velocity,
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
			    const SolveReal dx) const;

    //
    // Debug check
    //

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
