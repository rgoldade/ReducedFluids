#ifndef HDK_CONSTRAINT_UTILITIES_H
#define HDK_CONSTRAINT_UTILITIES_H

#include <Eigen/Sparse>

#include <SIM/SIM_FieldUtils.h>
#include <SIM/SIM_RawIndexField.h>

#include <UT/UT_ParallelUtil.h>

#include "../External/GeometricMultigridPressureSolver/Source/HDK_GeometricMultigridOperators.h"
#include "../External/GeometricMultigridPressureSolver/Source/HDK_Utilities.h"

class SIM_VectorField;
class SIM_ScalarField;

namespace HDK::Utilities
{
    //
    // Common material label builder for dense and constraint bubbles
    //

    void
    overwriteIndices(SIM_RawIndexField &indices,
			const exint searchValue,
			const exint replaceValue);

    template<typename Vector>
    void
    applyInitialGuess(Vector &solutionVector,
			const SIM_RawField &pressure,
			const SIM_RawIndexField &activeCellIndices)
    {
	using SIM::FieldUtils::getFieldValue;

	UT_Interrupt *boss = UTgetInterrupt();

	UTparallelForEachNumber(activeCellIndices.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
	{
	    UT_VoxelArrayIteratorI vit;
	    vit.setConstArray(activeCellIndices.field());

	    if (boss->opInterrupt())
		return;    

	    UT_VoxelTileIteratorI vitt;

	    for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	    {
		vit.myTileStart = tileNumber;
		vit.myTileEnd = tileNumber + 1;
		vit.rewind();

		if (!vit.atEnd())
		{
		    if (!vit.isTileConstant())
		    {
			vitt.setTile(vit);

			for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
			{
			    exint activeIndex = vitt.getValue();
			    if (activeIndex >= 0)
			    {
				UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
				solutionVector(activeIndex) = getFieldValue(pressure, cell);
			    }
			}
		    }
		}
	    }
	});
    }
} // namespace HDK::Utilities
#endif