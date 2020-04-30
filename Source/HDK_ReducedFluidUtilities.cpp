#include "HDK_ReducedFluidUtilities.h"

namespace HDK::Utilities
{
    void
    overwriteIndices(SIM_RawIndexField &indices,
			const exint searchValue,
			const exint replaceValue)
    {
	UT_Interrupt *boss = UTgetInterrupt();
	UTparallelForEachNumber(indices.field()->numTiles(), [&](const UT_BlockedRange<int> &range)
	{
	    UT_VoxelArrayIteratorI vit(indices.fieldNC());

	    if (boss->opInterrupt())
		return;    

	    for (int tileNumber = range.begin(); tileNumber != range.end(); ++tileNumber)
	    {
		vit.myTileStart = tileNumber;
		vit.myTileEnd = tileNumber + 1;
		vit.rewind();

		if (vit.isTileConstant() && vit.getValue() == searchValue)
		{
		    UT_VoxelTile<exint> *tile = vit.getTile();
		    tile->makeConstant(replaceValue);
		}
		else
		{
		    for (; !vit.atEnd(); vit.advance())
		    {
			if (vit.getValue() == searchValue)
			    vit.setValue(replaceValue);
		    }
		}
	    }
	});
    }
} // namespace HDK::Utilities