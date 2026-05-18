#pragma once


#include "./TetMesh.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./Fiber.h"
#include "./io.h"
#include "src/Arrangement.h"


namespace augmentation
{
    std::vector<int> computeRegularVertexSheets(TetMesh &, Arrangement &, ReebSpace2 &); 
}
