#pragma once


#include "./TetMesh.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./Fiber.h"
#include "./io.h"
#include "src/Arrangement.h"


namespace augmentation
{
    std::map<int, std::vector<int>> computeRegularVerticesSheets(TetMesh &, Arrangement &, ReebSpace2 &); 
}
