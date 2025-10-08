#pragma once

#include "./CGALTypedefs.h"

#include "./ReebSpace.h"
#include "./ReebSpace2.h"
#include "./Arrangement.h"
#include "./TetMesh.h"

#include <map>


namespace unitTests
{

    inline bool testAreSheetsIdentical(const TetMesh &tetMesh, const Arrangement &regularArrangement, const Arrangement &singularArrangement, ReebSpace &reebSpace, ReebSpace2 &reebSpaceS)
    {

        if (reebSpace.correspondenceGraph.getComponentRepresentatives().size() != reebSpaceS.numberOfSheets)
        {
            return false;
        }

        // Determine which points are singular, we don't need the regular ones
        std::set<Point_2> isPointSingular;
        for (const auto &[edge, singularType] : tetMesh.edgeSingularTypes)
        {
            if (singularType == 2 || singularType == 0)
            {
                isPointSingular.insert(singularArrangement.arrangementPoints[edge[0]]);
                isPointSingular.insert(singularArrangement.arrangementPoints[edge[1]]);
            }
        }


        std::map<int, std::set<Point_2>> pointsPerSheetR;

        for (auto v = regularArrangement.arr.vertices_begin(); v != regularArrangement.arr.vertices_end(); ++v)
        {
            // We only want the singular points
            if (false == isPointSingular.contains(v->point()))
            {
                continue;
            }

            // Get all incident face
            Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
            first = curr = v->incident_halfedges();
            do 
            {
                if (curr->face()->is_unbounded())
                {
                    continue;
                }

                const int currentFaceID = regularArrangement.arrangementFacesIdices.at(curr->face());

                for (const auto &[componentId, triangleId] : reebSpace.preimageGraphs[currentFaceID].getComponentRepresentatives())
                {
                    const int sheetId = reebSpace.correspondenceGraph.findElement({currentFaceID, componentId});

                    pointsPerSheetR[sheetId].insert(v->point());
                }

                curr->face();
            } while (++curr != first);

            //std::cout << std::endl;
        }

        //printf("\n\n--------------------------------------------\n");
        //printf("Printing points per sheet in AT\n");
        //printf("--------------------------------------------\n");
        //for (const auto &[sheetId, points] : pointsPerSheetR)
        //{
            //printf("Sheet with ID = %d has these ponts:\n", sheetId);

            //for (const auto &point : points)
            //{
                //std::cout << point << std::endl;
            //}
        //}



        std::map<int, std::set<Point_2>> pointsPerSheetS;

        for (auto v = singularArrangement.arr.vertices_begin(); v != singularArrangement.arr.vertices_end(); ++v)
        {
            // Get all incident face
            Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
            first = curr = v->incident_halfedges();
            do 
            {
                if (curr->face()->is_unbounded())
                {
                    continue;
                }

                const int currentFaceID = curr->face()->data();

                for (const auto &componentId : reebSpaceS.correspondenceGraph[currentFaceID])
                {
                    const int sheetId = reebSpaceS.correspondenceGraphDS.find(componentId);

                    pointsPerSheetS[sheetId].insert(v->point());
                }

                curr->face();
            } while (++curr != first);

            //std::cout << std::endl;
        }

        //printf("\n\n--------------------------------------------\n");
        //printf("Printing points per sheet in SAT\n");
        //printf("--------------------------------------------\n");
        //for (const auto &[sheetId, points] : pointsPerSheetS)
        //{
            //printf("Sheet with ID = %d has these ponts:\n", sheetId);

            //for (const auto &point : points)
            //{
                //std::cout << point << std::endl;
            //}
        //}




        std::map<int, int> matching;

        // To make sure the values are matched only once
        std::set<int> matchedS;

        for (const auto &[sheetIdR, pointsR] : pointsPerSheetR)
        {
            bool foundMatch = false;

            for (const auto &[sheetIdS, pointsS] : pointsPerSheetS)
            {
                if (pointsR == pointsS)
                {

                    //printf("\n\n--------------------------------------------\n");
                    //printf("Here's the current sheet matching:\n");
                    //printf("--------------------------------------------\n");
                    //for (int i = 0 ; i < matching.size() ; i++)
                    //for (const auto &[sheetIdR, sheetIdS] : matching)
                    //{
                        //printf("Pair : [%d, %d]\n", sheetIdR, sheetIdS);
                    //}

                    //printf("------------ Trying to match sheet %d with sheet %d\n", sheetIdR, sheetIdS);

                    // Something went wrong a sheet should only be matched once
                    if (matching.contains(sheetIdR) || matchedS.contains(sheetIdS))
                    {
                        //printf("Sheet already matched!\n");
                        return false;
                    }

                    matching[sheetIdR] = sheetIdS;
                    matchedS.insert(sheetIdS);
                    //sheetMatching.insert({sheetIdR, sheetIdS});
                    foundMatch = true;
                    break;
                }
            }
        }

        //printf("\n\n--------------------------------------------\n");
        //printf("Here's the sheet matching:\n");
        //printf("--------------------------------------------\n");
        ////for (int i = 0 ; i < matching.size() ; i++)
        //for (const auto &[sheetIdR, sheetIdS] : matching)
        //{
            //printf("Pair : [%d, %d]\n", sheetIdR, sheetIdS);
        //}

        return true;
    }

}
