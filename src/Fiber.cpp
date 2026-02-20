#include "./CGALTypedefs.h"

#include "./Fiber.h"
#include "./DisjointSet.h"
#include "./FiberPoint.h"
#include "./FiberGraph.h"
#include "./Timer.h"
#include "src/ReebSpace2.h"
#include "src/TetMesh.h"

#include <queue>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>



std::tuple<std::vector<int>, std::vector<int>, bool> getActiveTrianglesInPath(const std::vector<int> &pathA, const std::vector<int> &pathB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet)
{
    //std::cerr << "-------------------------------------------------\n";

    //std::cerr << "\n\nPath A: " << std::endl;

    //for (const int &triangleId : pathA)
    //{
        //std::cerr << triangleId << " ";
    //}

    //std::cerr << "\nPath B: " << std::endl;

    //for (const int &triangleId : pathB)
    //{
        //std::cerr << triangleId << " ";
    //}

    //std::cerr << "\n\nMinus triangles: " << std::endl;

    //for (const int &triangleId : minusTrianglesSet)
    //{
        //std::cerr << triangleId << std::endl;
    //}

    //std::cerr << "\nPlus triangles: " << std::endl;

    //for (const int &triangleId : plusTrianglesSet)
    //{
        //std::cerr << triangleId << std::endl;
    //}


    bool isSecondPathFlipped = false;


    std::vector<int> pathActiveTrianglesA;
    for (int i = 0 ; i < pathA.size() ; i++)
    {
        if (minusTrianglesSet.contains(pathA[i]))
        {
            pathActiveTrianglesA.push_back(i);
        }
    }

    std::vector<int> pathActiveTrianglesB;
    for (int i = 0 ; i < pathB.size() ; i++)
    {
        if (plusTrianglesSet.contains(pathB[i]))
        {
            pathActiveTrianglesB.push_back(i);
        }
    }

    int prevA = pathActiveTrianglesA[0] - 1;
    int nextA = pathActiveTrianglesA.back() + 1;

    int prevB = pathActiveTrianglesB[0] - 1;
    int nextB = pathActiveTrianglesB.back() + 1;

    // If the substring starts at the left endpint
    if (prevA == -1)
    {
        // But at the right endpoint of B
        if (prevB != -1)
        {
            isSecondPathFlipped = true;
        }

    }
    // If the substring starts at the right endpoint
    else if (prevA == pathA.size())
    {
        // But at the left endpoint of B
        if (prevB != pathB.size())
        {
            isSecondPathFlipped = true;
        }
    }
    // If the substring is in the middle
    else 
    {
        if (pathA[prevA] != pathB[prevB])
        {
            isSecondPathFlipped = true;
        }

    }

    if (isSecondPathFlipped)
    {
        //std::cerr << "We reversed!\n";
        std::reverse(pathActiveTrianglesB.begin(), pathActiveTrianglesB.end());
    }

    //std::cerr << "\n\nSubstring A: ";

    //for (const int &index : pathActiveTrianglesA)
    //{
        //std::cerr << pathA[index] << " ";
    //}

    //std::cerr << "\nSubstring B: ";

    //for (const int &index : pathActiveTrianglesB)
    //{
        //std::cerr << pathB[index] << " ";
    //}

    //std::cerr << std::endl;

    return {pathActiveTrianglesA, pathActiveTrianglesB, isSecondPathFlipped};
}


std::tuple<std::vector<int>, std::vector<int>, bool> getActiveTrianglesInCycle(const std::vector<int> &cycleA, const std::vector<int> &cycleB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet)
{
    //std::cerr << "It's a cycle " << std::endl;

    bool isSecondCycleFlipped = false;

    const int start = 0;
    int current = start;
    int startingPointA = -1;


    // Find the starting point
    //
    do
    {
        const int next = (current + 1) % cycleA.size();

        if (
                false  == minusTrianglesSet.contains(cycleA[current]) &&
                true == minusTrianglesSet.contains(cycleA[next])
           )
        {
            startingPointA = current;
            break;
        }

        current = next;
    } while (current != start);

    if (startingPointA == -1)
    {
        throw std::runtime_error("startingPointA not found in cycleA");
    }

    //std::cerr << "Starting point is " << cycleA[startingPointA] << std::endl;


    // Add the affected triangles
    //
    std::vector<int> cycleActiveTrianglesA;
    current = (startingPointA + 1) % cycleA.size();
    do
    {
        cycleActiveTrianglesA.push_back(current);
        const int next = (current + 1) % cycleA.size();
        current = next;

    } while (minusTrianglesSet.contains(cycleA[current]));





    // Find the starting point in cycleB
    current = 0;
    int startingPointB = -1;
    do
    {

        if (cycleB[current] == cycleA[startingPointA])
        {
            startingPointB = current;
            break;
        }

        const int next = (current + 1) % cycleB.size();
        current = next;
    } while (current != start);

    if (startingPointB == -1)
    {
        throw std::runtime_error("startingPointB not found in cycleB");
    }


    const int sizeA = static_cast<int>(cycleA.size());
    const int sizeB = static_cast<int>(cycleB.size());

    const int previousA = (startingPointA == 0) ? sizeA - 1 : startingPointA - 1;
    const int nextA = (startingPointA == sizeA - 1) ? 0 : startingPointA + 1;

    const int previousB = (startingPointB == 0) ? sizeB - 1 : startingPointB - 1;
    const int nextB = (startingPointB == sizeB - 1) ? 0 : startingPointB + 1;


    std::vector<int> cycleActiveTrianglesB;

    // Already have the same orientation
    //
    if (cycleA[previousA] == cycleB[previousB])
    {
        current = (startingPointB + 1) % cycleB.size();
        do
        {
            cycleActiveTrianglesB.push_back(current);
            const int next = (current + 1) % cycleB.size();
            current = next;

        } while (plusTrianglesSet.contains(cycleB[current]));
    }

    // Reverse orientation
    //
    else if (cycleA[previousA] == cycleB[nextB])
    {
        isSecondCycleFlipped = true;

        current = (startingPointB == 0) ? sizeB - 1 : startingPointB - 1;
        do
        {
            cycleActiveTrianglesB.push_back(current);
            const int previous = (current == 0) ? sizeB - 1 : current - 1;
            current = previous;

        } while (plusTrianglesSet.contains(cycleB[current]));


        //std::cerr << "\n\n\n-----------------------------------------------------------------\n";
        //std::cerr << "Cycle A : ";
        //for (const int &p : cycleA)
        //{
            //std::cerr << p << " ";

        //}

        //std::cerr << "\nCycle B : ";
        //for (const int &p : cycleB)
        //{
            //std::cerr << p << " ";
        //}

        //std::cerr << "\n\nThe active triangles in cycle A are : ";
        //for (const int &t : cycleActiveTrianglesA)
        //{
            //std::cerr << cycleA[t] << "  ";
        //}
        //std::cerr << "\nThe active triangles in cycle B are : ";
        //for (const int &t : cycleActiveTrianglesB)
        //{
            //std::cerr << cycleB[t] << "  ";
        //}

        //std::cerr << "\n\n\n";

        //std::cerr << "StartingPointA = " << cycleA[startingPointA] << " StartingPointB = " << cycleB[startingPointB] << " \n ";
        //std::cerr << "Next A = " << cycleA[nextA] << " Next B = " << cycleB[startingPointB] << " \n ";
        //std::cerr << "Prev A = " << cycleA[previousA] << " Prev B = " << cycleB[previousB] << " \n ";
        //std::cerr << "Cycle A : \n";
        //for (const int &p : cycleA)
        //{
        //std::cerr << p << " ";

        //}

        //std::cerr << "\n Cycle B : \n";
        //for (const int &p : cycleB)
        //{
        //std::cerr << p << " ";
        //}

        //throw std::runtime_error("Cycles orientation reversed.");
    }

    else
    {
        //std::cerr << "StartingPointA = " << cycleA[startingPointA] << " StartingPointB = " << cycleB[startingPointB] << " \n ";
        //std::cerr << "Next A = " << cycleA[nextA] << " Next B = " << cycleB[startingPointB] << " \n ";
        //std::cerr << "Prev A = " << cycleA[previousA] << " Prev B = " << cycleB[previousB] << " \n ";
        //std::cerr << "Cycle A : \n";
        //for (const int &p : cycleA)
        //{
        //std::cerr << p << " ";

        //}

        //std::cerr << "\n Cycle B : \n";
        //for (const int &p : cycleB)
        //{
        //std::cerr << p << " ";
        //}

        throw std::runtime_error("Cycles could not be orientted properly.");
    }



    //std::vector<int> cycleActiveTrianglesAReturn(cycleActiveTrianglesA.size());
    //for (int i = 0 ; i < cycleActiveTrianglesA.size() ; i++)
    //{
        //cycleActiveTrianglesAReturn[i] = cycleA[cycleActiveTrianglesA[i]];
    //}

    //std::vector<int> cycleActiveTrianglesBReturn(cycleActiveTrianglesB.size());
    //for (int i = 0 ; i < cycleActiveTrianglesB.size() ; i++)
    //{
        //cycleActiveTrianglesBReturn[i] = cycleB[cycleActiveTrianglesB[i]];
    //}

    //return {cycleActiveTrianglesAReturn, cycleActiveTrianglesBReturn};


    return {cycleActiveTrianglesA, cycleActiveTrianglesB, isSecondCycleFlipped};
}

std::vector<FiberPoint> computeTrianglesBetweenCorrespondingPaths(const std::vector<int> &pathA, const std::vector<int> &pathB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet, const std::vector<std::unordered_map<int, std::array<double, 3>>> &barycentricCoordinatesPerTriangle, const TetMesh &tetMesh, const int sheetId, const std::array<float, 3> &sheetColour, const int i, const std::array<float, 3> &edgePointDomain)
{
    std::vector<FiberPoint> allFiberPoints;

    const auto [pathActiveTrianglesA, pathActiveTrianglesB, isSecondPathFlipped] = getActiveTrianglesInPath(pathA, pathB, minusTrianglesSet, plusTrianglesSet);


    //
    //      Connect cycleA
    //
    //      a
    //      |\
    //      | \
    //      |  \
    //      b---o
    //      .  .
    //      . .
    //      ..
    //      c
    //      
    //


    for (int k = 0 ; k < pathActiveTrianglesA.size() - 1 ; k++)
    {
        std::vector<FiberPoint> componentFiberPoints;

        const int triangleIdA = pathA[pathActiveTrianglesA[k]];
        const int triangleIdB = pathA[pathActiveTrianglesA[k + 1]];

        const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i-1].at(triangleIdA);
        const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i-1].at(triangleIdB);

        //fprintf(stderr, "Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

        componentFiberPoints.push_back(FiberPoint(
                    a[0], 
                    a[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                    sheetColour,
                    sheetId,
                    triangleIdA
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    b[0], 
                    b[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                    sheetColour,
                    sheetId,
                    triangleIdB
                    ));


        componentFiberPoints.push_back(FiberPoint(
                    edgePointDomain,
                    sheetColour,
                    sheetId,
                    -1
                    ));


        allFiberPoints.insert(
                allFiberPoints.end(), 
                std::make_move_iterator(componentFiberPoints.begin()), 
                std::make_move_iterator(componentFiberPoints.end())
                );

    }




    //
    //
    //      Connect cycleB
    //
    //            a
    //           /|
    //          / | 
    //         /  |  
    //        o---b
    //         .  .
    //          . .
    //           ..
    //            c
    //      
    //
    for (int k = 0 ; k < pathActiveTrianglesB.size() - 1 ; k++)
    {
        std::vector<FiberPoint> componentFiberPoints;

        const int triangleIdA = pathB[pathActiveTrianglesB[k]];
        const int triangleIdB = pathB[pathActiveTrianglesB[k + 1]];

        const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i].at(triangleIdA);
        const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i].at(triangleIdB);

        //fprintf(stderr, "Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

        componentFiberPoints.push_back(FiberPoint(
                    a[0], 
                    a[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                    sheetColour,
                    sheetId,
                    triangleIdA
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    edgePointDomain,
                    sheetColour,
                    sheetId,
                    -1
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    b[0], 
                    b[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                    sheetColour,
                    sheetId,
                    triangleIdB
                    ));

        allFiberPoints.insert(
                allFiberPoints.end(), 
                std::make_move_iterator(componentFiberPoints.begin()), 
                std::make_move_iterator(componentFiberPoints.end())
                );
    }









    // The last two triangles (cap)
    //
    //      a-------b
    //       \     /.
    //        \   /  
    //         \ /  
    //          o
    //         / \
    //        /   \
    //       /     \
    //      c-------d


    const int triangleIdA = pathA[pathActiveTrianglesA[0]];
    const int triangleIdB = pathB[pathActiveTrianglesB[0]];

    const int triangleIdC = pathA[pathActiveTrianglesA.back()];
    const int triangleIdD = pathB[pathActiveTrianglesB.back()];

    const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i-1].at(triangleIdA);
    const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i].at(triangleIdB);
    const std::array<double, 3> &c = barycentricCoordinatesPerTriangle[i-1].at(triangleIdC);
    const std::array<double, 3> &d = barycentricCoordinatesPerTriangle[i].at(triangleIdD);



    std::vector<FiberPoint> componentFiberPoints;

    //fprintf(stderr, "Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

    // Add aob
    componentFiberPoints.push_back(FiberPoint(
                a[0], 
                a[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                sheetColour,
                sheetId,
                triangleIdA
                ));


    componentFiberPoints.push_back(FiberPoint(
                edgePointDomain,
                sheetColour,
                sheetId,
                -1
                ));

    componentFiberPoints.push_back(FiberPoint(
                b[0], 
                b[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                sheetColour,
                sheetId,
                triangleIdB
                ));


    //fprintf(stderr, "Adding triangle %d %d %d\n", triangleIdC, triangleIdD, -1);

    // Add cdo
    componentFiberPoints.push_back(FiberPoint(
                c[0], 
                c[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdC),
                sheetColour,
                sheetId,
                triangleIdC
                ));

    componentFiberPoints.push_back(FiberPoint(
                d[0], 
                d[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdD),
                sheetColour,
                sheetId,
                triangleIdD
                ));


    componentFiberPoints.push_back(FiberPoint(
                edgePointDomain,
                sheetColour,
                sheetId,
                -1
                ));


    //      u------v
    //      |     /|
    //      |    / |
    //      |   /  |
    //      |  /   |
    //      | /    |
    //      |/     |
    //      a------b
    //         .
    //         .
    //         .
    //      c------d
    //      |     /|
    //      |    / |
    //      |   /  |
    //      |  /   |
    //      | /    |
    //      |/     |
    //      w------z


    // @TODO I need to know if the direction of the second cycle is flipped!
    //
    const int index_a = pathActiveTrianglesA[0];
    const int index_b = pathActiveTrianglesB[0];

    const int index_c = pathActiveTrianglesA.back();
    const int index_d = pathActiveTrianglesB.back();

    const int sizeA = static_cast<int>(pathA.size());
    const int sizeB = static_cast<int>(pathB.size());


    const int previous_a = index_a - 1;
    const int next_c = index_c + 1;

    int previous_b, next_d;
    if (false == isSecondPathFlipped)
    {
        previous_b = index_b - 1;
        next_d = index_d + 1;
    }
    else
    {
        previous_b = index_b + 1;
        next_d = index_d - 1;
    }


    if (previous_a != -1 && previous_b != -1)
    {
        const int triangleIdU = pathA[previous_a];
        const int triangleIdV = pathB[previous_b];

        const std::array<double, 3> &u = barycentricCoordinatesPerTriangle[i-1].at(triangleIdU);
        const std::array<double, 3> &v = barycentricCoordinatesPerTriangle[i].at(triangleIdV);

        //fprintf(stderr, "avu - Adding triangle %d %d A) %d B)\n", triangleIdA, triangleIdV, triangleIdU);

        //avu
        componentFiberPoints.push_back(FiberPoint(
                    a[0], 
                    a[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                    sheetColour,
                    sheetId,
                    triangleIdA
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    v[0], 
                    v[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdV),
                    sheetColour,
                    sheetId,
                    triangleIdV
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    u[0], 
                    u[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdU),
                    sheetColour,
                    sheetId,
                    triangleIdU
                    ));


        //fprintf(stderr, "abv - Adding triangle %d %d %d\n", triangleIdA, triangleIdB, triangleIdV);

        //abv
        componentFiberPoints.push_back(FiberPoint(
                    a[0], 
                    a[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                    sheetColour,
                    sheetId,
                    triangleIdA
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    b[0], 
                    b[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                    sheetColour,
                    sheetId,
                    triangleIdB
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    v[0], 
                    v[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdV),
                    sheetColour,
                    sheetId,
                    triangleIdV
                    ));

    }

    if (next_c != pathA.size() && next_d != pathB.size())
    {
        const int triangleIdW = pathA[next_c];
        const int triangleIdZ = pathB[next_d];

        const std::array<double, 3> &w = barycentricCoordinatesPerTriangle[i-1].at(triangleIdW);
        const std::array<double, 3> &z = barycentricCoordinatesPerTriangle[i].at(triangleIdZ);




        //fprintf(stderr, "dcw - Adding triangle %d B) %d A) %d\n", triangleIdD, triangleIdC, triangleIdW);

        // dcw

        componentFiberPoints.push_back(FiberPoint(
                    d[0], 
                    d[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdD),
                    sheetColour,
                    sheetId,
                    triangleIdD
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    c[0], 
                    c[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdC),
                    sheetColour,
                    sheetId,
                    triangleIdC
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    w[0], 
                    w[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdW),
                    sheetColour,
                    sheetId,
                    triangleIdW
                    ));

        //fprintf(stderr, "dwz - Adding triangle %d %d %d\n", triangleIdD, triangleIdW, triangleIdZ);

        // dwz
        componentFiberPoints.push_back(FiberPoint(
                    d[0], 
                    d[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdD),
                    sheetColour,
                    sheetId,
                    triangleIdD
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    w[0], 
                    w[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdW),
                    sheetColour,
                    sheetId,
                    triangleIdW
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    z[0], 
                    z[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdZ),
                    sheetColour,
                    sheetId,
                    triangleIdZ
                    ));




    }


    allFiberPoints.insert(
            allFiberPoints.end(), 
            std::make_move_iterator(componentFiberPoints.begin()), 
            std::make_move_iterator(componentFiberPoints.end())
            );


    return allFiberPoints;
}

std::vector<FiberPoint> computeTrianglesBetweenCorrespondingCycles(const std::vector<int> &cycleA, const std::vector<int> &cycleB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet, const std::vector<std::unordered_map<int, std::array<double, 3>>> &barycentricCoordinatesPerTriangle, const TetMesh &tetMesh, const int sheetId, const std::array<float, 3> &sheetColour, const int i, const std::array<float, 3> &edgePointDomain)
{
    std::vector<FiberPoint> allFiberPoints;

    const auto [cycleActiveTrianglesA, cycleActiveTrianglesB, isSecondCycleFlipped] = getActiveTrianglesInCycle(
            cycleA, 
            cycleB,
            minusTrianglesSet, 
            plusTrianglesSet
            );



    //
    //      Connect cycleA
    //
    //      a
    //      |\
    //      | \
    //      |  \
    //      b---o
    //      .  .
    //      . .
    //      ..
    //      c
    //      
    //


    for (int k = 0 ; k < cycleActiveTrianglesA.size() - 1 ; k++)
    {
        std::vector<FiberPoint> componentFiberPoints;

        const int triangleIdA = cycleA[cycleActiveTrianglesA[k]];
        const int triangleIdB = cycleA[cycleActiveTrianglesA[k + 1]];

        const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i-1].at(triangleIdA);
        const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i-1].at(triangleIdB);

        //fprintf(stderr, "Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

        componentFiberPoints.push_back(FiberPoint(
                    a[0], 
                    a[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                    sheetColour,
                    sheetId,
                    triangleIdA
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    b[0], 
                    b[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                    sheetColour,
                    sheetId,
                    triangleIdB
                    ));


        componentFiberPoints.push_back(FiberPoint(
                    edgePointDomain,
                    sheetColour,
                    sheetId,
                    -1
                    ));


        allFiberPoints.insert(
                allFiberPoints.end(), 
                std::make_move_iterator(componentFiberPoints.begin()), 
                std::make_move_iterator(componentFiberPoints.end())
                );

    }




    //
    //
    //      Connect cycleB
    //
    //            a
    //           /|
    //          / | 
    //         /  |  
    //        o---b
    //         .  .
    //          . .
    //           ..
    //            c
    //      
    //
    for (int k = 0 ; k < cycleActiveTrianglesB.size() - 1 ; k++)
    {
        std::vector<FiberPoint> componentFiberPoints;

        const int triangleIdA = cycleB[cycleActiveTrianglesB[k]];
        const int triangleIdB = cycleB[cycleActiveTrianglesB[k + 1]];

        const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i].at(triangleIdA);
        const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i].at(triangleIdB);

        //fprintf(stderr, "Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

        componentFiberPoints.push_back(FiberPoint(
                    a[0], 
                    a[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                    sheetColour,
                    sheetId,
                    triangleIdA
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    edgePointDomain,
                    sheetColour,
                    sheetId,
                    -1
                    ));

        componentFiberPoints.push_back(FiberPoint(
                    b[0], 
                    b[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                    sheetColour,
                    sheetId,
                    triangleIdB
                    ));

        allFiberPoints.insert(
                allFiberPoints.end(), 
                std::make_move_iterator(componentFiberPoints.begin()), 
                std::make_move_iterator(componentFiberPoints.end())
                );
    }

    // The last two triangles (cap)
    //
    //      a-------b
    //       \     /.
    //        \   /  
    //         \ /  
    //          o
    //         / \
    //        /   \
    //       /     \
    //      c-------d


    const int triangleIdA = cycleA[cycleActiveTrianglesA[0]];
    const int triangleIdB = cycleB[cycleActiveTrianglesB[0]];

    const int triangleIdC = cycleA[cycleActiveTrianglesA.back()];
    const int triangleIdD = cycleB[cycleActiveTrianglesB.back()];

    const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i-1].at(triangleIdA);
    const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i].at(triangleIdB);
    const std::array<double, 3> &c = barycentricCoordinatesPerTriangle[i-1].at(triangleIdC);
    const std::array<double, 3> &d = barycentricCoordinatesPerTriangle[i].at(triangleIdD);



    std::vector<FiberPoint> componentFiberPoints;

    //fprintf(stderr, "Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

    // Add aob
    componentFiberPoints.push_back(FiberPoint(
                a[0], 
                a[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                sheetColour,
                sheetId,
                triangleIdA
                ));


    componentFiberPoints.push_back(FiberPoint(
                edgePointDomain,
                sheetColour,
                sheetId,
                -1
                ));

    componentFiberPoints.push_back(FiberPoint(
                b[0], 
                b[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                sheetColour,
                sheetId,
                triangleIdB
                ));


    //fprintf(stderr, "Adding triangle %d %d %d\n", triangleIdC, triangleIdD, -1);

    // Add cdo
    componentFiberPoints.push_back(FiberPoint(
                c[0], 
                c[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdC),
                sheetColour,
                sheetId,
                triangleIdC
                ));

    componentFiberPoints.push_back(FiberPoint(
                d[0], 
                d[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdD),
                sheetColour,
                sheetId,
                triangleIdD
                ));


    componentFiberPoints.push_back(FiberPoint(
                edgePointDomain,
                sheetColour,
                sheetId,
                -1
                ));









    //      u------v
    //      |     /|
    //      |    / |
    //      |   /  |
    //      |  /   |
    //      | /    |
    //      |/     |
    //      a------b
    //         .
    //         .
    //         .
    //      c------d
    //      |     /|
    //      |    / |
    //      |   /  |
    //      |  /   |
    //      | /    |
    //      |/     |
    //      w------z


    // @TODO I need to know if the direction of the second cycle is flipped!
    //
    const int index_a = cycleActiveTrianglesA[0];
    const int index_b = cycleActiveTrianglesB[0];

    const int index_c = cycleActiveTrianglesA.back();
    const int index_d = cycleActiveTrianglesB.back();

    const int sizeA = static_cast<int>(cycleA.size());
    const int sizeB = static_cast<int>(cycleB.size());


    const int previous_a = (index_a == 0) ? sizeA - 1 : index_a - 1;
    const int next_c = (index_c + 1) % sizeA;

    int previous_b, next_d;
    if (false == isSecondCycleFlipped)
    {
        previous_b = (index_b == 0) ? sizeB - 1 : index_b - 1;
        next_d = (index_d + 1) % sizeB;
    }
    else
    {
        previous_b = (index_b + 1) % sizeB;
        next_d = (index_d == 0) ? sizeB - 1 : index_d - 1;
    }


    const int triangleIdU = cycleA[previous_a];
    const int triangleIdV = cycleB[previous_b];

    const int triangleIdW = cycleA[next_c];
    const int triangleIdZ = cycleB[next_d];

    const std::array<double, 3> &u = barycentricCoordinatesPerTriangle[i-1].at(triangleIdU);
    const std::array<double, 3> &v = barycentricCoordinatesPerTriangle[i].at(triangleIdV);
    const std::array<double, 3> &w = barycentricCoordinatesPerTriangle[i-1].at(triangleIdW);
    const std::array<double, 3> &z = barycentricCoordinatesPerTriangle[i].at(triangleIdZ);



    //fprintf(stderr, "avu - Adding triangle %d %d A) %d B)\n", triangleIdA, triangleIdV, triangleIdU);

    //avu
    componentFiberPoints.push_back(FiberPoint(
                a[0], 
                a[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                sheetColour,
                sheetId,
                triangleIdA
                ));

    componentFiberPoints.push_back(FiberPoint(
                v[0], 
                v[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdV),
                sheetColour,
                sheetId,
                triangleIdV
                ));

    componentFiberPoints.push_back(FiberPoint(
                u[0], 
                u[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdU),
                sheetColour,
                sheetId,
                triangleIdU
                ));


    //fprintf(stderr, "abv - Adding triangle %d %d %d\n", triangleIdA, triangleIdB, triangleIdV);

    //abv
    componentFiberPoints.push_back(FiberPoint(
                a[0], 
                a[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                sheetColour,
                sheetId,
                triangleIdA
                ));

    componentFiberPoints.push_back(FiberPoint(
                b[0], 
                b[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                sheetColour,
                sheetId,
                triangleIdB
                ));

    componentFiberPoints.push_back(FiberPoint(
                v[0], 
                v[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdV),
                sheetColour,
                sheetId,
                triangleIdV
                ));

    //fprintf(stderr, "dcw - Adding triangle %d B) %d A) %d\n", triangleIdD, triangleIdC, triangleIdW);

    // dcw

    componentFiberPoints.push_back(FiberPoint(
                d[0], 
                d[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdD),
                sheetColour,
                sheetId,
                triangleIdD
                ));

    componentFiberPoints.push_back(FiberPoint(
                c[0], 
                c[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdC),
                sheetColour,
                sheetId,
                triangleIdC
                ));

    componentFiberPoints.push_back(FiberPoint(
                w[0], 
                w[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdW),
                sheetColour,
                sheetId,
                triangleIdW
                ));

    //fprintf(stderr, "dwz - Adding triangle %d %d %d\n", triangleIdD, triangleIdW, triangleIdZ);

    // dwz
    componentFiberPoints.push_back(FiberPoint(
                d[0], 
                d[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdD),
                sheetColour,
                sheetId,
                triangleIdD
                ));

    componentFiberPoints.push_back(FiberPoint(
                w[0], 
                w[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdW),
                sheetColour,
                sheetId,
                triangleIdW
                ));

    componentFiberPoints.push_back(FiberPoint(
                z[0], 
                z[1], 
                tetMesh.getTriangleVerticesCoordinates(triangleIdZ),
                sheetColour,
                sheetId,
                triangleIdZ
                ));


    allFiberPoints.insert(
            allFiberPoints.end(), 
            std::make_move_iterator(componentFiberPoints.begin()), 
            std::make_move_iterator(componentFiberPoints.end())
            );


    return allFiberPoints;
}



std::array<double, 3> computeBarycentricCoordinates(const TetMesh &tetMesh, const int &triangleId, const std::array<double, 2> &fiberPoint)
{
    // Unpack the actual numerical  coordinates of the vertices of the triangle
    std::vector<CartesianPoint> triangleVertexCoordinates;

    for (const int &vertexId : tetMesh.triangles[triangleId])
    {
        triangleVertexCoordinates.push_back(
                CartesianPoint(
                    tetMesh.vertexCoordinatesF[vertexId], 
                    tetMesh.vertexCoordinatesG[vertexId]
                    ));
    }

    std::array<double, 3> barycentricCoordinates;

    CartesianPoint P(fiberPoint[0], fiberPoint[1]);

    CGAL::Barycentric_coordinates::triangle_coordinates_2(
            triangleVertexCoordinates[0], 
            triangleVertexCoordinates[1], 
            triangleVertexCoordinates[2], 
            P, barycentricCoordinates.begin());

    return barycentricCoordinates;
}


std::vector<FiberPoint> fiber::computeFiberSurface(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::vector<std::array<double, 2>> &controlPoints)
{

    //const Point_2 startPoint(0.479988 , 0.215557);
    //const Point_2 endPoint(0.259815, 0.116635);

    //const Point_2 startPoint(0.215557, 0.479988);
    //const Point_2 endPoint(0.116635, 0.259815);

    const Point_2 startPoint(controlPoints[0][0], controlPoints[0][1]);
    const Point_2 endPoint(controlPoints[1][0], controlPoints[1][1]);


    const Segment_2 controlSegment(startPoint, endPoint);




    std::cout << "Start point : " << startPoint << " end point " << endPoint << std::endl;


    Timer::start();

    // Compute intersectinos with the AABB tree
    //
    Timer::start();
    std::vector<TreeAABB::Primitive_id> intersectedSegmentsAABB;
    singularArrangement.tree.all_intersected_primitives(controlSegment, std::back_inserter(intersectedSegmentsAABB));
    Timer::stop("Computed AABB intersections in         :");

    Timer::start();




    std::vector<std::pair<K::FT, int>> intersectedSegments;
    intersectedSegments.reserve(intersectedSegmentsAABB.size());

    if (intersectedSegmentsAABB.size() == 0)
    {
        std::cout << "No segments were intersectd!.\n";
        return {};
    }



    for (auto id : intersectedSegmentsAABB)
    {
        const Segment_2& s = *id;   // dereference iterator to get the original segment

        // Get the ID of the original segment.
        const int segmentIndex = id - singularArrangement.allSegments.begin();


        const int indexSource = singularArrangement.arrangementPointIndices[s.source()];
        const int indexTarget = singularArrangement.arrangementPointIndices[s.target()];

        const int edgeId = tetMesh.edgeIndices.at({indexSource, indexTarget});

        if (segmentIndex != edgeId)
        {
            throw std::runtime_error("Issue in AABB tree segments indices.");
        }

        // Double check we get the same segment back.
        Point_2 a = singularArrangement.arrangementPoints[tetMesh.edges[segmentIndex][0]];
        Point_2 b = singularArrangement.arrangementPoints[tetMesh.edges[segmentIndex][1]];

        const bool match =
            (s.source() == a && s.target() == b) ||
            (s.source() == b && s.target() == a);

        if (match == false)
        {
            throw std::runtime_error("Issue in AABB tree segments interation.");
        }

        // ---- Collinear check with your controlSegment ----
        if (CGAL::collinear(controlSegment.source(), controlSegment.target(), s.source()) &&
                CGAL::collinear(controlSegment.source(), controlSegment.target(), s.target()))
        {
            // The intersection segment is collinear with the control segment
            std::cerr << "------------------------------------------ Collinear overlap detected for segment " << segmentIndex << std::endl;
            continue;
        }


        //
        //
        //
        //                          s.source()
        //                              |
        //                              |
        //                              |
        // searchSegment.source() ------x---------- searchSegment.target()
        //                              |
        //                              |
        //                              |
        //                              |
        //                          s.target()
        //
        K::FT alpha = CGAL::Intersections::internal::s2s2_alpha(
                controlSegment.target().x(), controlSegment.target().y(),
                controlSegment.source().x(), controlSegment.source().y(),
                s.source().x(), s.source().y(),
                s.target().x(), s.target().y());

        intersectedSegments.emplace_back(alpha, segmentIndex);


        //std::cout << "---- Intersected segment with ID " << segmentIndex << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentIndex)) << " and alpha " << alpha << std::endl;
    }

    Timer::stop("Computed Alpha intersections           :");
    //Timer::stop("Computed AABB intersections in         :");


    Timer::start();
    std::sort(intersectedSegments.begin(), intersectedSegments.end());
    Timer::stop("Sorting alpha intersections            :");









    // Set up the orientations
    //

    std::vector<bool> typicalOrientation;
    typicalOrientation.reserve(intersectedSegmentsAABB.size());


    std::vector<double> intersectedSegmentsAlphas;
    intersectedSegmentsAlphas.reserve(intersectedSegmentsAABB.size());

    for (int i = 0 ; i <  intersectedSegments.size() ; i++)
    {
        const auto &[alpha1, edgeId] = intersectedSegments[i];

        // Change orientation in case we need to
        const std::array<int, 2> edge = tetMesh.edges.at(edgeId);
        const int segmentSourceId = edge[0];
        const int segmentTargetId = edge[1];
        const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];
        const Point_2 &d = singularArrangement.arrangementPoints[segmentTargetId];

        //std::cout << c << " - " << d << std::endl;

        //const Point_2 &a = startPoint;
        //const Point_2 &b = endPoint;
        //
        // If the orientation this way around, reverse it
        //   (regular segment)
        //          d 
        //           \
        //            \
        // a ----------\--------- b (singular segment)
        //              \
        //               \
        //                c
        //
        if (CGAL::orientation(startPoint, endPoint, c) == CGAL::RIGHT_TURN)
        {
            //typicalOrientation[i] = false;
            typicalOrientation.emplace_back(false);

            //std::cout << "Orienting " << startPoint << ", " << endPoint << ", " << c << std::endl;
        }
        else
        {
            typicalOrientation.emplace_back(true);
            //std::cout << "Orienting Else " << startPoint << ", " << endPoint << ", " << c << std::endl;
        }

        K::FT alpha = CGAL::Intersections::internal::s2s2_alpha(
                d.x(), d.y(),
                c.x(), c.y(),
                controlSegment.target().x(), controlSegment.target().y(),
                controlSegment.source().x(), controlSegment.source().y()
                );

        intersectedSegmentsAlphas.emplace_back(CGAL::to_double(alpha));
    }

















    // Fet up the fiber points in between the intersecting segments
    //
    Timer::start();
    std::vector<Point_2> fiberPoints;
    fiberPoints.reserve(intersectedSegments.size() + 1);


    fiberPoints.emplace_back(startPoint);


    for (int i = 1 ; i <  intersectedSegments.size() ; i++)
    {
        const auto &[alpha1, edgeId1] = intersectedSegments[i-1];
        const auto &[alpha2, edgeId2] = intersectedSegments[i];

        K::FT fiberPointAlpha = (alpha1 + alpha2) / 2.0;

        Point_2 fiberPoint = startPoint + fiberPointAlpha * (endPoint - startPoint);

        fiberPoints.emplace_back(fiberPoint);
    }

    fiberPoints.emplace_back(endPoint);
    Timer::stop("Finding the midpoint alphas            :");



    std::vector<FiberGraph> fiberGraphs(fiberPoints.size());
    std::vector<std::pair<std::map<int, std::vector<int>>, std::map<int, std::vector<int>>>> pathsAndCycles(fiberPoints.size());
    std::vector<std::unordered_map<int, std::array<double, 3>>> barycentricCoordinatesPerTriangle(fiberPoints.size());


    Timer::start();

    // Do the heavy lifting
    //
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0 ; i < fiberPoints.size() ; i++)
    {
        // Set up the fiber points
        const Point_2 &fiberPoint = fiberPoints[i];
        const std::array<double, 2> fiberPointDouble = { CGAL::to_double(fiberPoint.x()), CGAL::to_double(fiberPoint.y()) };

        // Compute the fiber graph and its paths and cycles
        fiberGraphs[i] = fiber::computeFiberGraph(tetMesh, singularArrangement, reebSpace, fiberPointDouble);
        pathsAndCycles[i] = fiber::buildFiberGraphPathsAndCycles(tetMesh, reebSpace, fiberGraphs[i]);


        // Compute barycentric coodrinates of all triangles in the current fiber graph
        for (const auto &[triangleId, componentId] : fiberGraphs[i].componentRoot)
        {
            barycentricCoordinatesPerTriangle[i][triangleId] = computeBarycentricCoordinates(tetMesh, triangleId, fiberPointDouble);
        }
    }

    Timer::stop("Computed fiber graphs and coords       :");










    // Compute all the triangles
    //

    Timer::start();

    std::vector<std::vector<FiberPoint>> allFiberPoints(fiberPoints.size());
    std::cout << "Computing " << fiberPoints.size() << " fibers...\n";

    for (int i = 1 ; i < fiberPoints.size() ; i++)
    {
        // Set up the fiber point
        //
        const Point_2 &fiberPoint = fiberPoints[i];
        const std::array<double, 2> fiberPointDouble = {CGAL::to_double(fiberPoint.x()), CGAL::to_double(fiberPoint.y())};

        // Set up the edge point
        //
        const int edgeId = intersectedSegments[i-1].second;

        const auto [vertexIdA, vertexIdB] = tetMesh.edges.at(edgeId); 
        const std::array<float, 3> vertexCoordinatesA = tetMesh.vertexDomainCoordinates[vertexIdA]; 
        const std::array<float, 3> vertexCoordinatesB = tetMesh.vertexDomainCoordinates[vertexIdB]; 
        const float alpha = intersectedSegmentsAlphas[i-1];
        const std::array<float, 3> edgePointDomain =
        {
            (1.0f - alpha) * vertexCoordinatesA[0] + alpha * vertexCoordinatesB[0],
            (1.0f - alpha) * vertexCoordinatesA[1] + alpha * vertexCoordinatesB[1],
            (1.0f - alpha) * vertexCoordinatesA[2] + alpha * vertexCoordinatesB[2],
        };



        // Extract the paths and cycles of both fibers
        // 
        const auto pathsA = pathsAndCycles[i-1].first;
        const auto cyclesA = pathsAndCycles[i-1].second;


        const auto pathsB = pathsAndCycles[i].first;
        const auto cyclesB = pathsAndCycles[i].second;



        //allFiberPoints[i] =
        //fiber::computeFiberSAT(
        //tetMesh,
        //singularArrangement,
        //reebSpace,
        //fiberPointDouble
        //);

        // Compute correspondence related stuff
        //
        const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(intersectedSegments[i-1].second, typicalOrientation[i-1]);
        const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(intersectedSegments[i-1].second, typicalOrientation[i-1]);

        const std::unordered_set<int> minusTrianglesSet(minusTriangles.begin(), minusTriangles.end());
        const std::unordered_set<int> plusTrianglesSet(plusTriangles.begin(), plusTriangles.end());


        //std::cerr << "\n\n\n\nMinus triangles: " << std::endl;

        //for (const int &triangleId : minusTriangles)
        //{
        //std::cerr << triangleId << std::endl;
        //}

        //std::cerr << "\nPlus triangles: " << std::endl;

        //for (const int &triangleId : plusTriangles)
        //{
        //std::cerr << triangleId << std::endl;
        //}




        std::unordered_set<int> affectedComponentsA;

        for (const int &triangle : minusTriangles)
        {
            affectedComponentsA.insert(fiberGraphs[i-1].componentRoot.at(triangle));
        }

        // Take 
        //
        const std::vector<std::pair<int, int>> correspondence = fiberGraphs[i-1].establishCorrespondence(tetMesh, {intersectedSegments[i-1].second, typicalOrientation[i-1] }, fiberGraphs[i]);




        //printf("\n\n--------------------------------------------------------------------------------------\n");



        //for (const auto &[componentId, path] : pathsA)
        //{
        //std::cout << "\n\nHere's the path fiber from correspondence A with component id " << componentId << "\n";
        //for (const int &triangleId : path)
        //{
        //std::cout << triangleId << " ";
        //}
        //}

        //for (const auto &[componentId, cycle] : cyclesA)
        //{
        //std::cout << "\n\nHere's the cycle fiber from correspondence A with component id " << componentId << "\n";
        //for (const int &triangleId : cycle)
        //{
        //std::cout << triangleId << " ";
        //}

        //}

        //for (const auto &[componentId, path] : pathsB)
        //{
        //std::cout << "\n\nHere's the path fiber from correspondence B with component id " << componentId << "\n";
        //for (const int &triangleId : path)
        //{
        //std::cout << triangleId << " ";
        //}
        //std::cout << std::endl;
        //}

        //for (const auto &[componentId, cycle] : cyclesB)
        //{
        //std::cout << "\n\nHere's the cycle fiber from correspondence B with component id " << componentId << "\n";
        //for (const int &triangleId : cycle)
        //{
        //std::cout << triangleId << " ";
        //}
        //std::cout << std::endl;

        //}



        // For interpolate the non-affected triangles between corresponding components
        //
        for (const auto &[componentA, componentB] : correspondence)
        {
            const int sheetId = reebSpace.correspondenceGraphDS.find(componentA);
            const std::array<float, 3> sheetColour = fiber::fiberColours[sheetId % fiber::fiberColours.size()];


            // If it's a path
            if (pathsA.contains(componentA))
            {
                std::vector<int> path = pathsA.at(componentA);



                for (int j = 0 ; j + 1 < path.size() ; j++)
                {

                    //std::cout << "Adding some triangles\n";


                    const int triangleIdA = path[j];
                    const int triangleIdB = path[j+1];


                    // Make sure these triangle are not active in the trasition
                    //
                    if (
                            true == minusTrianglesSet.contains(triangleIdA) ||
                            true == minusTrianglesSet.contains(triangleIdB)
                            //false == fiberGraphs[i].componentRoot.contains(triangleIdA) ||
                            //false == fiberGraphs[i].componentRoot.contains(triangleIdB)
                       )
                    {
                        continue;
                    }

                    const int triangleIdC = path[j];
                    const int triangleIdD = path[j+1];

                    const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i-1].at(triangleIdA);
                    const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i-1].at(triangleIdB);

                    const std::array<double, 3> &c = barycentricCoordinatesPerTriangle[i].at(triangleIdC);
                    const std::array<double, 3> &d = barycentricCoordinatesPerTriangle[i].at(triangleIdD);

                    std::vector<FiberPoint> componentFiberPoints;

                    //
                    //
                    //      a------c
                    //      |     /|
                    //      |    / |
                    //      |   /  |
                    //      |  /   |
                    //      | /    |
                    //      |/     |
                    //      b------d
                    //      
                    //


                    // Triangle abd
                    componentFiberPoints.push_back(FiberPoint(
                                a[0], 
                                a[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                                sheetColour,
                                sheetId,
                                triangleIdA
                                ));

                    componentFiberPoints.push_back(FiberPoint(
                                b[0], 
                                b[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                                sheetColour,
                                sheetId,
                                triangleIdB
                                ));

                    componentFiberPoints.push_back(FiberPoint(
                                c[0], 
                                c[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdC),
                                sheetColour,
                                sheetId,
                                triangleIdC
                                ));


                    // Triangle dbc
                    componentFiberPoints.push_back(FiberPoint(
                                d[0], 
                                d[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdD),
                                sheetColour,
                                sheetId,
                                triangleIdD
                                ));

                    componentFiberPoints.push_back(FiberPoint(
                                b[0], 
                                b[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                                sheetColour,
                                sheetId,
                                triangleIdB
                                ));

                    componentFiberPoints.push_back(FiberPoint(
                                c[0], 
                                c[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdC),
                                sheetColour,
                                sheetId,
                                triangleIdC
                                ));


                    allFiberPoints[i].push_back(componentFiberPoints[0]);
                    allFiberPoints[i].push_back(componentFiberPoints[1]);
                    allFiberPoints[i].push_back(componentFiberPoints[2]);
                    allFiberPoints[i].push_back(componentFiberPoints[3]);
                    allFiberPoints[i].push_back(componentFiberPoints[4]);
                    allFiberPoints[i].push_back(componentFiberPoints[5]);

                    //faceFibers.push_back(fb);
                    //faceFibers.push_back(fb2);
                }





                //if (pathsA.at(componentA) != pathsB.at(componentB))
                //{
                //printf("\n\n--------------------------------------------------------------------------------------\n");

                //std::cout << "\nHere's the path fiber from correspondence A\n";
                //for (const int &triangleId : pathsA.at(componentA))
                //{
                //std::cout << triangleId << " ";
                //}

                //std::cout << "\n\nHere's the path fiber from correspondence B\n";
                //for (const int &triangleId : pathsB.at(componentB))
                //{
                //std::cout << triangleId << " ";
                //}
                //std::cout << std::endl;

                //return {};

                //throw std::runtime_error("Corresponding paths are not equal.");
                //}

            }
            // If it's a cycle
            else if (cyclesA.contains(componentA))
            {

                std::vector<int> cycle = cyclesA.at(componentA);


                for (int j = 0 ; j  < cycle.size() ; j++)
                {
                    const int triangleIdA = cycle[j];
                    const int triangleIdB = cycle[(j + 1) % (cycle.size())];

                    // If either of the triangle is affected, skip for now
                    if (
                            true == minusTrianglesSet.contains(triangleIdA) ||
                            true == minusTrianglesSet.contains(triangleIdB)
                       )
                    {
                        continue;
                    }

                    // If A and B are not affect, then C and C must not be
                    const int triangleIdC = cycle[j];
                    const int triangleIdD = cycle[(j + 1) % (cycle.size())];

                    const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i-1].at(triangleIdA);
                    const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i-1].at(triangleIdB);

                    const std::array<double, 3> &c = barycentricCoordinatesPerTriangle[i].at(triangleIdC);
                    const std::array<double, 3> &d = barycentricCoordinatesPerTriangle[i].at(triangleIdD);

                    std::vector<FiberPoint> componentFiberPoints;

                    //
                    //
                    //      a------c
                    //      |     /|
                    //      |    / |
                    //      |   /  |
                    //      |  /   |
                    //      | /    |
                    //      |/     |
                    //      b------d
                    //      
                    //

                    //fprintf(stderr, "1) Adding triangle %d A) %d A) %d B)\n", triangleIdA, triangleIdB, triangleIdC);

                    // Triangle abd
                    componentFiberPoints.push_back(FiberPoint(
                                a[0], 
                                a[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                                sheetColour,
                                sheetId,
                                triangleIdA
                                ));

                    componentFiberPoints.push_back(FiberPoint(
                                b[0], 
                                b[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                                sheetColour,
                                sheetId,
                                triangleIdB
                                ));

                    componentFiberPoints.push_back(FiberPoint(
                                c[0], 
                                c[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdC),
                                sheetColour,
                                sheetId,
                                triangleIdC
                                ));


                    //fprintf(stderr, "2) Adding triangle %d B) %d B) %d A)\n", triangleIdD, triangleIdC, triangleIdB);


                    // Triangle dcb
                    componentFiberPoints.push_back(FiberPoint(
                                d[0], 
                                d[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdD),
                                sheetColour,
                                sheetId,
                                triangleIdD
                                ));

                    componentFiberPoints.push_back(FiberPoint(
                                c[0], 
                                c[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdC),
                                sheetColour,
                                sheetId,
                                triangleIdC
                                ));

                    componentFiberPoints.push_back(FiberPoint(
                                b[0], 
                                b[1], 
                                tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                                sheetColour,
                                sheetId,
                                triangleIdB
                                ));

                    allFiberPoints[i].insert(
                            allFiberPoints[i].end(), 
                            std::make_move_iterator(componentFiberPoints.begin()), 
                            std::make_move_iterator(componentFiberPoints.end())
                            );
                }
            }
            else
            {
                throw std::runtime_error("Neither a path nor a cycle.");
            }
        }

        const int intersectedSegmentType = tetMesh.edgeSingularTypes.at(tetMesh.edges.at(intersectedSegments[i-1].second));

        // If we are intersecting a regular segment
        //
        if (intersectedSegmentType == 1 || intersectedSegmentType == -1)
        {
            //std::cerr << "Regular segment " << std::endl;

            // Find the componentIds
            const int componentA = fiberGraphs[i-1].componentRoot.at(minusTriangles[0]);
            const int componentB = fiberGraphs[i].componentRoot.at(plusTriangles[0]);


            const int sheetId = reebSpace.correspondenceGraphDS.find(componentA);
            const std::array<float, 3> sheetColour = fiber::fiberColours[sheetId % fiber::fiberColours.size()];



            // If it's a path
            if (pathsA.contains(componentA))
            {
                //std::cerr << "It's a path " << std::endl;

                const std::vector<FiberPoint> cycleFiberPoints = computeTrianglesBetweenCorrespondingPaths(
                        pathsA.at(componentA), pathsB.at(componentB), 
                        minusTrianglesSet, plusTrianglesSet, 
                        barycentricCoordinatesPerTriangle, 
                        tetMesh, sheetId, sheetColour, i, edgePointDomain);

                allFiberPoints[i].insert(
                        allFiberPoints[i].end(), 
                        std::make_move_iterator(cycleFiberPoints.begin()), 
                        std::make_move_iterator(cycleFiberPoints.end())
                        );

            }


            // If it's a cycle
            else if (cyclesA.contains(componentA))
            {
                const std::vector<FiberPoint> cycleFiberPoints = computeTrianglesBetweenCorrespondingCycles(
                        cyclesA.at(componentA), cyclesB.at(componentB), 
                        minusTrianglesSet, plusTrianglesSet, 
                        barycentricCoordinatesPerTriangle, 
                        tetMesh, sheetId, sheetColour, i, edgePointDomain);

                allFiberPoints[i].insert(
                        allFiberPoints[i].end(), 
                        std::make_move_iterator(cycleFiberPoints.begin()), 
                        std::make_move_iterator(cycleFiberPoints.end())
                        );



            }
            else
            {
                throw std::runtime_error("Neither a path nor a cycle.");
            }
        }





        if (intersectedSegmentType == 0)
        {



            if (minusTrianglesSet.size() == 0)
            {
                const int componentB = fiberGraphs[i].componentRoot.at(plusTriangles[0]);
                const int sheetId = reebSpace.correspondenceGraphDS.find(componentB);
                const std::array<float, 3> sheetColour = fiber::fiberColours[sheetId % fiber::fiberColours.size()];

                if (pathsB.contains(componentB))
                {
                    const auto pathB =  pathsB.at(componentB);

                    for (int k = 0 ; k < pathB.size() - 1 ; k++)
                    {
                        std::vector<FiberPoint> componentFiberPoints;

                        const int triangleIdA = pathB[k];
                        const int triangleIdB = pathB[k + 1];

                        const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i].at(triangleIdA);
                        const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i].at(triangleIdB);

                        //fprintf(stderr, "1 - Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

                        componentFiberPoints.push_back(FiberPoint(
                                    a[0], 
                                    a[1], 
                                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                                    sheetColour,
                                    sheetId,
                                    triangleIdA
                                    ));

                        componentFiberPoints.push_back(FiberPoint(
                                    edgePointDomain,
                                    sheetColour,
                                    sheetId,
                                    -1
                                    ));

                        componentFiberPoints.push_back(FiberPoint(
                                    b[0], 
                                    b[1], 
                                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                                    sheetColour,
                                    sheetId,
                                    triangleIdB
                                    ));

                        allFiberPoints[i].insert(
                                allFiberPoints[i].end(), 
                                std::make_move_iterator(componentFiberPoints.begin()), 
                                std::make_move_iterator(componentFiberPoints.end())
                                );
                    }

                }
                else if (cyclesB.contains(componentB))
                {
                    const auto cycleB =  cyclesB.at(componentB);

                    for (int k = 0 ; k < cycleB.size() ; k++)
                    {
                        std::vector<FiberPoint> componentFiberPoints;

                        const int triangleIdA = cycleB[k];
                        const int triangleIdB = cycleB[(k + 1) % cycleB.size()];

                        const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i].at(triangleIdA);
                        const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i].at(triangleIdB);

                        //fprintf(stderr, "CycB Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

                        componentFiberPoints.push_back(FiberPoint(
                                    a[0], 
                                    a[1], 
                                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                                    sheetColour,
                                    sheetId,
                                    triangleIdA
                                    ));

                        componentFiberPoints.push_back(FiberPoint(
                                    edgePointDomain,
                                    sheetColour,
                                    sheetId,
                                    -1
                                    ));

                        componentFiberPoints.push_back(FiberPoint(
                                    b[0], 
                                    b[1], 
                                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                                    sheetColour,
                                    sheetId,
                                    triangleIdB
                                    ));

                        allFiberPoints[i].insert(
                                allFiberPoints[i].end(), 
                                std::make_move_iterator(componentFiberPoints.begin()), 
                                std::make_move_iterator(componentFiberPoints.end())
                                );
                    }

                }
                else
                {
                    throw std::runtime_error("Neither a path nor a cycle.");
                }

            }

            else
            {
                const int componentA = fiberGraphs[i-1].componentRoot.at(minusTriangles[0]);
                const int sheetId = reebSpace.correspondenceGraphDS.find(componentA);
                const std::array<float, 3> sheetColour = fiber::fiberColours[sheetId % fiber::fiberColours.size()];

                if (pathsA.contains(componentA))
                {
                    const auto pathA =  pathsA.at(componentA);

                    for (int k = 0 ; k < pathA.size() - 1 ; k++)
                    {
                        std::vector<FiberPoint> componentFiberPoints;

                        const int triangleIdA = pathA[k];
                        const int triangleIdB = pathA[k + 1];

                        const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i-1].at(triangleIdA);
                        const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i-1].at(triangleIdB);

                        //fprintf(stderr, "Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

                        componentFiberPoints.push_back(FiberPoint(
                                    a[0], 
                                    a[1], 
                                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                                    sheetColour,
                                    sheetId,
                                    triangleIdA
                                    ));

                        componentFiberPoints.push_back(FiberPoint(
                                    edgePointDomain,
                                    sheetColour,
                                    sheetId,
                                    -1
                                    ));

                        componentFiberPoints.push_back(FiberPoint(
                                    b[0], 
                                    b[1], 
                                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                                    sheetColour,
                                    sheetId,
                                    triangleIdB
                                    ));

                        allFiberPoints[i].insert(
                                allFiberPoints[i].end(), 
                                std::make_move_iterator(componentFiberPoints.begin()), 
                                std::make_move_iterator(componentFiberPoints.end())
                                );
                    }

                }
                else if (cyclesA.contains(componentA))
                {
                    const auto cycleA =  cyclesA.at(componentA);

                    for (int k = 0 ; k < cycleA.size() ; k++)
                    {
                        std::vector<FiberPoint> componentFiberPoints;

                        const int triangleIdA = cycleA[k];
                        const int triangleIdB = cycleA[(k + 1) % cycleA.size()];

                        const std::array<double, 3> &a = barycentricCoordinatesPerTriangle[i-1].at(triangleIdA);
                        const std::array<double, 3> &b = barycentricCoordinatesPerTriangle[i-1].at(triangleIdB);

                        //fprintf(stderr, "CycA Adding triangle %d %d %d\n", triangleIdA, triangleIdB, -1);

                        componentFiberPoints.push_back(FiberPoint(
                                    a[0], 
                                    a[1], 
                                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                                    sheetColour,
                                    sheetId,
                                    triangleIdA
                                    ));

                        componentFiberPoints.push_back(FiberPoint(
                                    edgePointDomain,
                                    sheetColour,
                                    sheetId,
                                    -1
                                    ));

                        componentFiberPoints.push_back(FiberPoint(
                                    b[0], 
                                    b[1], 
                                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                                    sheetColour,
                                    sheetId,
                                    triangleIdB
                                    ));

                        allFiberPoints[i].insert(
                                allFiberPoints[i].end(), 
                                std::make_move_iterator(componentFiberPoints.begin()), 
                                std::make_move_iterator(componentFiberPoints.end())
                                );
                    }

                }
                else
                {
                    throw std::runtime_error("Neither a path nor a cycle.");
                }


            }

        }


        //fprintf(stderr, "\n-----------------------------------------------------------\n");
        //printf("\n\n");
    }

    Timer::stop("Computing fiber surface triangles      :");
    std::cout << std::endl << std::flush;

    std::vector<FiberPoint> result;

    for (const auto& fiberPointsVector : allFiberPoints)
    {
        result.insert(result.end(), fiberPointsVector.begin(), fiberPointsVector.end());
    }


    int totalFiberSize = 0;

    for (int i = 1 ; i < fiberPoints.size() ; i++)
    {
        totalFiberSize += fiberGraphs[i].componentRoot.size();
    }

    std::cerr << "Number of triangles : " << result.size() << std::endl;
    std::cerr << "Number of total fibers : " << totalFiberSize << std::endl;

    return result;
}






// Given a point p within a face F, find a vertex v in the polygon such that the segment pv is entirely within F
// This relies on that fact given a point in a simple, non-intersecting polygon there at least one visible vertex from p.
//
// In practise we use the number of intersected halfEdge, which sould be zero, and only one vertex should be intersected.
//
Point_2 findVisibleVertex(const Face_const_handle activeFace, const Point_2 p, Arrangement &singularArrangement)
{
    auto circ = activeFace->outer_ccb();
    auto start = circ;
    do
    {
        Arrangement_2::X_monotone_curve_2 segmentMonotoneCurve(
                p,
                circ->source()->point());

            std::vector<Arrangement_2::Vertex_handle> vertices;
            std::vector<Arrangement_2::Halfedge_handle> halfEdges;

            CGAL::zone(
                    singularArrangement.arr, 
                    segmentMonotoneCurve, 
                    CGAL::dispatch_or_drop_output<Arrangement_2::Vertex_handle, Arrangement_2::Halfedge_handle>(
                        std::back_inserter(vertices), 
                        std::back_inserter(halfEdges))
                    );

            if (halfEdges.size() == 0 && vertices.size() == 1)
            {
                return circ->source()->point();
            }


        ++circ;
    } while (circ != start);

    throw std::runtime_error("No arrangement vertices are visible from the control point.");
}

FiberGraph fiber::computeFiberGraph(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint)
{
    // Segfautl for ~/Projects/data/reeb-space-test-data/nana/trajectories/State_2/step_00720_d_4.vtu
    //controlPoint = {-0.0566964, 0.107433};

    //controlPoint = {-0.0326639, 0.114234}; // pinch point!

    const Point_2 controlPointEPEC(controlPoint[0], controlPoint[1]);

    //Timer::start();


    // Get the active face
    //
    //Timer::start();
    Face_const_handle activeFace = singularArrangement.getActiveFace(controlPoint);

    if (activeFace->is_unbounded())
    {
        return {};
    }

    //printf("---------------------------------------------------------------------------\n");
    //std::cout << "The active face ID is " << activeFace->data() << std::endl;

    const int activeFaceId = singularArrangement.arrangementFacesIdices[activeFace];
    //Timer::stop("Computed active face in                :");


    // Get a visible point from within the face
    //
    //Timer::start();
    Point_2 closestHalfEdgeVertexPoint = findVisibleVertex(activeFace, controlPointEPEC, singularArrangement);
    //Timer::stop("Found visible point  in                :");

    //std::cout << "The fiber   point is " << controlPointEPEC << std::endl;
    //std::cout << "The visible point is " << closestHalfEdgeVertexPoint << std::endl;

    // Set up the control Segment
    //
    const Segment_2 controlSegment(controlPointEPEC, closestHalfEdgeVertexPoint);

    if (controlSegment.squared_length() == 0.0)
    {
        std::cerr << "The segment has zero lenght!" << std::endl;

    }


    // Compute intersectinos with the AABB tree
    //
    //Timer::start();
    std::vector<TreeAABB::Primitive_id> intersectedSegmentsAABB;
    singularArrangement.tree.all_intersected_primitives(controlSegment, std::back_inserter(intersectedSegmentsAABB));
    //Timer::stop("Computed AABB intersections in         :");


    // Compute the alpha for the intersected segments
    // 
    //Timer::start();

    std::vector<std::pair<K::FT, int>> intersectedSegments;
    intersectedSegments.reserve(intersectedSegmentsAABB.size());

    for (auto id : intersectedSegmentsAABB)
    {
        const Segment_2& s = *id;   // dereference iterator to get the original segment

        // Get the ID of the original segment.
        const int segmentIndex = id - singularArrangement.allSegments.begin();


        const int indexSource = singularArrangement.arrangementPointIndices[s.source()];
        const int indexTarget = singularArrangement.arrangementPointIndices[s.target()];

        const int edgeId = tetMesh.edgeIndices.at({indexSource, indexTarget});

        if (segmentIndex != edgeId)
        {
            throw std::runtime_error("Issue in AABB tree segments indices.");
        }


        // Double check we get the same segment back.
        Point_2 a = singularArrangement.arrangementPoints[tetMesh.edges[segmentIndex][0]];
        Point_2 b = singularArrangement.arrangementPoints[tetMesh.edges[segmentIndex][1]];

        const bool match =
            (s.source() == a && s.target() == b) ||
            (s.source() == b && s.target() == a);

        if (match == false)
        {
            throw std::runtime_error("Issue in AABB tree segments interation.");
        }

        // ---- Collinear check with your controlSegment ----
        if (CGAL::collinear(controlSegment.source(), controlSegment.target(), s.source()) &&
                CGAL::collinear(controlSegment.source(), controlSegment.target(), s.target()))
        {
            // The intersection segment is collinear with the control segment
            std::cerr << "------------------------------------------ Collinear overlap detected for segment " << segmentIndex << std::endl;
            continue;
        }


        //
        //
        //
        //                          s.source()
        //                              |
        //                              |
        //                              |
        // searchSegment.source() ------x---------- searchSegment.target()
        //                              |
        //                              |
        //                              |
        //                              |
        //                          s.target()
        //
        K::FT alpha = CGAL::Intersections::internal::s2s2_alpha(
                controlSegment.target().x(), controlSegment.target().y(),
                controlSegment.source().x(), controlSegment.source().y(),
                s.source().x(), s.source().y(),
                s.target().x(), s.target().y());

        intersectedSegments.emplace_back(alpha, segmentIndex);
    }

    //Timer::stop("Computed Alpha intersections           :");
    //Timer::stop("Computed AABB intersections in         :");


    //Timer::start();
    std::sort(intersectedSegments.begin(), intersectedSegments.end());
    //Timer::stop("Sorting alpha intersections            :");


    for (const auto &[alpha, edgeId] : intersectedSegments)
    {
        const int singularType = tetMesh.edgeSingularTypes.at(tetMesh.edges.at(edgeId));


        if (alpha < 1 && singularType != 1)
        {
            throw std::runtime_error("Control point segment intersects other singular segments.");
        }

        //std::cout << "Intersected segment with ID " << edgeId << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(edgeId)) << " and alpha " << alpha << std::endl;
    }






    //Timer::start();

    // Compute the fiber graph
    //

    // Go up to closestHalfEdgeVertexPoint
    auto currentHalfEdge = activeFace->outer_ccb();

    //std::cout << "Starting half-edge is " << currentHalfEdge->source()->point() << " -> " << currentHalfEdge->target()->point() << std::endl;

    FiberGraph pg = reebSpace.representativeFiberGraphs[activeFace->data()];
    //pg.printByRoot();
    pg.updateComponentsRegular(tetMesh, reebSpace.edgeRegionSegments[currentHalfEdge->data().id]);
    //pg.printByRoot();

    ++currentHalfEdge;
    do
    {
        //std::cout << "Next half-edge is " << currentHalfEdge->source()->point() << " -> " << currentHalfEdge->target()->point() << std::endl;
        pg.updateComponentsRegular(tetMesh, reebSpace.vertexRegionSegments[currentHalfEdge->prev()->data().id]);
        //pg.printByRoot();
        pg.updateComponentsRegular(tetMesh, reebSpace.edgeRegionSegments[currentHalfEdge->data().id]);
        //pg.printByRoot();

        // If we have reached the visible vertex AND the control segment is in its region (not always the case with poinch points)
        if (currentHalfEdge->target()->point() == closestHalfEdgeVertexPoint && reebSpace.ifSegmentInHalfEdgeRegion(currentHalfEdge, controlSegment))
        {
            break;
        }
        else
        {
            ++currentHalfEdge;
        }


    } while (true);


    //std::cout << "Final half-edge is " << currentHalfEdge->source()->point() << " -> " << currentHalfEdge->target()->point() << std::endl;


    // Go to contorlSegment in the region of the vertex
    //
    const auto &vertexRegion = reebSpace.vertexRegionSegments[currentHalfEdge->data().id];

    const int vertexMeshId = singularArrangement.arrangementPointIndices[closestHalfEdgeVertexPoint];

    for (int i = 0 ; i < vertexRegion.size() ; i++)
    {
        const std::pair<int, int> &intersectingSegment = vertexRegion[i];

        const int segmentId = intersectingSegment.first;

        const std::array<int, 2> edge = tetMesh.edges[segmentId];

        Point_2 b;
        if (edge[0] == vertexMeshId)
        {
            b = singularArrangement.arrangementPoints[edge[1]];
        }
        else
        {
            b = singularArrangement.arrangementPoints[edge[0]];
        }

        if (reebSpace.compareRegularSegments(currentHalfEdge, b, controlPointEPEC))
        {
            pg.updateComponentsRegular(tetMesh, {intersectingSegment});
            //std::cout << "It's happening!" << std::endl;
        }
        else
        {
            break;
        }
    }


    //
    // Got along the control segment 


    //std::cout << std::endl;
    //std::cout << std::endl;

    for (int i = intersectedSegments.size() - 1 ; i >= 0 ; i--)
    {
        if (intersectedSegments[i].first == 1.0)
        {
            continue;
        }


        const int segmentId = intersectedSegments[i].second;
        bool typicalOrientation = true;

        //std::cout << "Intersected segment with ID " << segmentId << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentId)) << " and alpha " << intersectedSegments[i].first << std::endl;

        // Change orientation in case we need to
        const std::array<int, 2> edge = tetMesh.edges.at(segmentId);
        const int segmentSourceId = edge[0];
        const int segmentTargetId = edge[1];
        const Point_2 &c = singularArrangement.arrangementPoints[segmentSourceId];
        const Point_2 &d = singularArrangement.arrangementPoints[segmentTargetId];

        //std::cout << c << " - " << d << std::endl;

        const Point_2 &a = closestHalfEdgeVertexPoint;
        const Point_2 &b = controlPointEPEC;
        //
        // If the orientation this way around, reverse it
        //   (regular segment)
        //          d 
        //           \
        //            \
        // a ----------\--------- b (singular segment)
        //              \
        //               \
        //                c
        //
        if (CGAL::orientation(a, b, c) == CGAL::RIGHT_TURN)
        {
            typicalOrientation = false;
        }


        //pg.printByRoot();


        //const std::vector<int> &minusTriangles = tetMesh.getMinusTriangles(segmentId, typicalOrientation);
        //const std::vector<int> &plusTriangles = tetMesh.getPlusTriangles(segmentId, typicalOrientation);

        //std::cout << "\n\n\n\nMinus triangles: " << std::endl;

        //for (const int &triangleId : minusTriangles)
        //{
        //std::cout << triangleId << std::endl;
        //}

        //std::cout << "\nPlus triangles: " << std::endl;

        //for (const int &triangleId : plusTriangles)
        //{
        //std::cout << triangleId << std::endl;
        //}


        pg.updateComponentsRegular(tetMesh, {{segmentId, typicalOrientation}});

    }

    return pg;
}

std::vector<FiberPoint> fiber::computeFiberSAT(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, std::array<double, 2> controlPoint)
{
    //Timer::start();
    FiberGraph fg = fiber::computeFiberGraph(tetMesh, singularArrangement, reebSpace, controlPoint);
    //Timer::stop("Computing fiber graph                  :");

    //Timer::start();
    //const std::vector<FiberPoint> fiber = fiber::processFiberGraph(tetMesh, singularArrangement, reebSpace, controlPoint, fg, {});
    //Timer::stop("Computing fiber                        :");

    //Timer::start();
    const std::vector<FiberPoint> fiber2 = fiber::processFiberGraph2(tetMesh, singularArrangement, reebSpace, controlPoint, fg, {});
    //Timer::stop("Computing fiber graph paths cyc        :");

    return fiber2;
}














std::vector<FiberPoint> fiber::processFiberGraph2(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, FiberGraph &pg, const std::set<int> activeSheets)
{
    // We first need to sort the triangles in the fiber graph
    std::vector<FiberPoint> faceFibers;


    std::unordered_map<int, std::array<double, 3>> barycentricCoordinatesPerTriangle;

    //std::cout << std::endl;
    //std::cout << std::endl;

    //std::cout << "The fiber has " << pg.componentRoot.size() << std::endl;

    //Timer::start();
    for (const auto &[triangleId, componentId] : pg.componentRoot)
    {
        barycentricCoordinatesPerTriangle[triangleId] = computeBarycentricCoordinates(tetMesh, triangleId, fiberPoint);
    }
    //Timer::stop("Computing barycentric coordinates      :");

    const auto [paths, cycles] = fiber::buildFiberGraphPathsAndCycles(tetMesh, reebSpace, pg);

    for (const auto &[componentId, path] : paths)
    {
        const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);
        const std::array<float, 3> sheetColour = fiber::fiberColours[sheetId % fiber::fiberColours.size()];

        for (int i = 0 ; i < path.size() - 1 ; i++)
        {
            const int triangleIdA = path[i];
            const int triangleIdB = path[i+1];

            const std::array<double, 3> &triangleABarycentricCoordinates = barycentricCoordinatesPerTriangle.at(triangleIdA);
            const std::array<double, 3> &triangleBBarycentricCoordinates = barycentricCoordinatesPerTriangle.at(triangleIdB);

            FiberPoint fb(
                    triangleABarycentricCoordinates[0], 
                    triangleABarycentricCoordinates[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                    sheetColour,
                    sheetId,
                    triangleIdA
                    );

            FiberPoint fb2(
                    triangleBBarycentricCoordinates[0], 
                    triangleBBarycentricCoordinates[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                    sheetColour,
                    sheetId,
                    triangleIdB
                    );

            faceFibers.push_back(fb);
            faceFibers.push_back(fb2);
        }
    }

    for (const auto &[componentId, cycle] : cycles)
    {
        const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);
        const std::array<float, 3> sheetColour = fiber::fiberColours[sheetId % fiber::fiberColours.size()];

        for (int i = 0 ; i < cycle.size() ; i++)
        {
            const int triangleIdA = cycle[i];
            const int triangleIdB = cycle[(i + 1) % (cycle.size())];

            const std::array<double, 3> &triangleABarycentricCoordinates = barycentricCoordinatesPerTriangle.at(triangleIdA);
            const std::array<double, 3> &triangleBBarycentricCoordinates = barycentricCoordinatesPerTriangle.at(triangleIdB);

            FiberPoint fb(
                    triangleABarycentricCoordinates[0], 
                    triangleABarycentricCoordinates[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                    sheetColour,
                    sheetId,
                    triangleIdA
                    );

            FiberPoint fb2(
                    triangleBBarycentricCoordinates[0], 
                    triangleBBarycentricCoordinates[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                    sheetColour,
                    sheetId,
                    triangleIdB
                    );

            faceFibers.push_back(fb);
            faceFibers.push_back(fb2);

        }
    }

    return faceFibers;
}

std::vector<FiberPoint> fiber::computeFiberFromFiberGraph(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint)
{
    // The output of this function
    std::vector<FiberPoint> faceFibers;

    // Determine which face has been clicked on
    Face_const_handle activeFace = arrangement.getActiveFace(fiberPoint);
    const int activeFaceId = arrangement.arrangementFacesIdices[activeFace];

    // Determine which sheets are active (contained the active face)
    std::set<int> activeSheets;
    std::cout << "The active faces are : ";
    for (const int &componentId : reebSpace.correspondenceGraph[activeFaceId])
    {
        const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);
        std::cout << sheetId << " ";

        activeSheets.insert(sheetId);
    }
    std::cout << std::endl;

    // Go through all faces and pull their representative fibers for their sheets
    for (auto face = arrangement.arr.faces_begin(); face != arrangement.arr.faces_end(); ++face) 
    {
        if (face->is_unbounded() || !face->has_outer_ccb())
        {
            continue;
        }

        const int faceId = face->data();

        FiberGraph &pg = reebSpace.representativeFiberGraphs[faceId];

        for (const int &componentId : reebSpace.correspondenceGraph[faceId])
        {
            const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);

            if (activeSheets.contains(sheetId))
            {
                std::vector<FiberPoint> newFaceFibers = processFiberGraph(tetMesh, arrangement, reebSpace, fiberPoint, pg, activeSheets);

                faceFibers.insert(
                        faceFibers.end(),
                        newFaceFibers.begin(),
                        newFaceFibers.end()
                        );
            }
        }

    }


    return faceFibers;
}
















std::vector<FiberPoint> fiber::computeFiber(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace &reebSpace, const std::array<double, 2> &fiberPoint, const int reebSheetIdOnly = -1)
{
    Face_const_handle activeFace = arrangement.getActiveFace(fiberPoint);
    const int activeFaceId = arrangement.arrangementFacesIdices[activeFace];

    //Timer::start();

    // The sizes of these data structures are linear in the size of the fiber, not an issue
    std::queue<int> bfsQueue;
    // This also acts as the visited array
    std::unordered_map<int, int> triangleSheetId;
    // Needed to close loops closed fibers, hack because we are not visiting tets, but triangles
    std::unordered_set<std::pair<int, int>, MyHash<std::pair<int, int>>> activeAdjacentTrianglesConnected;
    // Cache barycentric coordintes, they are expensive to compute
    std::unordered_map<int, std::array<double, 3>> triangleBarycentricCoordinates;

    if (reebSheetIdOnly == -1)
    {
        std::cout << "There are " << reebSpace.fiberSeeds[activeFaceId].size() << " fiber components with (sheet IDs, sorted IDs): ";
    }

    //vector<int> sheetIds;
    for (const auto &[fiberComponentId, triangleId] : reebSpace.fiberSeeds[activeFaceId])
    {
        const int sheetId = reebSpace.correspondenceGraph.findElement({activeFaceId, fiberComponentId});

        if (reebSheetIdOnly == -1 || sheetId == reebSheetIdOnly)
        {
            bfsQueue.push(triangleId);
            std::cout << "Adding seed for sheet ID " << sheetId << "\n";
        }

        //const int sheetId = reebSpace.findTriangle({activeFaceId, fiberComponentId});
        //triangleColour[triangleId] = sheetConsequitiveIndices[sheetId] % fiberColours.size();
        triangleSheetId[triangleId] = sheetId;

        //sheetIds.push_back(sheetId);

        if (reebSheetIdOnly == -1)
        {
            printf("(%d, %d) ", sheetId, reebSpace.sheetConsequitiveIndices[sheetId]);
        }
    }
    //std::cout << std::endl;

    CartesianPoint P(fiberPoint[0], fiberPoint[1]);
    std::vector<FiberPoint> faceFibers;

    while (false == bfsQueue.empty())
    {
        const int currentTriangleId = bfsQueue.front();
        const int currentSheeId = triangleSheetId[currentTriangleId];
        bfsQueue.pop();

        const std::array<float, 3> sheetColour = fiber::fiberColours[reebSpace.sheetConsequitiveIndices[currentSheeId] % fiber::fiberColours.size()];

        const std::set<int> triangleUnpacked = tetMesh.triangles[currentTriangleId];
        const std::vector<int> triangleIndices = std::vector<int>(triangleUnpacked.begin(), triangleUnpacked.end());

        std::array<double, 3> barycentricCoordinatesCurrent;

        if (triangleBarycentricCoordinates.contains(currentTriangleId))
        {
            barycentricCoordinatesCurrent = triangleBarycentricCoordinates[currentTriangleId];
        }
        else
        {
            const CartesianPoint A(tetMesh.vertexCoordinatesF[triangleIndices[0]], tetMesh.vertexCoordinatesG[triangleIndices[0]]);
            const CartesianPoint B(tetMesh.vertexCoordinatesF[triangleIndices[1]], tetMesh.vertexCoordinatesG[triangleIndices[1]]);
            const CartesianPoint C(tetMesh.vertexCoordinatesF[triangleIndices[2]], tetMesh.vertexCoordinatesG[triangleIndices[2]]);
            CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, barycentricCoordinatesCurrent.begin());
            triangleBarycentricCoordinates[currentTriangleId] = barycentricCoordinatesCurrent;
        }

        // Sanity check
        assert(barycentricCoordinatesCurrent[0] > 0 && barycentricCoordinatesCurrent[1] > 0 && barycentricCoordinatesCurrent[2] > 0);

        // Look at the neighbours
        for (const int &neighbourTriagleId : tetMesh.tetIncidentTriangles[currentTriangleId])
        {
            if (neighbourTriagleId == currentTriangleId)
            {
                continue;
            }

            // We can't skip visited neighbours, because there may be a fiber between us (completing a circle)
            //if (triangleColour.contains(neighbourTriagle)) { continue; }

            const std::set<int> triangle2Unpacked = tetMesh.triangles[neighbourTriagleId];
            const std::vector<int> triangle2Indices = std::vector<int>(triangle2Unpacked.begin(), triangle2Unpacked.end());



            // The neighbour is active if we've already seen it
            bool isActive = triangleSheetId.contains(neighbourTriagleId);

            // Or if the image of the triangle contains the fiber points
            // We use a fast test to avoid having to use barycentric coordinates all the time
            if (false == isActive)
            {
                CartesianPoint A(tetMesh.vertexCoordinatesF[triangle2Indices[0]], tetMesh.vertexCoordinatesG[triangle2Indices[0]]);
                CartesianPoint B(tetMesh.vertexCoordinatesF[triangle2Indices[1]], tetMesh.vertexCoordinatesG[triangle2Indices[1]]);
                CartesianPoint C(tetMesh.vertexCoordinatesF[triangle2Indices[2]], tetMesh.vertexCoordinatesG[triangle2Indices[2]]);

                std::vector<CartesianPoint> triangle = {A, B, C};
                const auto result = CGAL::bounded_side_2(triangle.begin(), triangle.end(), P);

                isActive = (result == CGAL::ON_BOUNDED_SIDE);
            }

            // Determine if the triangle is active
            if (isActive)
            {
                // Only add the neighbour if we have not already visited it
                if (false == triangleSheetId.contains(neighbourTriagleId))
                {
                    // BFS things
                    bfsQueue.push(neighbourTriagleId);
                    triangleSheetId[neighbourTriagleId] = currentSheeId;
                }

                // Even if we have aleady added a neighbour, maybe there still isn't a fiber between us (for finishing loops)

                // At this point, we know that both us and we neighbour are active, is there already a fiber between us? Then skip
                if (activeAdjacentTrianglesConnected.contains({currentTriangleId, neighbourTriagleId}))
                {
                    continue;
                }
                else
                {
                    activeAdjacentTrianglesConnected.insert({currentTriangleId, neighbourTriagleId});
                    activeAdjacentTrianglesConnected.insert({neighbourTriagleId, currentTriangleId});
                }




                // Compute barycentric coordinates for drawing
                std::array<double, 3> barycentricCoordinatesNeighbour;
                if (triangleBarycentricCoordinates.contains(neighbourTriagleId))
                {
                    barycentricCoordinatesNeighbour = triangleBarycentricCoordinates[neighbourTriagleId];
                }
                else
                {
                    CartesianPoint A(tetMesh.vertexCoordinatesF[triangle2Indices[0]], tetMesh.vertexCoordinatesG[triangle2Indices[0]]);
                    CartesianPoint B(tetMesh.vertexCoordinatesF[triangle2Indices[1]], tetMesh.vertexCoordinatesG[triangle2Indices[1]]);
                    CartesianPoint C(tetMesh.vertexCoordinatesF[triangle2Indices[2]], tetMesh.vertexCoordinatesG[triangle2Indices[2]]);

                    CGAL::Barycentric_coordinates::triangle_coordinates_2(A, B, C, P, barycentricCoordinatesNeighbour.begin());
                    triangleBarycentricCoordinates[neighbourTriagleId] = barycentricCoordinatesNeighbour;
                }

                


                //
                // Add a fiber segment
                //
                FiberPoint fb(
                        barycentricCoordinatesCurrent[0], 
                        barycentricCoordinatesCurrent[1], 
                        {
                            tetMesh.vertexDomainCoordinates[triangleIndices[0]],
                            tetMesh.vertexDomainCoordinates[triangleIndices[1]],
                            tetMesh.vertexDomainCoordinates[triangleIndices[2]],
                        },
                        sheetColour);
                fb.sheetId = currentSheeId;
                fb.triangleId = currentTriangleId;
                faceFibers.push_back(fb);

                FiberPoint fb2(barycentricCoordinatesNeighbour[0], barycentricCoordinatesNeighbour[1], {
                        tetMesh.vertexDomainCoordinates[triangle2Indices[0]],
                        tetMesh.vertexDomainCoordinates[triangle2Indices[1]],
                        tetMesh.vertexDomainCoordinates[triangle2Indices[2]],
                        },
                        sheetColour);
                fb2.sheetId = currentSheeId;
                fb2.triangleId = neighbourTriagleId;
                faceFibers.push_back(fb2);

                //printf("Adding fiber between %d -> %d\n", currentTriangleId, neighbourTriagleId);
            }
        }
    }

    return faceFibers;

    //Timer::stop("Computed fiber in                      :");
}



//
// Brute force fiber computation
//

//void Data::computeTetExitPoints(const GLfloat u, const GLfloat v, const std::vector<float> color)
//{
    //this->faceFibers.clear();
    //this->tetsWithFibers = vector<bool>(this->tetMesh.tetrahedra.size(), false);

    ////
    //// Get the ID of the face we are intersecting
    ////

    //// The query point (u, v)
    //Point_2 query_point(u, v);

    //// Locate the point in the arrangement
    //CGAL::Object result = this->arrangement.pl->locate(query_point);

    //// Try to assign to a face, edge or a vertex
    //Arrangement_2::Face_const_handle face;
    //Arrangement_2::Halfedge_const_handle edge;
    //Arrangement_2::Vertex_const_handle vertex;

    //int currentFaceID = 0;

    //if (CGAL::assign(face, result)) 
    //{
        //currentFaceID = this->arrangement.arrangementFacesIdices[face];
    //} 
    //// If we are on an edge, just grad an adjacent face
    //else if (CGAL::assign(edge, result)) 
    //{
        //face = edge->face();
        //currentFaceID = this->arrangement.arrangementFacesIdices[face];
    //} 
    //// If we are on a vertex grab an indicent edge and get its face
    //else if (CGAL::assign(vertex, result)) 
    //{
        //edge = vertex->incident_halfedges();
        //face = edge->face();
        //currentFaceID = this->arrangement.arrangementFacesIdices[face];
    //} else 
    //{
        //assert(false);
    //}

    //// For every tet, compute the two exit points
    //for(size_t tetId = 0 ; tetId < this->tetMesh.tetrahedra.size(); tetId++)
    //{
        //const auto tet = this->tetMesh.tetrahedra[tetId];

        //// For every triangle in every tet, get the fiber point in it
        //for(int i = 0 ; i < 4 ; i++)
        //{
            //for(int j = i + 1 ; j < 4 ; j++)
            //{
                //for(int k = j + 1 ; k < 4 ; k++)
                //{
                    //float x1 = this->tetMesh.vertexCoordinatesF[tet[i]];
                    //float y1 = this->tetMesh.vertexCoordinatesG[tet[i]];

                    //float x2 = this->tetMesh.vertexCoordinatesF[tet[j]];
                    //float y2 = this->tetMesh.vertexCoordinatesG[tet[j]];

                    //float x3 = this->tetMesh.vertexCoordinatesF[tet[k]];
                    //float y3 = this->tetMesh.vertexCoordinatesG[tet[k]];


                    //const float xmin = std::min({x1, x2, x3});
                    //const float xmax = std::max({x1, x2, x3});
                    //const float ymin = std::min({y1, y2, y3});
                    //const float ymax = std::max({y1, y2, y3});

                    //// This triangle is not relevant, point is outside the bounding box
                    //if (u < xmin || u > xmax || v < ymin || v > ymax) 
                    //{
                        //continue;
                    //}


                    //float det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

                    //float alpha = ((y2 - y3) * (u - x3) + (x3 - x2) * (v - y3)) / det;
                    //float beta = ((y3 - y1) * (u - x3) + (x1 - x3) * (v - y3)) / det;
                    //float gamma = 1 - alpha - beta;

                    //// Are we inside the triangle. We exclude the 0 and 1 because weird things happen there
                    //if (alpha > 0 && beta > 0 && gamma > 0 && alpha < 1 && beta < 1 && gamma < 1)
                    //{
                        ////printf("In triangle %ld, %ld, %ld in tet %ld.\n", tet[i], tet[j], tet[k], tetId);
                        ////printf("In triangle (%f, %f) | (%f, %f) | (%f, %f) comparing with point (%f, %f) and alpha = %f, betta = %f, gamma = %f.\n", x1, y1, x2, y2, x3, y3, x, y, alpha, betta, gamma);

                        ////const std::set<int> triangle({tet[i], tet[j], tet[k]});

                        //const int triangleVertexA = tet[i];
                        //const int triangleVertexB = tet[j];
                        //const int triangleVertexC = tet[k];

                        ////printf("Current triangle (%d, %d, %d)\n", triangleVertexA, triangleVertexB, triangleVertexC);


                        ////for (const auto [key, value] : this->faceDisjointSets[currentFaceID].data)
                        ////{
                        ////cout << "Triangle ";
                        ////for (const auto v : key)
                        ////{
                        ////cout << v << " ";
                        ////}

                        ////cout << " with root " << this->faceDisjointSets[currentFaceID].find(value) << " ( " << this->faceDisjointSets[currentFaceID].findTriangle(key) << ") " << endl;

                        ////}


                        ////for (const auto &[key, value] : this->reebSpace.data)
                        ////{
                        ////printf("Face ID = %d, fiber component root = %d, SheetID = %d\n", key.first, key.second, this->reebSpace.findTriangle(key));

                        ////}


                        //// Which sheets does this fiber belong to?
                        //// 1. Triangle -> Face ComponentID
                        //const int triangleID = this->tetMesh.triangleIndices[std::set<int>({triangleVertexA, triangleVertexB, triangleVertexC})];
                        //const int componentID = this->reebSpace.preimageGraphs[currentFaceID].findElement(triangleID);
                        ////printf("The face ID is %d and the component ID is = %d\n", currentFaceID, componentID);

                        //// 2. Fac ComponentID -> Reeb Space Sheet
                        ////const int pairToHIndex = this->vertexHtoIndex[{currentFaceID, componentID}];
                        ////const int sheetID = this->reebSpace.findTriangle(pairToHIndex);
                        //const int sheetID = this->reebSpace.correspondenceGraph.findElement({currentFaceID, componentID});

                        ////printf("The Sheet ID is = %d\n", sheetID);

                        //const int sheetColourID = this->reebSpace.sheetConsequitiveIndices[sheetID];

                        //// 3. Get the colou of the sheet
                        //const array<float, 3> sheetColour = fiber::fiberColours[sheetColourID];


                        //FiberPoint fb(alpha, beta, {
                                //this->tetMesh.vertexDomainCoordinates[tet[i]],
                                //this->tetMesh.vertexDomainCoordinates[tet[j]],
                                //this->tetMesh.vertexDomainCoordinates[tet[k]],
                                //},
                                //sheetColour);

                        //this->faceFibers.push_back(fb);
                        //this->tetsWithFibers[tetId] = true;
                    //}
                //}
            //}
        //}
    //}
//}




std::vector<std::pair<int, int>> BFSFiberSearch(const int &seedTriangleId, const TetMesh &tetMesh, const FiberGraph &pg, const std::array<double, 2> &fiberPoint, const int sheetId, std::vector<int> &parent)
{
    std::vector<FiberPoint> faceFibers;

    // Set up the BFs
    std::queue<int> bfsQueue;
    bfsQueue.push(seedTriangleId);

    parent[seedTriangleId] = seedTriangleId;

    // Set up the array to hold all the discovered segments
    std::vector<std::pair<int, int>> fiberSegments;

    while (false == bfsQueue.empty())
    {
        const int currentTriangleId = bfsQueue.front();
        bfsQueue.pop();

        for (const int &neighbourTriangleId : tetMesh.tetIncidentTriangles[currentTriangleId])
        {

            // If the neighbour is not active, skip it
            if (false == pg.componentRoot.contains(neighbourTriangleId))
            {
                continue;
            }

            // If the neighbour has NOT been visited, visited and add fiber segment 
            if (parent[neighbourTriangleId] == -1)
            {
                parent[neighbourTriangleId] = currentTriangleId;
                bfsQueue.push(neighbourTriangleId);

                fiberSegments.push_back({
                        std::min(currentTriangleId, neighbourTriangleId),
                        std::max(currentTriangleId, neighbourTriangleId)
                        });

            }

            // If the neighbour has been visited, but it's not our parent, then it completes a cycle, add it
            else if (neighbourTriangleId != parent[currentTriangleId])
            {
                fiberSegments.push_back({
                        std::min(currentTriangleId, neighbourTriangleId),
                        std::max(currentTriangleId, neighbourTriangleId)
                        });

            }


        }
    }


    return fiberSegments;
}


std::vector<FiberPoint> fiber::processFiberGraph(const TetMesh &tetMesh, Arrangement &arrangement, ReebSpace2 &reebSpace, const std::array<double, 2> &fiberPoint, FiberGraph &pg, const std::set<int> activeSheets)
{
    // We first need to sort the triangles in the fiber graph
    std::vector<FiberPoint> faceFibers;


    std::unordered_map<int, std::array<double, 3>> barycentricCoordinatesPerTriangle;

    //std::cout << std::endl;
    //std::cout << std::endl;

    //std::cout << "The fiber has " << pg.componentRoot.size() << std::endl;

    //Timer::start();
    for (const auto &[triangleId, componentId] : pg.componentRoot)
    {
        barycentricCoordinatesPerTriangle[triangleId] = computeBarycentricCoordinates(tetMesh, triangleId, fiberPoint);
    }
    //Timer::stop("Computing barycentric coordinates      :");


    //Timer::start();

    std::vector<int> parent(tetMesh.triangles.size(), -1);
    for (const auto &[triangleId, componentId] : pg.componentRoot)
    {
        if (parent[triangleId] != -1)
        {
            continue;
        }

        const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);
        const std::array<float, 3> sheetColour = fiber::fiberColours[sheetId % fiber::fiberColours.size()];

        //std::cout << "Sheet Id is : " << sheetId << std::endl;

        const std::vector<std::pair<int, int>> fiberSegments = BFSFiberSearch(triangleId, tetMesh, pg, fiberPoint, sheetId, parent);
        //const std::set<std::pair<int, int>> fiberSegments = {};

        // Set up the fiber points for each segment
        for (const auto &[triangleIdA, triangleIdB] : fiberSegments)
        {
            const std::array<double, 3> &triangleABarycentricCoordinates = barycentricCoordinatesPerTriangle.at(triangleIdA);
            const std::array<double, 3> &triangleBBarycentricCoordinates = barycentricCoordinatesPerTriangle.at(triangleIdB);

            FiberPoint fb(
                    triangleABarycentricCoordinates[0], 
                    triangleABarycentricCoordinates[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdA),
                    sheetColour,
                    sheetId,
                    triangleIdA
                    );

            FiberPoint fb2(
                    triangleBBarycentricCoordinates[0], 
                    triangleBBarycentricCoordinates[1], 
                    tetMesh.getTriangleVerticesCoordinates(triangleIdB),
                    sheetColour,
                    sheetId,
                    triangleIdB
                    );

            faceFibers.push_back(fb);
            faceFibers.push_back(fb2);

        }

    }
    //Timer::stop("Computing fiber connectivity           :");

    //std::cout << "There are this many face fibers " << faceFibers.size() << std::endl;

    return faceFibers;

}

std::vector<int> fiber::extractPath(const int &start, const std::unordered_map<int, std::vector<int>> &fgAdj, std::vector<bool> &visited)
{
    std::vector<int> path;

    path.push_back(start);
    visited[start] = true;

    int current = start;

    do
    {
        for (const int &neighbour : fgAdj.at(current))
        {
            if (false == visited[neighbour])
            {
                path.push_back(neighbour);
                visited[neighbour] = true;
                current = neighbour;
            }
        }


    } while (fgAdj.at(current).size() == 2);


    return path;
}

std::vector<int> fiber::extractCycle(const int &start, const std::unordered_map<int, std::vector<int>> &fgAdj, std::vector<bool> &visited)
{
    std::vector<int> cycle;

    cycle.push_back(start);
    visited[start] = true;

    int current = start;

    do
    {
        for (const int &neighbour : fgAdj.at(current))
        {
            if (false == visited[neighbour])
            {
                cycle.push_back(neighbour);
                visited[neighbour] = true;
                current = neighbour;
                break;
            }

            // Finish the loop when we find the cycle end
            if (neighbour == start && cycle.size() > 2)
            {
                current = neighbour;
            }
        }


    } while (current != start);


    return cycle;
}


std::pair<std::map<int, std::vector<int>>, std::map<int, std::vector<int>>> fiber::buildFiberGraphPathsAndCycles(const TetMesh &tetMesh, ReebSpace2 &reebSpace, FiberGraph fg)
{
    // Build the adjacency list
    //
    std::unordered_map<int, std::vector<int>> fgAdj;

    //std::cout << "Building the adjacency list ...\n";


    //std::cout << "The fiber graphs is : \n";
    //fg.printByRoot();


    for (const auto &[triangleId, componentId] : fg.componentRoot)
    {
        for (const int &neighbourTriangleId : tetMesh.tetIncidentTriangles[triangleId])
        {
            if (true == fg.componentRoot.contains(neighbourTriangleId))
            {
                fgAdj[triangleId].push_back(neighbourTriangleId);
            }
        }
    }

    std::vector<bool> visited(tetMesh.triangles.size(), false);

    // Search from the endpoints of paths
    //
    std::map<int, std::vector<int>> paths;   

    for (const auto &[triangleId, neighbours] : fgAdj)
    {
        // If this is an endpoint that has not been visited
        if (neighbours.size() == 1 && visited[triangleId] == false)
        {
            const int componentId = fg.componentRoot.at(triangleId);
            //const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);

            if (paths.contains(componentId))
            {
                throw std::runtime_error("There is already a path with this component Id.");
            }

            paths[componentId] = extractPath(triangleId, fgAdj, visited);

            //std::cout << "Found path: \n";

            //for (const int &triangleId : path)
            //{
                //std::cout << triangleId << " ";
            //}
        }
    }



    // Search from cycles
    //
    std::map<int, std::vector<int>> cycles;   

    for (const auto &[triangleId, neighbours] : fgAdj)
    {
        if (visited[triangleId] == false)
        {
            const int componentId = fg.componentRoot.at(triangleId);
            //const int sheetId = reebSpace.correspondenceGraphDS.find(componentId);

            if (cycles.contains(componentId))
            {
                throw std::runtime_error("There is already a cycle with this component Id.");
            }

            cycles[componentId] = extractCycle(triangleId, fgAdj, visited);

            //std::cout << "\n\nFound cycle: \n";

            //for (const int &triangleId : cycle)
            //{
                //std::cout << triangleId << " ";
            //}
        }
    }

    //printf("\n\n");

    return {paths, cycles};
}
