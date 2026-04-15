#include "./Fiber.h"

#include "./FiberStitching.h"
#include "./FiberLabeling.h"
#include "./ColourTable.h"

#include "./Timer.h"
#include "src/ColourTable.h"


std::tuple<std::vector<int>, std::vector<int>, bool> fiber::stitching::getActiveTrianglesInPath(const std::vector<int> &pathA, const std::vector<int> &pathB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet)
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


std::tuple<std::vector<int>, std::vector<int>, bool> fiber::stitching::getActiveTrianglesInCycle(const std::vector<int> &cycleA, const std::vector<int> &cycleB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet)
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

std::vector<FiberPoint> fiber::stitching::computeTrianglesBetweenCorrespondingPaths(const std::vector<int> &pathA, const std::vector<int> &pathB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet, const std::vector<std::unordered_map<int, std::array<double, 3>>> &barycentricCoordinatesPerTriangle, const TetMesh &tetMesh, const int sheetId, const std::array<float, 3> &sheetColour, const int i, const std::array<float, 3> &edgePointDomain)
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

std::vector<FiberPoint> fiber::stitching::computeTrianglesBetweenCorrespondingCycles(const std::vector<int> &cycleA, const std::vector<int> &cycleB, const std::unordered_set<int> &minusTrianglesSet, const std::unordered_set<int> &plusTrianglesSet, const std::vector<std::unordered_map<int, std::array<double, 3>>> &barycentricCoordinatesPerTriangle, const TetMesh &tetMesh, const int sheetId, const std::array<float, 3> &sheetColour, const int i, const std::array<float, 3> &edgePointDomain)
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



std::array<double, 3> fiber::stitching::computeBarycentricCoordinates(const TetMesh &tetMesh, const int &triangleId, const std::array<double, 2> &fiberPoint)
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



std::vector<FiberPoint> fiber::stitching::computeFiberSurface(TetMesh &tetMesh, Arrangement &singularArrangement, ReebSpace2 &reebSpace, const std::vector<std::array<double, 2>> &controlPoints, int _sheetId)
{

    //const Point_2 startPoint(0.479988 , 0.215557);
    //const Point_2 endPoint(0.259815, 0.116635);

    //const Point_2 startPoint(-0.031168186836691678, 0.14228954860228807);
    //const Point_2 endPoint(0.038491389921685215, -0.15976729170485587);

    const Point_2 startPoint(controlPoints[0][0], controlPoints[0][1]);
    const Point_2 endPoint(controlPoints[1][0], controlPoints[1][1]);


    const Segment_2 controlSegment(startPoint, endPoint);

    std::cout << "Start point : " << std::setprecision(17) << startPoint << " end point " << std::setprecision(17) << endPoint << std::endl;


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


        if (tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentIndex)) == 2)
        {
            //std::cout << "---- Intersected segment with ID " << segmentIndex << " and type " << tetMesh.edgeSingularTypes.at(tetMesh.edges.at(segmentIndex)) << " and alpha " << alpha << std::endl;
        }
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

    std::cout << "Computing " << fiberPoints.size() << " fiber graphs...\n";

    // Do the heavy lifting
    //
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0 ; i < fiberPoints.size() ; i++)
    {
        // Set up the fiber points
        const Point_2 &fiberPoint = fiberPoints[i];
        const std::array<double, 2> fiberPointDouble = { CGAL::to_double(fiberPoint.x()), CGAL::to_double(fiberPoint.y()) };

        // Compute the fiber graph and its paths and cycles
        fiberGraphs[i] = fiber::labeling::computeFiberGraph(tetMesh, singularArrangement, reebSpace, fiberPointDouble);
        pathsAndCycles[i] = fiber::stitching::buildFiberGraphPathsAndCycles(tetMesh, reebSpace, fiberGraphs[i]);


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

    std::set<int> uniqueSheetIds;

    //#pragma omp parallel for schedule(dynamic)
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
            const int sheetSortId = reebSpace.sheetOrder.at(sheetId);
            const std::array<float, 3> sheetColour = colours::getColour(sheetSortId);

            if (_sheetId != -1 && sheetId != _sheetId)
            {
                continue;
            }

            uniqueSheetIds.insert(sheetId);

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
            const int sheetSortId = reebSpace.sheetOrder.at(sheetId);

            if (_sheetId != -1 && sheetId != _sheetId)
            {
                continue;
            }

            uniqueSheetIds.insert(sheetId);


            const std::array<float, 3> sheetColour = colours::getColour(sheetSortId);



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
                const int sheetSortId = reebSpace.sheetOrder.at(sheetId);
                const std::array<float, 3> sheetColour = colours::getColour(sheetSortId);

                if (_sheetId != -1 && sheetId != _sheetId)
                {
                    continue;
                }

                uniqueSheetIds.insert(sheetId);

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
                const int sheetSortId = reebSpace.sheetOrder.at(sheetId);
                const std::array<float, 3> sheetColour = colours::getColour(sheetSortId);

                if (_sheetId != -1 && sheetId != _sheetId)
                {
                    continue;
                }

                uniqueSheetIds.insert(sheetId);

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
    std::cerr << "Number of sheets : " << uniqueSheetIds.size() << std::endl;

    return result;
}

std::vector<int> fiber::stitching::extractPath(const int &start, const std::unordered_map<int, std::vector<int>> &fgAdj, std::vector<bool> &visited)
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

std::vector<int> fiber::stitching::extractCycle(const int &start, const std::unordered_map<int, std::vector<int>> &fgAdj, std::vector<bool> &visited)
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

std::pair<std::map<int, std::vector<int>>, std::map<int, std::vector<int>>> fiber::stitching::buildFiberGraphPathsAndCycles(const TetMesh &tetMesh, ReebSpace2 &reebSpace, FiberGraph fg)
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
