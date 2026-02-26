#include <CommandLineParser.h>
#include<unordered_map>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkPLYWriter.h>


#include<ttkFiberSurface.h>
#include<ttkRangePolygon.h>

using namespace ttk;

using namespace std;

string home("/Users/mohit/Desktop/testttk/data/");

#define VTU

#ifdef VTI

#define VTKGRIDTYPE vtkImageData

#define TTKGRIDTYPE ttkImageData

#define DATAREADER vtkXMLImageDataReader

#endif



#ifdef VTU

#define VTKGRIDTYPE vtkUnstructuredGrid

#define TTKGRIDTYPE ttkUnstructuredGrid

#define DATAREADER vtkXMLUnstructuredGridReader

#endif





string dataFilePath("/Users/mohit/Desktop/phd/fiberSurfacepaper/applications/ethane/ethane.vtu");

string fscpFilePath("/Users/mohit/Desktop/phd/fiberSurfacepaper/applications/ethane/");

vtkSmartPointer<DATAREADER> dataReader = vtkSmartPointer<DATAREADER>::New();

vtkSmartPointer<vtkXMLUnstructuredGridReader> cpReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

TTKGRIDTYPE *dataTTKGrid = TTKGRIDTYPE::New();

VTKGRIDTYPE *dataVTKGrid = VTKGRIDTYPE::New();

vtkUnstructuredGrid *cpVTKGrid = vtkUnstructuredGrid::New();;

vtkSmartPointer<ttkTriangulationFilter> triangulationFilter

  = vtkSmartPointer<ttkTriangulationFilter>::New();

Triangulation *triangulation;

long int totalPoints;

string field1Name;

string field2Name;



string dataFiles[] = {"0tq", "0t", "0q"};






void initTriangulation()

{

    triangulationFilter->SetInputData(dataVTKGrid);

    triangulationFilter->Update();

   

    dataTTKGrid->DeepCopy(dynamic_cast<VTKGRIDTYPE *>(triangulationFilter->GetOutput()));

    

    triangulation = dataTTKGrid->getTriangulation();

}

void loadData()

{

    dataReader->SetFileName(dataFilePath.data());

    dataReader->Update();

    dataVTKGrid->DeepCopy(dynamic_cast<VTKGRIDTYPE *>(dataReader->GetOutput()));

    field1Name = dataVTKGrid->GetPointData()->GetArrayName(0);

    field2Name = dataVTKGrid->GetPointData()->GetArrayName(1);

    totalPoints = dataVTKGrid->GetPointData()->GetArray(0)->GetSize();

    initTriangulation();

}



void fiberSurfaceUsingTTK(bool isOctTree, bool isCompleteData, vtkDataObject *data, const char* outputFileName)

{

    cout<<"\n\n\nTime taken using TTK, isOctTree: "<<isOctTree<<", isCompleteData: "<<isCompleteData<<"\n";

    ttkRangePolygon *cp = ttkRangePolygon::New();

    cp->SetInputData(cpReader->GetOutput());

    cp->Update();

    ttkFiberSurface *fs = ttkFiberSurface::New();

    fs->SetInputData(0, data);

    fs->SetInputData(1, cp->GetOutput());

    fs->SetDataUcomponent(field1Name.c_str());

    fs->SetDataVcomponent(field2Name.c_str());

    fs->SetPolygonUcomponent(field1Name.c_str());

    fs->SetPolygonVcomponent(field2Name.c_str());

    fs->SetUseAllCores(false);

    fs->setThreadNumber(1);

    fs->SetRangeOctree(isOctTree);

    fs->setDebugLevel(3);

    //start = clock();

    fs->Update();

    vtkSmartPointer<vtkPLYWriter> sepWriter

          = vtkSmartPointer<vtkPLYWriter>::New();

    sepWriter->SetInputConnection(fs->GetOutputPort());

    sepWriter->SetFileName(outputFileName);

    sepWriter->Write();

}





int main(int argc, char **argv) {

    ttk::CommandLineParser parser;

   

    fiberSurface();

    

  return 0;

}
