import sys
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

def perturb_vtu(input_file, epsilon, output_file):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(input_file)
    reader.Update()
    mesh = reader.GetOutput()

    point_data = mesh.GetPointData()
    for i in range(point_data.GetNumberOfArrays()):
        arr = point_data.GetArray(i)
        np_arr = vtk_to_numpy(arr)
        noise = np.random.uniform(-epsilon, epsilon, np_arr.shape)
        perturbed = np_arr + noise
        vtk_arr = numpy_to_vtk(perturbed)
        vtk_arr.SetName(arr.GetName())
        point_data.GetArray(i).DeepCopy(vtk_arr)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(mesh)
    writer.Write()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python perturb.py <input.vtu> <epsilon> <output.vtu>")
        sys.exit(1)

    perturb_vtu(sys.argv[1], float(sys.argv[2]), sys.argv[3])
    print("Done.")
