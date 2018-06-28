# Plot a line based on vtkTable
pdi = self.GetPolyDataInput()
pdo =  self.GetPolyDataOutput()
numPoints = pdi.GetNumberOfPoints()
pdo.Allocate()
for i in range(0, numPoints-1):
    points = [i, i+1]
    # VTK_LINE is 3
    pdo.InsertNextCell(3, 2, points)
