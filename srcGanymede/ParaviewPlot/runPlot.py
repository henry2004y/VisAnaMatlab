#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
testvtk = LegacyVTKReader(FileNames=['/Users/hyzhou/Documents/research/Ganymede/scripts/test.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1928, 1172]

# show data in view
testvtkDisplay = Show(testvtk, renderView1)
# trace defaults for the display properties.
testvtkDisplay.Representation = 'Outline'
testvtkDisplay.ColorArrayName = ['POINTS', '']
testvtkDisplay.OSPRayScaleArray = 'Pe'
testvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
testvtkDisplay.SelectOrientationVectors = 'magnetic_field'
testvtkDisplay.ScaleFactor = 0.4
testvtkDisplay.SelectScaleArray = 'Pe'
testvtkDisplay.GlyphType = 'Arrow'
testvtkDisplay.GlyphTableIndexArray = 'Pe'
testvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
testvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
testvtkDisplay.ScalarOpacityUnitDistance = 0.06493873527970921

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Slice'
slice1 = Slice(Input=testvtk)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [-1.8100000023841858, 0.0, 0.0]


# get color transfer function/color map for 'Pe'
peLUT = GetColorTransferFunction('Pe')

# show data in view
slice1Display = Show(slice1, renderView1)
# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'Pe']
slice1Display.LookupTable = peLUT
slice1Display.OSPRayScaleArray = 'Pe'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'magnetic_field'
slice1Display.ScaleFactor = 0.4
slice1Display.SelectScaleArray = 'Pe'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'Pe'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'
slice1Display.GaussianRadius = 0.2
slice1Display.SetScaleArray = ['POINTS', 'Pe']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'Pe']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'CSV Reader'
traj_G8csv = CSVReader(FileName=['/Users/hyzhou/Documents/research/Ganymede/scripts/Traj_G8.csv'])

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# get layout
layout1 = GetLayout()

# place view in the layout
layout1.AssignView(2, spreadSheetView1)

# show data in view
traj_G8csvDisplay = Show(traj_G8csv, spreadSheetView1)
# trace defaults for the display properties.
traj_G8csvDisplay.FieldAssociation = 'Row Data'

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=traj_G8csv)
tableToPoints1.XColumn = 'yr month day hr min sec X Y Z Bx By Bz'
tableToPoints1.YColumn = 'yr month day hr min sec X Y Z Bx By Bz'
tableToPoints1.ZColumn = 'yr month day hr min sec X Y Z Bx By Bz'

# set active source
SetActiveSource(traj_G8csv)

# destroy tableToPoints1
Delete(tableToPoints1)
del tableToPoints1

# destroy traj_G8csv
Delete(traj_G8csv)
del traj_G8csv

# create a new 'CSV Reader'
test2csv = CSVReader(FileName=['/Users/hyzhou/Documents/research/Ganymede/scripts/test2.csv'])

# show data in view
test2csvDisplay = Show(test2csv, spreadSheetView1)
# trace defaults for the display properties.
test2csvDisplay.FieldAssociation = 'Row Data'

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=test2csv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'x'
tableToPoints1.ZColumn = 'x'

# Properties modified on tableToPoints1
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# hide data in view
Hide(test2csv, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(Input=tableToPoints1)
programmableFilter1.Script = ''
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# Properties modified on programmableFilter1
programmableFilter1.Script = 'pdi = self.GetPolyDataInput()\npdo =  self.GetPolyDataOutput()\nnumPoints = pdi.GetNumberOfPoints()\npdo.Allocate()\nfor i in range(0, numPoints-1):\n    points = [i, i+1]\n    # VTK_LINE is 3                                                                                                                                   \n    pdo.InsertNextCell(3, 2, points)'
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# show data in view
programmableFilter1Display = Show(programmableFilter1, spreadSheetView1)

# hide data in view
Hide(tableToPoints1, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(programmableFilter1)

# show data in view
programmableFilter1Display_1 = Show(programmableFilter1, renderView1)
# trace defaults for the display properties.
programmableFilter1Display_1.Representation = 'Surface'
programmableFilter1Display_1.ColorArrayName = [None, '']
programmableFilter1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
programmableFilter1Display_1.SelectOrientationVectors = 'None'
programmableFilter1Display_1.ScaleFactor = 0.8832000000000001
programmableFilter1Display_1.SelectScaleArray = 'None'
programmableFilter1Display_1.GlyphType = 'Arrow'
programmableFilter1Display_1.GlyphTableIndexArray = 'None'
programmableFilter1Display_1.DataAxesGrid = 'GridAxesRepresentation'
programmableFilter1Display_1.PolarAxes = 'PolarAxesRepresentation'
programmableFilter1Display_1.GaussianRadius = 0.44160000000000005
programmableFilter1Display_1.SetScaleArray = [None, '']
programmableFilter1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
programmableFilter1Display_1.OpacityArray = [None, '']
programmableFilter1Display_1.OpacityTransferFunction = 'PiecewiseFunction'

# set active source
SetActiveSource(slice1)

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
spreadSheetView1.Update()

# Rescale transfer function
peLUT.RescaleTransferFunction(0.0, 2.43000006676)

# get opacity transfer function/opacity map for 'Pe'
pePWF = GetOpacityTransferFunction('Pe')

# Rescale transfer function
pePWF.RescaleTransferFunction(0.0, 2.43000006676)

# set active source
SetActiveSource(programmableFilter1)

# Properties modified on programmableFilter1Display_1
programmableFilter1Display_1.LineWidth = 4.0

# set active view
SetActiveView(spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
layout1.Collapse(2)

# set active view
SetActiveView(renderView1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-9.224774453475051, -8.192566983981928, 2.1061598712938174]
renderView1.CameraFocalPoint = [-1.8100000023841882, 1.6689388944688312e-15, -1.4654097609970226e-16]
renderView1.CameraViewUp = [0.06400420904926449, 0.1937455266996361, 0.9789617623318384]
renderView1.CameraParallelScale = 2.9113742453882194

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
