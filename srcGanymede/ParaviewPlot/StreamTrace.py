#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
test2csv = CSVReader(FileName=['/Users/hyzhou/Documents/research/Ganymede/scripts/test2.csv'])

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
test2csvDisplay = Show(test2csv, spreadSheetView1)
# trace defaults for the display properties.
test2csvDisplay.FieldAssociation = 'Row Data'

# update the view to ensure updated data information
spreadSheetView1.Update()

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
layout1.Collapse(2)

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1934, 1172]

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=test2csv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)
# trace defaults for the display properties.
tableToPoints1Display.Representation = 'Surface'
tableToPoints1Display.ColorArrayName = [None, '']
tableToPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display.SelectOrientationVectors = 'None'
tableToPoints1Display.ScaleFactor = 0.8806000000000002
tableToPoints1Display.SelectScaleArray = 'None'
tableToPoints1Display.GlyphType = 'Arrow'
tableToPoints1Display.GlyphTableIndexArray = 'None'
tableToPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display.PolarAxes = 'PolarAxesRepresentation'
tableToPoints1Display.GaussianRadius = 0.4403000000000001
tableToPoints1Display.SetScaleArray = [None, '']
tableToPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.OpacityArray = [None, '']
tableToPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Legacy VTK Reader'
testvtk = LegacyVTKReader(FileNames=['/Users/hyzhou/Documents/research/Ganymede/scripts/test.vtk'])

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

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on testvtkDisplay.DataAxesGrid
testvtkDisplay.DataAxesGrid.GridAxesVisibility = 1
testvtkDisplay.DataAxesGrid.XTitle = 'X $[R_G]$'
testvtkDisplay.DataAxesGrid.YTitle = 'Y $[R_G]$'
testvtkDisplay.DataAxesGrid.ZTitle = 'Z $[R_G]$'

# current camera placement for renderView1
renderView1.CameraPosition = [-8.556874081783633, -5.214463231217096, -0.02753961213872662]
renderView1.CameraFocalPoint = [-1.5400963859839703, 0.48187528919379985, -0.02753961213872662]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 4.895565094367408

