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

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=test2csv)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# hide data in view
Hide(test2csv, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Legacy VTK Reader'
testvtk = LegacyVTKReader(FileNames=['/Users/hyzhou/Documents/research/Ganymede/scripts/test.vtk'])

# show data in view
testvtkDisplay = Show(testvtk, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# create a new 'Stream Tracer With Custom Source'
streamTracerWithCustomSource1 = StreamTracerWithCustomSource(Input=testvtk)

# close an empty frame
layout1.Collapse(2)

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1934, 1172]

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(testvtk)

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

# set active source
SetActiveSource(streamTracerWithCustomSource1)

# get color transfer function/color map for 'Pe'
peLUT = GetColorTransferFunction('Pe')

# show data in view
streamTracerWithCustomSource1Display = Show(streamTracerWithCustomSource1, renderView1)
# trace defaults for the display properties.
streamTracerWithCustomSource1Display.Representation = 'Surface'
streamTracerWithCustomSource1Display.ColorArrayName = ['POINTS', 'Pe']
streamTracerWithCustomSource1Display.LookupTable = peLUT
streamTracerWithCustomSource1Display.OSPRayScaleArray = 'Pe'
streamTracerWithCustomSource1Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracerWithCustomSource1Display.SelectOrientationVectors = 'Normals'
streamTracerWithCustomSource1Display.ScaleFactor = 0.3999990463256836
streamTracerWithCustomSource1Display.SelectScaleArray = 'Pe'
streamTracerWithCustomSource1Display.GlyphType = 'Arrow'
streamTracerWithCustomSource1Display.GlyphTableIndexArray = 'Pe'
streamTracerWithCustomSource1Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracerWithCustomSource1Display.PolarAxes = 'PolarAxesRepresentation'
streamTracerWithCustomSource1Display.GaussianRadius = 0.1999995231628418
streamTracerWithCustomSource1Display.SetScaleArray = ['POINTS', 'Pe']
streamTracerWithCustomSource1Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracerWithCustomSource1Display.OpacityArray = ['POINTS', 'Pe']
streamTracerWithCustomSource1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
streamTracerWithCustomSource1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(streamTracerWithCustomSource1)

# reset view to fit data
renderView1.ResetCamera()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-12.97173297578864, -1.395944261579842, 0.0]
renderView1.CameraFocalPoint = [-1.8100000023841893, 2.7373500376594672e-15, 0.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 2.9113742453882194

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
