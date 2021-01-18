Mesh.MshFileVersion = 2.2;
///Parametrization:
//
// It is assumed that the XY-data arrays are flat and parallel to the surface at a given height and
// corresponding data do not change vertically across a thin layer. 
//
// .... these are the horizontal coordinates of the lower-left (south-west) end of the data array in the mesh [m]:
//
DataRefX={DataRefX};
DataRefY={DataRefY};
//
// ... this is the height of the grav and magnetic data above ground [m]
// 
DataHeightAboveGround={DataHeightAboveGround};
//
// .... this is the spacing of the data array [m]:
//
DataSpacingX={DataSpacingX};
DataSpacingY={DataSpacingY};
//
// ... number of data points in east-west (X) and north-south (Y) direction:
//
DataNumX={DataNumX};
DataNumY={DataNumY};
//
// Note: the resolution specified here should roughly match the resolution of the actual data as input data are interpolated to the resolution in the mesh
//
// ... this is the "thickness" of the data array = the thickness of the vertical layer. 
//
DataMeshSizeVertical={DataMeshSizeVertical};
//
// ... this is the thickness of region below the data area. We call this the core region which is the focus of the inversion
//
CoreThickness={CoreThickness};
//
// ... there is also an air layer and this is its thickness [m] (no updates for density and magnetization here)
//
AirLayerThickness={AirLayerThickness};
//
// ... there is padding around the core and air layer. For the subsurface there will be updates in padding region but not for the air layer  
//
PaddingX={PaddingX};
PaddingY={PaddingY};
PaddingZ={PaddingZ};
PaddingAir={PaddingAir};
//
// ... these are factors by which the DataMeshSizeVertical is raised in the air layer and in the core. 
//
MeshSizeAirFactor={MeshSizeAirFactor};
MeshSizeCoreFactor={MeshSizeCoreFactor};
//
// ... these are factors by which the core and air layer mesh size are raised for the padding zone. 
//
MeshSizePaddingFactor={MeshSizePaddingFactor};
// 
ProgressionType=1;
DataLX=DataSpacingX*(DataNumX-1);
DataLY=DataSpacingY*(DataNumY-1);
Depth=CoreThickness+PaddingY;


MeshSizeAir={MeshSizeAirFactor}*DataMeshSizeVertical;
MeshSizeCore={MeshSizeCoreFactor}*DataMeshSizeVertical;

MeshSizeBase={MeshSizePaddingFactor}*MeshSizeCore;
MeshSizePaddingAir={MeshSizePaddingFactor}*MeshSizeAir;
MeshSizePaddingSurface={MeshSizePaddingFactor}*MeshSizeCore;

DataMeshGround=DataMeshSizeVertical;

//Core region:
Point(1)={{DataRefX,DataRefY, -CoreThickness, MeshSizeCore}};
Point(2)={{DataRefX,DataRefY, 0., DataMeshGround}};
Point(4)={{DataRefX,DataRefY, AirLayerThickness, MeshSizeAir}};

Point(5)={{DataRefX+DataLX,DataRefY, -CoreThickness, MeshSizeCore}};
Point(6)={{DataRefX+DataLX,DataRefY, 0., DataMeshGround}};
Point(8)={{DataRefX+DataLX,DataRefY, AirLayerThickness, MeshSizeAir}};

Point(9)={{DataRefX,DataRefY+DataLY, -CoreThickness, MeshSizeCore}};
Point(10)={{DataRefX,DataRefY+DataLY, 0., DataMeshGround}};
Point(12)={{DataRefX,DataRefY+DataLY, AirLayerThickness, MeshSizeAir}};

Point(13)={{DataRefX+DataLX,DataRefY+DataLY, -CoreThickness, MeshSizeCore}};
Point(14)={{DataRefX+DataLX,DataRefY+DataLY, 0., DataMeshGround}};
Point(16)={{DataRefX+DataLX,DataRefY+DataLY, AirLayerThickness, MeshSizeAir}};

Line(1) = {{6,14}};
Line(2) = {{14,10}};
Line(3) = {{10,2}};
Line(4) = {{2,6}};

Line Loop (1) = {{1:4}};
Plane Surface(1)={{1}};
out[]=Extrude {{0,0,DataMeshSizeVertical}}{{Surface{{1}};Layers{{1}};}};

out1[]=Extrude {{0,0,DataMeshSizeVertical}}{{Surface{{out[0]}};Layers{{1}};}};
Physical Volume("DataArea")={{out1[1]}};

Line(130) = {{5,13}};
Line(131) = {{13,14}};
Line(132) = {{6,5}};
Line Loop(130)={{130,131,-1,132}};
Line(133) = {{5,1}};
Line(134) = {{1,2}};
Line Loop(133) = {{-132,-4,-134,-133}};
Line(135) = {{1,9}};
Line(136) = {{9,10}};
Line Loop(135) = {{135,136,3,-134}};
Line(138) = {{9,13}};
Line Loop(138) = {{-138,136,-2,-131}};
Line Loop(139) = {{138,-130,133,135}};
Plane Surface(130)={{130}};
Plane Surface(133)={{133}};
Plane Surface(135)={{135}};
Plane Surface(138)={{138}};
Plane Surface(139)={{139}};
Surface Loop(140) = {{1,130,133,135,138,139}};
Volume(3) = {{140}};
Physical Volume("Base")={{3}};

Line(140) = {{27,8}};
Line(141) = {{8,16}};
Line(142) = {{16,28}};
Line Loop(140)={{-140,28,-142,-141}};
Line(143) = {{4,8}};
Line(144) = {{4,36}};
Line Loop(143) = {{-143,144,31,140}};
Line(145) = {{4,12}};
Line(146) = {{12,32}};
Line Loop(145) = {{145,146,30,-144}};
Line(148) = {{12,16}};
Line Loop(148) = {{148,142,29,-146}};
Line Loop(149) = {{141,-148,-145,143}};
Plane Surface(140)={{140}};
Plane Surface(143)={{143}};
Plane Surface(145)={{145}};
Plane Surface(148)={{148}};
Plane Surface(149)={{149}};
Surface Loop(150) = {{48,140,143,145,148,149}};
Volume(4) = {{150}};
Physical Volume("Air")={{4, out[1]}};

//Padding:
Point(201)={{DataRefX-PaddingX,DataRefY-PaddingY, -Depth, MeshSizeBase}};
Point(202)={{DataRefX-PaddingX,DataRefY-PaddingY, 0, MeshSizePaddingSurface}};
Point(203)={{DataRefX-PaddingX,DataRefY-PaddingY, AirLayerThickness+PaddingAir, MeshSizePaddingAir}};

Point(204)={{DataRefX+DataLX+PaddingX,DataRefY-PaddingY, -Depth, MeshSizeBase}};
Point(205)={{DataRefX+DataLX+PaddingX,DataRefY-PaddingY, 0, MeshSizePaddingSurface}};
Point(206)={{DataRefX+DataLX+PaddingX,DataRefY-PaddingY, AirLayerThickness+PaddingAir, MeshSizePaddingAir}};

Point(207)={{DataRefX-PaddingX,DataRefY+DataLY+PaddingY, -Depth, MeshSizeBase}};
Point(208)={{DataRefX-PaddingX,DataRefY+DataLY+PaddingY, 0, MeshSizePaddingSurface}};
Point(209)={{DataRefX-PaddingX,DataRefY+DataLY+PaddingY, AirLayerThickness+PaddingAir, MeshSizePaddingAir}};

Point(210)={{DataRefX+DataLX+PaddingX,DataRefY+DataLY+PaddingY, -Depth, MeshSizeBase}};
Point(211)={{DataRefX+DataLX+PaddingX,DataRefY+DataLY+PaddingY, 0, MeshSizePaddingSurface}};
Point(212)={{DataRefX+DataLX+PaddingX,DataRefY+DataLY+PaddingY, AirLayerThickness+PaddingAir, MeshSizePaddingAir}};

Line(229) = {{201, 204}};
Line(230) = {{204, 205}};
Line(231) = {{204, 210}};
Line(232) = {{205, 211}};
Line(233) = {{211, 210}};
Line(234) = {{210, 207}};
Line(235) = {{211, 208}};
Line(236) = {{207, 208}};
Line(237) = {{208, 202}};
Line(238) = {{207, 201}};
Line(239) = {{201, 202}};
Line(240) = {{202, 205}};
Line(241) = {{203, 202}};
Line(242) = {{203, 206}};
Line(243) = {{206, 212}};
Line(244) = {{212, 209}};
Line(245) = {{203, 209}};
Line(246) = {{209, 208}};
Line(247) = {{212, 211}};
Line(248) = {{205, 206}};
Line Loop(217) = {{237, 240, 232, 235}};
Plane Surface(217) = {{1, 217}};
Line Loop(218) = {{242, 243, 244, -245}};
Plane Surface(218) = {{218}};
Line Loop(219) = {{246, 237, -241, 245}};
Plane Surface(219) = {{219}};
Line Loop(220) = {{246, -235, -247, 244}};
Plane Surface(220) = {{220}};
Line Loop(221) = {{243, 247, -232, 248}};
Plane Surface(221) = {{221}};
Line Loop(222) = {{232, 233, -231, 230}};
Plane Surface(222) = {{222}};
Line Loop(223) = {{230, -240, -239, 229}};
Plane Surface(223) = {{223}};
Line Loop(224) = {{238, 239, -237, -236}};
Plane Surface(224) = {{224}};
Line Loop(225) = {{236, -235, 233, 234}};
Plane Surface(225) = {{225}};
Line Loop(226) = {{229, 231, 234, 238}};
Plane Surface(226) = {{226}};
Line Loop(227) = {{248, -242, 241, 240}};
Plane Surface(227) = {{227}};
Physical Surface("TopFace") = {{218}};
Physical Surface("BottomFace") = {{226}};
Surface Loop(6) = {{225, 224, 226, 223, 222, 130, 138, 135, 133,139,217}};
Volume(6) = {{6}};
Physical Volume("PaddingBase") = {{6}};
Surface Loop(7) = {{218, 221,227,220,219,149,148,140,143,145,21,25,13,17,47,39,35,43,217}};
Volume(7) = {{7}};
Physical Volume("PaddingAir") = {{7}};

Field[1]=Box;




