#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <algorithm>

#include <Windows.h>
#include <direct.h>

#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSTLReader.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkMassProperties.h"
#include "vtkBooleanOperationPolyDataFilter.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkAssembly.h"
#include "vtkLookupTable.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#include "vtkIdList.h"
#include "vtkCell.h"
#include "vtkPoints.h"

#include "PolygonMap.h"

#include "vtkAutoInit.h" 

VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);

const int X = 0;
const int Y = 1;
const int Z = 2;


using RWindow = vtkRenderWindow;

void setEnvironment(int argc, char **argv);

void startTest();
void base(SPoint<vtkSTLReader> reader, std::vector<std::vector<int>> &groups);
void baseTest();

int main(int argc, char **argv)
{
	//setEnvironment(argc, argv);
	srand(time(nullptr));
	startTest();
	
	return 0;

}

void setEnvironment(int argc, char **argv)
{
	putenv("path=..\\dll\\");
	std::cout << getcwd(new char[_MAX_DIR + 1], _MAX_DIR) << std::endl;
}

void startTest()
{
	SPoint<vtkSTLReader> reader = SPoint<vtkSTLReader>::New();
	reader->SetFileName("stl/1.stl");
	reader->Update();

	//Points
	SPoint<vtkPolyData> mypolyData = reader->GetOutput();
	double	pointsNum = mypolyData->GetNumberOfPoints();
	mypolyData->Print(std::cout);

	std::cout << "graph construct start" << std::endl;
	PolygonMap plateGraph(mypolyData);
	std::cout << "graph construct finish" << std::endl;

	std::vector<int> cellLabels;
	int labelNum;
	std::cout << "start segment" << std::endl;
	labelNum = segment(plateGraph, cellLabels);
	std::cout << "segment finish" << std::endl;

	std::vector<std::vector<int>> groups;
	std::vector<int> labelsData;
	std::cout << "group" << std::endl;
	groupVertice(groups, labelsData, cellLabels, labelNum, mypolyData);
	std::cout << "splice" << std::endl;
	flagmentsSplice(groups, labelsData, mypolyData);
	std::cout << "finish" << std::endl;

	base(reader, groups);
}

void base(SPoint<vtkSTLReader> reader, std::vector<std::vector<int>> &groups)
{
	SPoint<vtkTriangleFilter> Filter = SPoint<vtkTriangleFilter>::New();
	Filter->SetInputConnection(reader->GetOutputPort());
	SPoint<vtkMassProperties> pro = SPoint<vtkMassProperties>::New();
	pro->SetInputConnection(Filter->GetOutputPort());
	pro->Update();

	SPoint<vtkFloatArray> MyPointsArray = SPoint<vtkFloatArray>::New();

	//盘子中的点
	SPoint<vtkPolyData> mypolyData = reader->GetOutput();
	double	pointsNum = mypolyData->GetNumberOfPoints();
	//std::fstream groupsData("segment.txt");

	SPoint<vtkCellArray> polys = mypolyData->GetPolys();
	int polysNum = polys->GetNumberOfCells();

	polys->Print(std::cout);

	SPoint<vtkIdList> temp = SPoint<vtkIdList>::New();

	for (int i = 0; i < pointsNum; i++)
	{
		MyPointsArray->InsertTuple1(i, i);
	}

	mypolyData->GetPointData()->SetScalars(MyPointsArray);

	vtkLookupTable *MyColor = vtkLookupTable::New();
	MyColor->SetNumberOfColors(pointsNum);

	double rgb[3]{};

	for (int i = 0; i < pointsNum; i++)
	{
		MyColor->SetTableValue(i, 0.0, 0.0, 0.0);
	}

	for (int i0 = 0; i0 < groups.size(); i0++)
	{
			rgb[0] = (rand() % 256) / 255.0;
			rgb[1] = (rand() % 256) / 255.0;
			rgb[2] = (rand() % 256) / 255.0;
			//groupsData << i0 << std::endl;

			for (int i1 = 0; i1 < groups[i0].size(); i1++)
			{
				MyColor->SetTableValue(groups[i0][i1], rgb[0], rgb[1], rgb[2]);
				//groupsData << groups[i0][i1] << " ";
			}

			//groupsData << '\n' << std::endl;
	}
	//groupsData.close();

	MyColor->Build();

	//Rendering
	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
	mapper->SetInputData(mypolyData);
	mapper->SetScalarRange(0, pointsNum - 1);
	mapper->SetLookupTable(MyColor);


	vtkActor *act = vtkActor::New();
	act->SetMapper(mapper);

	vtkRenderer *render = vtkRenderer::New();
	render->SetGradientBackground(false);
	render->AddActor(act);
	render->SetBackground(0.0, 0.0, 0.0);

	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(render);
	renWin->SetSize(800, 600);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);

	iren->Initialize();
	iren->Start();
}

void baseTest()
{
	vtkSphereSource *sph = vtkSphereSource::New();
	sph->SetThetaResolution(30);
	sph->SetPhiResolution(30);
	sph->SetCenter(0.25, 0, 0.25);
	sph->SetRadius(0.499);
	sph->Update();

	vtkCubeSource *cube = vtkCubeSource::New();
	cube->SetCenter(0, 0, 0);
	cube->SetXLength(0.5);
	cube->SetYLength(0.5);
	cube->SetZLength(0.5);
	cube->Update();

	vtkTriangleFilter *Tf = vtkTriangleFilter::New();
	Tf->SetInputConnection(cube->GetOutputPort());
	vtkMassProperties *polygonProperties = vtkMassProperties::New();
	polygonProperties->SetInputConnection(Tf->GetOutputPort());
	polygonProperties->Update();

	vtkTriangleFilter *Tf1 = vtkTriangleFilter::New();
	Tf1->SetInputConnection(sph->GetOutputPort());
	vtkMassProperties *polygonProperties1 = vtkMassProperties::New();
	polygonProperties1->SetInputConnection(Tf1->GetOutputPort());
	polygonProperties1->Update();


	vtkBooleanOperationPolyDataFilter *add = vtkBooleanOperationPolyDataFilter::New();
	add->SetInputConnection(0, Tf->GetOutputPort());
	add->SetInputConnection(1, Tf1->GetOutputPort());
	add->SetOperationToIntersection();

	vtkMassProperties *polygonProperties2 = vtkMassProperties::New();
	polygonProperties2->SetInputConnection(add->GetOutputPort());
	polygonProperties2->Update();

	double vol = polygonProperties->GetVolume();
	double vol1 = polygonProperties1->GetVolume();
	double vol2 = polygonProperties2->GetVolume();
	double result = vol2 / vol1;

	std::cout << vol << endl;
	std::cout << vol1 << endl;
	std::cout << vol2 << endl;
	std::cout << result << endl;


	vtkPolyDataMapper *sphMapper = vtkPolyDataMapper::New();
	sphMapper->SetInputConnection(sph->GetOutputPort());

	vtkPolyDataMapper *cubeMapper = vtkPolyDataMapper::New();
	cubeMapper->SetInputConnection(cube->GetOutputPort());

	vtkActor *sphActor = vtkActor::New();
	sphActor->SetMapper(sphMapper);

	vtkActor *cubeActor = vtkActor::New();
	cubeActor->SetMapper(cubeMapper);

	vtkAssembly *asse = vtkAssembly::New();
	asse->AddPart(sphActor);
	asse->AddPart(cubeActor);

	vtkRenderer *ren = vtkRenderer::New();
	ren->AddActor(asse);
	ren->SetBackground(0.1, 0.2, 0.4);

	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(ren);
	renWin->SetSize(300, 300);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	vtkInteractorStyleTrackballCamera *style =
		vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);

	iren->Initialize();
	iren->Start();

	sph->Delete();
	sphMapper->Delete();
	sphActor->Delete();
	ren->Delete();
	renWin->Delete();
	iren->Delete();
	style->Delete();
}