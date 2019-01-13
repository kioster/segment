#include"vtkRenderWindowInteractor.h"
#include"vtkRenderWindow.h"
#include"vtkSmartPointer.h"
#include"vtkRenderer.h"
#include"vtkInteractorStyleTrackballCamera.h"
#include"vtkCylinderSource.h"
#include"vtkPolyDataMapper.h"
#include"vtkActor.h"
#include"vtkSphereSource.h"
#include"vtkSpherePuzzle.h"
#include"vtkSphereWidget.h"
#include"vtkConeSource.h"
#include"vtkCubeSource.h"
#include"vtkTriangleFilter.h"
#include"vtkMassProperties.h"
#include"vtkAssembly.h"
#include"vtkDataArray.h"
#include"vtkAppendPolyData.h"
#include"vtkUnstructuredGridBunykRayCastFunction.h"
#include"vtkBooleanOperationPolyDataFilter.h"
using namespace std;
int main()
{

	//得到盘子的总体积
	vtkSTLReader *reader = vtkSTLReader::New();
	reader->SetFileName("E:\\Program Files (x86)\\stl\\1.stl");
	reader->Update();
	vtkTriangleFilter *Filter = vtkTriangleFilter::New();
	Filter->SetInputConnection(reader->GetOutputPort());
	vtkMassProperties *pro = vtkMassProperties::New();
	pro->SetInputConnection(Filter->GetOutputPort());
	pro->Update();
	vtkFloatArray *MyPointsArray = vtkFloatArray::New();
	
	double vol;
	vol = pro->GetVolume();
	//cout << vol << endl;
	//盘子中的点
	vtkPolyData *MypolyData = reader->GetOutput();
	double	PointsNumber = MypolyData->GetNumberOfPoints();
	cout << PointsNumber << endl;
	for (int i = 0; i < PointsNumber; i++)
	{
		MyPointsArray->InsertTuple1(i, i);
	}

	MypolyData->GetPointData()->SetScalars(MyPointsArray);

	//double Result[5002];
	//ofstream in;
	//in.open("ResultArray.txt", ios::trunc);
	ifstream fin("ResultArray1.txt");
	double Result[5002];
	for (int i = 0; i < 5002; i++)
	{
		fin >> Result[i];
	}

	vtkLookupTable *MyColor = vtkLookupTable::New();
	MyColor->SetNumberOfColors(PointsNumber);
	for (int i = 0; i < PointsNumber; i++)
	{
		MyColor->SetTableValue(i, 1, 1, 1 + (0.5 - Result[i])*(0.5 - Result[i]));
	}
	MyColor->Build();

	//渲染显示
	vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
	mapper->SetInputData(MypolyData);
	mapper->SetScalarRange(0, 5001);
	mapper->SetLookupTable(MyColor);


	vtkActor *act = vtkActor::New();
	act->SetMapper(mapper);

	vtkRenderer *render = vtkRenderer::New();
	render->AddActor(act);
	render->SetBackground(0.0, 0.0, 0.0);

	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(render);
	renWin->SetSize(300, 300);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	vtkInteractorStyleTrackballCamera *style = vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);

	iren->Initialize();
	iren->Start();

	return 0;

}