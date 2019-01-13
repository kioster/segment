#ifndef POLYGONVERTICE_H
#define POLYGONVERTICE_H

#include <map>
#include <vector>

#include "vtkCell.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

template<class T>
using SPoint = vtkSmartPointer<T>;

class PolygonMap
{
private:
	std::map<int, std::vector<int>> _polygonMap;
	std::vector<std::vector<double>> _feature;
	SPoint<vtkPolyData> _polyData;

	void getCellPointsId(int cellId, std::vector<int>& points);
	bool isBeside(std::vector<int> &point1, std::vector<int> &point2);
	void turnIdToPosition(std::vector<std::vector<double>> &cellPointsPos, std::vector<int> &cellPointsId);

	void initialMap(std::vector<std::vector<int>> &cells);
	void initialFeature(std::vector<std::vector<int>> &cells);
	void computeFeature(std::vector<double>& result, std::vector<std::vector<double>> &cellPointsPos);

public:
	PolygonMap(SPoint<vtkPolyData> initialData = nullptr);

	void getBesideNodes(int id, std::vector<int> &neighbours) { neighbours = this->_polygonMap[id]; }
	void getFeature(int cellId, std::vector<double>& result) { result = this->_feature[cellId]; }
	int getNodesNum() { return this->_polygonMap.size(); }
};

class VectorSizeCompare
{
private:
	std::vector<std::vector<int>> &data;
public:
	VectorSizeCompare(std::vector<std::vector<int>> &data) : data(data){}
	bool operator() (int left, int right) { return data[left].size() < data[right].size();
}
};

int segment(PolygonMap &graph, std::vector<int> &labels);
void groupVertice(std::vector<std::vector<int>> &groups, std::vector<int> &pointsLabel, std::vector<int> &labels, int labelNum, SPoint<vtkPolyData> polygon);
void flagmentsSplice(std::vector<std::vector<int>> &pointsGroups, std::vector<int> &pointsLabels, SPoint<vtkPolyData> polygon);

#endif