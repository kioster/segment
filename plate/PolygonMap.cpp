#include "PolygonMap.h"

#include <cmath>
#include <queue>
#include <stack>
#include <set>
#include <iostream>
#include <algorithm>

#include "vtkCell.h"
#include "vtkIdList.h"
#include "vtkPoints.h"

const double DECLINATION = 13.0;
const int MIN_GROUP_SIZE = 100;

void getCellPointsId(SPoint<vtkPolyData> polygon, int cellId, std::vector<int>& points);
void sortGroups(std::vector<std::vector<int>> &groups, std::vector<int> &groupMap);
void getNearPoints(SPoint<vtkPolyData> polygon, int cellId, std::vector<int>& points);
double getDeclination(std::vector<double> &feature1, std::vector<double> &feature2);
int findRoot(std::vector<int> &tree, int node);
void groupAdjust(std::vector<std::vector<int>> &groups, std::vector<int> &pointsLabel);

void PolygonMap::getCellPointsId(int cellId, std::vector<int>& points)
{
	SPoint<vtkCell> result = this->_polyData->GetCell(cellId);
	SPoint<vtkIdList> resultList = result->GetPointIds();

	for (int i = 0; i < resultList->GetNumberOfIds(); i++)
	{
		points[i] = resultList->GetId(i);
	}
}

bool PolygonMap::isBeside(std::vector<int> &cell1, std::vector<int> &cell2)
{
	int intersectionSize = 0;
	for (int i1 = 0; i1 < cell1.size(); i1++)
	{
		for (int i2 = 0; i2 < cell2.size(); i2++)
		{
			intersectionSize += (cell1[i1] == cell2[i2]);
		}
	}

	return intersectionSize == 2;
}

void PolygonMap::turnIdToPosition(std::vector<std::vector<double>> &cellPointsPos, std::vector<int> &cellPointsId)
{
	double *pos = nullptr;
	for (int i = 0; i < cellPointsPos.size(); i++)
	{
		pos = this->_polyData->GetPoint(cellPointsId[i]);
		cellPointsPos[i][0] = pos[0];
		cellPointsPos[i][1] = pos[1];
		cellPointsPos[i][2] = pos[2];
	}
}

void PolygonMap::computeFeature(std::vector<double>& result, std::vector<std::vector<double>> &cellPointsPos)
{
	std::vector<std::vector<double>> vectors(2);
	vectors[0] = { cellPointsPos[1][0] - cellPointsPos[0][0],
		cellPointsPos[1][1] - cellPointsPos[0][1],
		cellPointsPos[1][2] - cellPointsPos[0][2] };
	vectors[1] = { cellPointsPos[2][0] - cellPointsPos[1][0],
		cellPointsPos[2][1] - cellPointsPos[1][1],
		cellPointsPos[2][2] - cellPointsPos[1][2] };

	result = { vectors[0][1] * vectors[1][2] - vectors[0][2] * vectors[1][1],
		vectors[0][2] * vectors[1][0] - vectors[0][0] * vectors[1][2],
		vectors[0][0] * vectors[1][1] - vectors[0][1] * vectors[1][0] };
}

void PolygonMap::initialMap(std::vector<std::vector<int>> &cells)
{
	for (int i0 = 0; i0 < cells.size(); i0++)
	{
		for (int i1 = i0 + 1; i1 < cells.size(); i1++)
		{
			if (this->_polygonMap[i0].size() == 3)
			{
				break;
			}

			if (this->isBeside(cells[i0], cells[i1]))
			{
				this->_polygonMap[i0].push_back(i1);
				this->_polygonMap[i1].push_back(i0);
			}
		}
	}
}

void PolygonMap::initialFeature(std::vector<std::vector<int>> &cells)
{
	std::vector<int> cellPointsId(3);
	std::vector<std::vector<double>> cellPointsPos(3, std::vector<double>(3));
	std::vector<double> feature(3);

	for (int i = 0; i < cells.size(); i++)
	{
		this->getCellPointsId(i, cellPointsId);
		this->turnIdToPosition(cellPointsPos, cellPointsId);
		this->computeFeature(feature, cellPointsPos);
		this->_feature.push_back(feature);
	}
}

PolygonMap::PolygonMap(SPoint<vtkPolyData> initialData)
{
	this->_polyData = initialData;

	if (this->_polyData != nullptr)
	{
		std::vector<int> cellPoints(3);
		std::vector<std::vector<int>> cells;

		for (int i = 0; i < this->_polyData->GetNumberOfCells(); i++)
		{
			this->getCellPointsId(i, cellPoints);
			cells.push_back(cellPoints);
		}

		//std::cout << "initial map" << std::endl;
		this->initialMap(cells);
		//std::cout << "initial feature" << std::endl;
		this->initialFeature(cells);
	}
}

double getDeclination(std::vector<double> &feature1, std::vector<double> &feature2)
{
	double distance1 = std::sqrt(feature1[0] * feature1[0] +
		feature1[1] * feature1[1] +
		feature1[2] * feature1[2]);
	double distance2 = std::sqrt(feature2[0] * feature2[0] +
		feature2[1] * feature2[1] +
		feature2[2] * feature2[2]);

	double dotProduct = feature1[0] * feature2[0] +
		feature1[1] * feature2[1] +
		feature1[2] * feature2[2];

	double result = std::acos(dotProduct / (distance1 * distance2)) * 180 / 3.14;

	if (result > 90.0)
	{
		result = 180.0 - result;
	}

	return result;
}

void getCellPointsId(SPoint<vtkPolyData> polygon, int cellId, std::vector<int>& points)
{
	SPoint<vtkCell> result = polygon->GetCell(cellId);
	SPoint<vtkIdList> resultList = result->GetPointIds();

	for (int i = 0; i < resultList->GetNumberOfIds(); i++)
	{
		points[i] = resultList->GetId(i);
	}
}

void getNearPoints(SPoint<vtkPolyData> polygon, int cellId, std::vector<int>& points)
{
	SPoint<vtkIdList> cellList = SPoint<vtkIdList>::New();
	std::vector<int> cellPointList(3);
	polygon->GetPointCells(cellId, cellList);
	std::set<int> pointsSet;

	for (int i0 = 0; i0 < cellList->GetNumberOfIds(); i0++)
	{
		getCellPointsId(polygon, cellList->GetId(i0), cellPointList);

		for (int i1 = 0; i1 < cellPointList.size(); i1++)
		{
			pointsSet.insert(cellPointList[i1]);
		}
	}

	points.resize(pointsSet.size());
	std::copy(pointsSet.begin(), pointsSet.end(), points.begin());
}

int segment(PolygonMap &graph, std::vector<int> &labels)
{
	int nodesNum = graph.getNodesNum();
	labels.resize(nodesNum, -1);
	std::set<int> edge;

	std::vector<int> neighbours;
	std::stack<int> waitingCells;


	std::vector<double> features[2];
	int firstCell = 0;
	int labelTag = 0;
	labels[firstCell] = labelTag;
	waitingCells.push(firstCell);
	graph.getBesideNodes(firstCell, neighbours);
	double declination = 0.0;

	do
	{
		while (!waitingCells.empty())
		{
			graph.getFeature(firstCell, features[0]);

			for (int i = 0; i < 3; i++)
			{
				graph.getFeature(neighbours[i], features[1]);

				if (labels[neighbours[i]] == -1)
				{
					if (getDeclination(features[0], features[1]) > DECLINATION)
					{
						edge.insert(neighbours[i]);
					}
					else
					{
						firstCell = neighbours[i];
						labels[firstCell] = labels[waitingCells.top()];
						waitingCells.push(firstCell);
						break;
					}
				}

				if (i == 2)
				{
					waitingCells.pop();

					if (waitingCells.empty())
					{
						break;
					}

					firstCell = waitingCells.top();
				}
			}

			graph.getBesideNodes(firstCell, neighbours);
		}

		for (auto it = edge.begin(); it != edge.end();)
		{
			if (labels[*it] != -1)
			{
				it = edge.erase(it);
				continue;
			}

			++it;
		}

		if (!edge.empty())
		{
			firstCell = *edge.begin();
			++labelTag;
			labels[firstCell] = labelTag;
			waitingCells.push(firstCell);
			graph.getBesideNodes(firstCell, neighbours);
		}
		else
		{
			break;
		}
	} while (true);

	return labelTag + 1;
}

void groupAdjust(std::vector<std::vector<int>> &groups, std::vector<int> &pointsLabel)
{
	std::map<int, std::vector<int>> groupMap;

	for (int i = 0; i < pointsLabel.size(); i++)
	{
		groupMap[pointsLabel[i]].push_back(i);
	}

	groups.clear();

	for (auto &item : groupMap)
	{
		groups.push_back(item.second);
	}

	for (int i0 = 0; i0 < groups.size(); i0++)
	{
		for (int i1 = 0; i1 < groups[i0].size(); i1++)
		{
			pointsLabel[groups[i0][i1]] = i0;
		}
	}
}

void groupVertice(std::vector<std::vector<int>> &groups, std::vector<int> &pointsLabel, std::vector<int> &labels, int labelNum, SPoint<vtkPolyData> polygon)
{
	std::vector<std::vector<int>> labelsGroup(labelNum);
	std::vector<int> labelsGroupMap(labelNum, 0);
	VectorSizeCompare compare(labelsGroup);

	for (int i = 0; i < labels.size(); i++)
	{
		labelsGroup[labels[i]].push_back(i);
	}

	for (int i = 0; i < labelsGroupMap.size(); i++)
	{
		labelsGroupMap[i] = i;
	}

	std::sort(labelsGroupMap.begin(), labelsGroupMap.end(), compare);
	pointsLabel.resize(polygon->GetNumberOfPoints());
	std::vector<int> points(3);

	for (int i0 = 0; i0 < labelsGroupMap.size(); i0++)
	{
		for (int i1 = 0; i1 < labelsGroup[labelsGroupMap[i0]].size(); i1++)
		{
			getCellPointsId(polygon, labelsGroup[labelsGroupMap[i0]][i1], points);
			pointsLabel[points[0]] = labelsGroupMap[i0];
			pointsLabel[points[1]] = labelsGroupMap[i0];
			pointsLabel[points[2]] = labelsGroupMap[i0];
		}
	}

	groupAdjust(groups, pointsLabel);
}

int findRoot(std::vector<int> &tree, int node)
{
	int root = node;

	while (tree[root] != -1)
	{
		root = tree[root];

		if (root == node)
		{
			break;
		}
	}

	return root;
}

void flagmentsSplice(std::vector<std::vector<int>> &pointsGroups, std::vector<int> &pointsLabel, SPoint<vtkPolyData> polygon)
{

	VectorSizeCompare compare(pointsGroups);
	std::vector<int> pointsGroupMap(pointsGroups.size());

	for (int i = 0; i < pointsGroupMap.size(); i++)
	{
		pointsGroupMap[i] = i;
	}

	std::sort(pointsGroupMap.begin(), pointsGroupMap.end(), compare);

	std::vector<int> nearPoints;
	int smallGroupId = 0;

	for (int i0 = pointsGroups.size() - 1; i0 >= 0 ; i0--)
	{
		for (int i1 = 0; i1 < pointsGroups[pointsGroupMap[i0]].size(); i1++)
		{
			getNearPoints(polygon, pointsGroups[pointsGroupMap[i0]][i1], nearPoints);

			for (int i2 = 0; i2 < nearPoints.size(); i2++)
			{
				if (pointsLabel[nearPoints[i2]] != pointsGroupMap[i0]
					&& pointsGroups[pointsLabel[nearPoints[i2]]].size() < MIN_GROUP_SIZE)
				{
					smallGroupId = pointsLabel[nearPoints[i2]];

					for (int i3 = 0; i3 < pointsGroups[smallGroupId].size(); i3++)
					{
						pointsGroups[pointsGroupMap[i0]].push_back(pointsGroups[smallGroupId][i3]);
						pointsLabel[pointsGroups[smallGroupId][i3]] = pointsGroupMap[i0];
					}

					pointsGroups[smallGroupId].clear();
				}

			}
		}
	}

	groupAdjust(pointsGroups, pointsLabel);
}

/*void flagmentsSplice(std::vector<std::vector<int>> &pointsGroups, std::vector<int> &pointsLabels, SPoint<vtkPolyData> polygon)
{
	std::vector<int> points;
	std::vector<int> maxNeighbour(pointsGroups.size(), -1);
	int largestGroup, root;

	for (int i0 = 0; i0 < polygon->GetNumberOfPoints(); i0++)
	{
		if (pointsGroups[pointsLabels[i0]].size() > 0 && pointsGroups[pointsLabels[i0]].size() < MIN_GROUP_SIZE)
		{
			getNearPoints(polygon, i0, points);
			largestGroup = findRoot(maxNeighbour, pointsLabels[points[0]]);

			for (int i1 = 1; i1 < points.size(); i1++)
			{
				root = findRoot(maxNeighbour, pointsLabels[points[i1]]);

				if (pointsGroups[largestGroup].size() < pointsGroups[root].size() && pointsLabels[i0] != root)
				{
					largestGroup = pointsLabels[root];
				}
			}

			if (largestGroup != pointsLabels[i0]
				&& (maxNeighbour[pointsLabels[i0]] == -1
				|| pointsGroups[maxNeighbour[pointsLabels[i0]]].size() < pointsGroups[largestGroup].size()))
			{
				maxNeighbour[pointsLabels[i0]] = largestGroup;
			}
		}
	}

	for (int i = 0; i < pointsGroups.size(); i++)
	{
		root = findRoot(maxNeighbour, i);

		if (i != root)
		{
			pointsGroups[root].insert(pointsGroups[root].end(), pointsGroups[i].begin(), pointsGroups[i].end());
			pointsGroups[i].clear();
		}
	}

	for (auto it = pointsGroups.begin(); it != pointsGroups.end();)
	{
		if (it->size() == 0)
		{
			it = pointsGroups.erase(it);
			continue;
		}

		++it;
	}
}*/