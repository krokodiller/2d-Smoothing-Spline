#pragma once
#include <vector>


class Point
{
public:
    Point();
    Point(double x, double y);
    bool operator==(const Point &other);

    double x, y;
};

using Node = Point;
class FiniteElement
{
public:
    FiniteElement();
    FiniteElement(Node leftbottom, Node rightbottom, Node lefttop, Point righttop);
    Node leftBottom();
    Node rightBottom();
    Node leftTop();
    Node rightTop();

private:
    std::vector<Node> nodes;
};


enum Points
{
    leftbottom,
    rightbottom,
    lefttop,
    righttop
};