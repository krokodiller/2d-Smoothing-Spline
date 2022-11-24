#include "FiniteElement.h"
Point::Point() : x(0), y(0)
{

}

Point::Point(double x, double y) : x(x), y(y)
{

}

bool Point::operator==(const Point& other)
{
    return this->x == other.x && this->y == other.y;
}

FiniteElement::FiniteElement() : nodes(4)
{

}

FiniteElement::FiniteElement(Node leftbottom, Node rightbottom, Node lefttop, Node righttop) : nodes(4)
{
    nodes[0] = leftbottom;
    nodes[1] = rightbottom;
    nodes[2] = lefttop;
    nodes[3] = righttop;
}

Node FiniteElement::leftBottom()
{
    return nodes[Points::leftbottom];
}

Node FiniteElement::rightBottom()
{
    return nodes[Points::rightbottom];
}

Node FiniteElement::leftTop()
{
    return nodes[Points::lefttop];
}

Node FiniteElement::rightTop()
{
    return nodes[Points::righttop];
}
