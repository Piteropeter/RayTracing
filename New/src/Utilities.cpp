#include "Utilities.h"

Point& Point::operator+(const Point& p)
{
	x = x + p.x;
	y = y + p.y;
	z = z + p.z;
	return *this;
}

Point& Point::operator*(const Point& p)
{
	x = x * p.x;
	y = y * p.y;
	z = z * p.z;
	return *this;
}

bool Point::operator==(const Point& p) const {
	if (x == p.x && y == p.y && z == p.z) return true;
	return false;
}

Color& Color::operator+(const Color& c)
{
	r = r + c.r;
	g = g + c.g;
	b = b + c.b;
	return *this;
}

Color& Color::operator*(const Color& c)
{
	r = r * c.r;
	g = g * c.g;
	b = b * c.b;
	return *this;
}

Color & Color::operator/(unsigned int x)
{
	r = r / x;
	g = g / x;
	b = b / x;
	return *this;
}

Color& Color::operator/(const Color& c)
{
	r = r / c.r;
	g = g / c.g;
	b = b / c.b;
	return *this;
}
