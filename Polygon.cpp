#include "Leonard.hpp"

void Polygon::close_clockwise(
    Polygon &polygon) noexcept
{
    Point tail = polygon.front();
    Point head = polygon.back();
    int edge; // which edge to start with
    if (head.j == 0)
        edge = 0;
    if (head.i == 0)
        edge = 1;
    if (head.j == 256)
        edge = 2;
    if (head.i == 256)
        edge = 3;
    while (!(head == tail))
    {
        switch (edge++ & 3) // rotate through edges clockwise
        {
        case 0: // bottom
            if (tail.j == 0 && tail.i < head.i) // is the tail on this edge
                head = tail;                    // and clockwise from head
            else                                // no take the whole thing
                head = Point(0, 0);             // to the corner
            break;
        case 1: // left
            if (tail.i == 0 && tail.j > head.j)
                head = tail;
            else
                head = Point(0, 256);
            break;
        case 2: // top
            if (tail.j == 256 && tail.i > head.i)
                head = tail;
            else
                head = Point(256, 256);
            break;
        case 3: // right
            if (tail.i == 256 && tail.j < head.j)
                head = tail;
            else
                head = Point(256, 0);
            break;
        }
        polygon.push_back(head);
    }
}

// Ramer–Douglas–Peucker algorithm

Polygon Polygon::simplify(
    const double tolerance) noexcept
{
    Polygon polygon;
    polygon.reserve(size());

    polygon.push_back(front());
    simplify(tolerance * 256, polygon, 0, size() - 1);
    polygon.push_back(back());

    return polygon;
};

void Polygon::simplify(
    const double tolerance,
    Polygon &simplified,
    const size_t lo,
    const size_t hi) const noexcept
{
    if (lo + 1 == hi)   // empty
        return;

    double x = (*this)[hi].i - (*this)[lo].i;
    double y = (*this)[hi].j - (*this)[lo].j;
    double theta = atan2(y, x); // angle of the line from lo to hi
    double c = cos(theta);
    double s = sin(theta);
    double max = 0;
    size_t keep;
    for (size_t mid = lo + 1; mid < hi; mid++)
    {   // rotate each vector lo to mid around lo then the Y coordinate is the
        // perpendicular distance of mid from the line from lo to hi
        double d = std::abs(c * ((*this)[mid].j - (*this)[lo].j) - s * ((*this)[mid].i - (*this)[lo].i));
        if (d > max)
        {
            keep = mid;
            max = d;
        }
    }
    if (max > tolerance) // significant, so some point we must keep in here
    {
        simplify(tolerance, simplified, lo, keep);
        simplified.push_back((*this)[keep]);
        simplify(tolerance, simplified, keep, hi);
    }
    // but if not, we don't need any of the points between lo and hi
}

// convenience routines
// Polygon Polygon::simplify(
//     const double tolerance) const noexcept
// {
//     Polygon polygon;
//     polygon.reserve((*this).size());

//     polygon.push_back((*this)[0]);
//     simplify(tolerance * Parameters::N, polygon, 0, (*this).size() - 1);
//     polygon.push_back((*this)[(*this).size() - 1]);

//     return polygon;
// }

// Polygon Candidate::in_polygon() const noexcept
// {
//     Polygon polygon;
//     polygon.reserve(separatrix.size() + 4);

//     for (auto &point : separatrix)
//         polygon.push_back(point);

//     close_clockwise(polygon);
//     return polygon;
// }

// Polygon Candidate::in_polygon(
//     double tolerance) const noexcept
// {
//     Polygon polygon;
//     polygon.reserve(separatrix.size() + 4);

//     polygon.push_back(separatrix[0]);
//     simplify(tolerance * Parameters::N, polygon, 0, separatrix.size() - 1);
//     polygon.push_back(separatrix[separatrix.size() - 1]);

//     close_clockwise(polygon);
//     return polygon;
// }

// Polygon Candidate::out_polygon() const noexcept
// {
//     Polygon polygon;
//     polygon.reserve(separatrix.size() + 4);

//     for (auto point = separatrix.rbegin(); point != separatrix.rend(); point++)
//         polygon.push_back(*point);

//     close_clockwise(polygon);
//     return polygon;
// }

// Polygon Candidate::out_polygon(
//     double tolerance) const noexcept
// {
//     Polygon polygon;
//     polygon.reserve(separatrix.size() + 4);

//     polygon.push_back(separatrix[0]);
//     simplify(tolerance * Parameters::N, polygon, 0, separatrix.size() - 1);
//     polygon.push_back(separatrix[separatrix.size() - 1]);
//     std::reverse(polygon.begin(), polygon.end());

//     close_clockwise(polygon);
//     return polygon;
// }
