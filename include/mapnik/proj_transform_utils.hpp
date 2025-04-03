#include <proj.h>
#include <iostream>
#include <cmath>
#include <string>

class bbox
{
  public:
    double minX;
    double minY;
    double maxX;
    double maxY;
    bool initialized;

    bbox(double minX, double minY, double maxX, double maxY)
        : minX(minX)
        , minY(minY)
        , maxX(maxX)
        , maxY(maxY)
        , initialized(true)
    {}
    bbox()
        : minX(0)
        , minY(0)
        , maxX(0)
        , maxY(0)
        , initialized(false)
    {}
    ~bbox() {}

    void fix()
    {
        if (minX > maxX)
        {
            std::swap(minX, maxX);
        }
        if (minY > maxY)
        {
            std::swap(minY, maxY);
        }
    }

    void expand(double x, double y)
    {
        if (!initialized)
        {
            minX = x;
            minY = y;
            maxX = x;
            maxY = y;
            initialized = true;
        }
        if (x < minX)
            minX = x;
        if (y < minY)
            minY = y;
        if (x > maxX)
            maxX = x;
        if (y > maxY)
            maxY = y;
    }
    std::string to_string()
    {
        return "[" + std::to_string(minX) + ", " + std::to_string(minY) + ", " + std::to_string(maxX) + ", " +
               std::to_string(maxY) + "]";
    }
};

struct point
{
    double x;
    double y;

    point(double x, double y)
        : x(x)
        , y(y)
    {}
    ~point() {}
};

bool forward(PJ* trans, double* x, double* y)
{
    double* z = 0;
    PJ_COORD coord = proj_trans(trans, PJ_FWD, proj_coord(*x, *y, 0, 0));
    *x = coord.xy.x; // Longitude, something is not clear for me
    *y = coord.xy.y; // Latitude
    return true;
}

point forward(PJ* proj_const, point p)
{

    point copy = point(p.x, p.y);
    if (!forward(proj_const, &copy.x, &copy.y))
    {
        std::cout << "Error: forward" << std::endl;
        return point(NAN, NAN);
    }
    return copy;
}

void add_to_projected(bbox* projected, point* border)
{
    // Check if projected is empty

    if (border == nullptr)
    {
        // Nothing to do
        return;
    }

    // Append the border to projected
    projected->expand(border->x, border->y);
};

bool is_valid_trans(point coord)
{
    // Check if the coordinates are valid (not NaN)
    return !std::isnan(coord.x) && !std::isnan(coord.y) && !std::isinf(coord.x) && !std::isinf(coord.y);
}

point mid_point(point a, point b)
{
    return point((a.x + b.x) / 2, (a.y + b.y) / 2);
}

/**
 * Binary search for the border
 * A(x,y) is the leftmost point
 * B(x,y) is the rightmost point
 *
 * If A is outside the projecting and B is inside, we check for the the boundary change
 * We do a binary search for the Y axis
 */
bool binary_search(PJ* proj, point* border, point ua, point ub)
{
    int max_iter = 20;

    // Sanity check:
    double minX = std::min(ua.x, ub.x);
    double maxX = std::max(ua.x, ub.x);
    double minY = std::min(ua.y, ub.y);
    double maxY = std::max(ua.y, ub.y);

    ua = point(minX, minY);
    ub = point(maxX, maxY);

    point ta = forward(proj, ua);
    point tb = forward(proj, ub);

    // Initialize the border
    if (!is_valid_trans(ta) && !is_valid_trans(tb))
    {
        return false;
    }

    // If all are valid, this method won't work
    if (is_valid_trans(ta) && is_valid_trans(tb))
    {
        // Set tb
        border->x = tb.x;
        border->y = tb.y;
        return true;
    }

    while (max_iter--)
    {
        point uc = mid_point(ua, ub);

        ta = forward(proj, ua);
        tb = forward(proj, ub);
        point tc = forward(proj, uc);

        if (is_valid_trans(tc))
        {
            // Check what was the last valid
            if (is_valid_trans(ta))
            {
                ua = uc;
            }

            if (is_valid_trans(tb))
            {
                ub = uc;
            }
        }
        else
        {
            if (!is_valid_trans(ta))
            {
                ua = uc;
            }

            if (!is_valid_trans(tb))
            {
                ub = uc;
            }
        }
    }

    // Check which one is valid
    if (is_valid_trans(ta))
    {
        border->x = ta.x;
        border->y = ta.y;
        return true;
    }

    if (is_valid_trans(tb))
    {
        border->x = tb.x;
        border->y = tb.y;
        return true;
    }

    return false;
}

bbox find_border(PJ *proj, bbox boundingbox)
{
    int points = 4;

    bbox projected = bbox();

    double startX = boundingbox.minX;
    double startY = boundingbox.minY;
    double endX = boundingbox.maxX;
    double endY = boundingbox.maxY;

    double stepX = (endX - startX) / points;
    double stepY = (endY - startY) / points;

    if (stepX == 0)
    {
        stepX = 1;
    }

    if (stepY == 0)
    {
        stepY = 1;
    }

    point* border = new point(0, 0);

    // Iterate over x
    for (double y = startY; y <= endY; y += stepY)
    {
        for (double x = startX; x <= endX; x += stepX)
        {
            double lastX = x - stepX;
            if (lastX < startX)
            {
                lastX = startX;
            }
            point ua = point(lastX, y);
            point ub = point(x, y);
            if (binary_search(proj, border, ua, ub))
            {
                add_to_projected(&projected, border);
            }
        }
    }

    // Iterate over y
    for (double x = startX; x <= endX; x += stepX)
    {
        for (double y = startY; y <= endY; y += stepY)
        {
            double lastY = y - stepY;
            if (lastY < startY)
            {
                lastY = startY;
            }
            point ua = point(x, lastY);
            point ub = point(x, y);
            if (binary_search(proj, border, ua, ub))
            {
                add_to_projected(&projected, border);
            }
        }
    }

    delete border;

    return projected;
}
