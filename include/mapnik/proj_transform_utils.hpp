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

    double x1;
    double y1;

    point(double x, double y)
        : x(x)
        , y(y)
        , x1(NAN)
        , y1(NAN)
    {}

    point(double x, double y, double x1, double y1)
        : x(x)
        , y(y)
        , x1(x1)
        , y1(y1)
    {}

    bool hasBorder() { return !std::isnan(x1) && !std::isnan(y1); }

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
    if (border->hasBorder())
    {
        projected->expand(border->x1, border->y1);

        // Invalidate border
        border->x1 = NAN;
        border->y1 = NAN;
    }
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

bool is_consecutive_sign_change(point a, point b, point c)
{
    /*
        Para que haya un cambio consecutivo, el punto c debe estar entre a y b, es decir:
        abs(abs(a.x - c.x) - abs(a.x - b.x)) < 1e-6 && abs(abs(a.y - c.y) - abs(a.y - b.y)) < 1e-6
    */

    // Comprobacion rapida de signo:
    if (a.x * b.x < 0 || a.y * b.y < 0 || c.x * b.x < 0 || c.y * b.y < 0)
    {
        return false;
    }

    // Comprobacion de que c esta entre a y b:
    return abs(abs(a.x - c.x) - abs(a.x - b.x)) < 1e-6 && abs(abs(a.y - c.y) - abs(a.y - b.y)) < 1e-6;
}

/**
 * Binary search for the border
 * A(x,y) is the leftmost | uppermost point
 * B(x,y) is the rightmost | lowermost point
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
        // Check if the middle is valid
        point uc = mid_point(ua, ub);
        point tc = forward(proj, uc);
        if (is_valid_trans(tc))
        {
            border->x = tc.x;
            border->y = tc.y;
            return true;
        }
        return false;
    }

    // If all are valid, check whether there is a non consecutive sign change: A(-10) B(-12) C(10) vs A(-10) B(-8) C(10)
    if (is_valid_trans(ta) && is_valid_trans(tb))
    {
        // Check the middle for non consecutive sign change
        point uc = mid_point(ua, ub);
        point tc = forward(proj, uc);
        if (!is_consecutive_sign_change(ta, tb, tc))
        {
            // Find again which one is non consecutive.
            // We can not call us recursively because it will end up in an infinite loop
            while (max_iter--)
            {
                // Check the same as later, but now we are checking which part is non consecutive

                uc = mid_point(ua, ub);

                ta = forward(proj, ua);
                tb = forward(proj, ub);
                tc = forward(proj, uc);

                if (!is_consecutive_sign_change(ta, tb, tc))
                {
                    // Check which part is non consecutive: [AC] or [BC]
                    // We know it by checking the sign of the middle point
                    if ((ta.x < 0 && tc.x > 0) || (ta.y < 0 && tc.y > 0) || (ta.x > 0 && tc.x < 0) ||
                        (ta.y > 0 && tc.y < 0))
                    {
                        ub = uc; // We want to iterate over [BC]
                    }
                    else if ((tb.x < 0 && tc.x > 0) || (tb.y < 0 && tc.y > 0) || (tb.x > 0 && tc.x < 0) ||
                             (tb.y > 0 && tc.y < 0))
                    {
                        ua = uc; // We want to iterate over [AC]
                    }
                    else
                    {
                        break;
                    }
                }
            }

            if (is_valid_trans(tc))
            {
                border->x = ta.x;
                border->y = ta.y;
                border->x1 = tb.x;
                border->y1 = tb.y;
                return true;
            }
        }
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

bbox find_border(PJ* proj, bbox boundingbox)
{
    int points = 9;

    bbox projected = bbox();

    double startX = boundingbox.minX;
    double startY = boundingbox.minY;
    double endX = boundingbox.maxX;
    double endY = boundingbox.maxY;

    double stepX = (endX - startX) / points;
    double stepY = (endY - startY) / points;

    // If proj has stere, directly (-180 to 180) in x
    std::cout << "proj->*: " << proj << std::endl;
    // Dump all of proj for debugging

    // Substract the minimum double distance, to avoid double precision problems

    if (stepX == 0)
    {
        stepX = 1;
    }

    if (stepY == 0)
    {
        stepY = 1;
    }

    stepX = stepX - std::numeric_limits<double>::epsilon();
    stepY = stepY - std::numeric_limits<double>::epsilon();

    point* border = new point(0, 0);

    // Iterate over x
    for (int iy = 0; iy <= points; iy++)
    {
        double y = startY + iy * stepY;
        for (int ix = 0; ix <= points; ix++)
        {
            double x = startX + ix * stepX;
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

    for (int ix = 0; ix <= points; ix++)
    {
        double x = startX + ix * stepX;
        for (int iy = 0; iy <= points; iy++)
        {
            double y = startY + iy * stepY;
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
