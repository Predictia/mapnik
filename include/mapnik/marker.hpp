/*****************************************************************************
 *
 * This file is part of Mapnik (c++ mapping toolkit)
 *
 * Copyright (C) 2014 Artem Pavlenko
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *****************************************************************************/

#ifndef MAPNIK_MARKER_HPP
#define MAPNIK_MARKER_HPP

// mapnik
#include <mapnik/image.hpp>
#include <mapnik/image_util.hpp>
#include <mapnik/svg/svg_storage.hpp>
#include <mapnik/svg/svg_path_attributes.hpp>
#include <mapnik/svg/svg_path_adapter.hpp>
#include <mapnik/util/noncopyable.hpp>
#include <mapnik/util/variant.hpp>

// agg
#include "agg_array.h"

// boost
#include <boost/optional.hpp>

// stl
#include <memory>

namespace mapnik
{

using attr_storage = agg::pod_bvector<mapnik::svg::path_attributes>;
using svg_storage_type = mapnik::svg::svg_storage<mapnik::svg::svg_path_storage,attr_storage>;
using svg_path_ptr = std::shared_ptr<svg_storage_type>;
using image_ptr = std::shared_ptr<image_any>;

struct marker_rgba8
{
public:
    marker_rgba8() 
        : bitmap_data_(4,4,true,true)
    {
        // create default OGC 4x4 black pixel
        bitmap_data_.set(0xff000000);
    }

    marker_rgba8(image_rgba8 const & data)
        : bitmap_data_(data) {}
    
    marker_rgba8(image_rgba8 && data)
        : bitmap_data_(std::move(data)) {}

    marker_rgba8(marker_rgba8 const& rhs)
        : bitmap_data_(rhs.bitmap_data_) {}

    marker_rgba8(marker_rgba8 && rhs) noexcept
        : bitmap_data_(std::move(rhs.bitmap_data_)) {}

    box2d<double> bounding_box() const
    {
        double width = bitmap_data_.width();
        double height = bitmap_data_.height();
        return box2d<double>(0, 0, width, height);
    }

    inline std::size_t width() const
    {
        return bitmap_data_.width();
    }

    inline std::size_t height() const
    {
        return bitmap_data_.height();
    }

    image_rgba8 const& get_data() const
    {
        return bitmap_data_;
    }

private:
    image_rgba8 bitmap_data_;
};

struct marker_svg
{
public:
    marker_svg() { }

    marker_svg(mapnik::svg_path_ptr data)
        : vector_data_(data) {}

    marker_svg(marker_svg const& rhs)
        : vector_data_(rhs.vector_data_) {}

    marker_svg(marker_svg && rhs) noexcept
        : vector_data_(rhs.vector_data_) {}

    box2d<double> bounding_box() const
    {
        return vector_data_->bounding_box();
    }

    inline double width() const
    {
        return vector_data_->bounding_box().width();
    }
    inline double height() const
    {
        return vector_data_->bounding_box().height();
    }

    mapnik::svg_path_ptr get_data() const
    {
        return vector_data_;
    }

private:
    mapnik::svg_path_ptr vector_data_;

};

struct marker_null 
{
public:
    box2d<double> bounding_box() const
    {
        return box2d<double>();
    }
    inline double width() const
    {
        return 0;
    }
    inline double height() const
    {
        return 0;
    }
};

using marker_base = util::variant<marker_svg, 
                                  marker_rgba8,
                                  marker_null>;
namespace detail {

struct get_marker_bbox_visitor
{
    template <typename T>
    box2d<double> operator()(T & data) const
    {
        return data.bounding_box();
    }
};

struct get_marker_width_visitor
{
    template <typename T>
    double operator()(T const& data) const
    {
        return static_cast<double>(data.width());
    }
};

struct get_marker_height_visitor
{
    template <typename T>
    double operator()(T const& data) const
    {
        return static_cast<double>(data.height());
    }
};

} // end detail ns

struct marker : marker_base
{
    marker() = default;

    template <typename T>
    marker(T && data) noexcept
        : marker_base(std::move(data)) {}

    double width() const
    {
        return util::apply_visitor(detail::get_marker_width_visitor(),*this);
    }

    double height() const
    {
        return util::apply_visitor(detail::get_marker_height_visitor(),*this);
    }

    box2d<double> bounding_box() const
    {
        return util::apply_visitor(detail::get_marker_bbox_visitor(),*this);
    }
};

} // end mapnik ns

#endif // MAPNIK_MARKER_HPP
