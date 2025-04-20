/*
Copyright (c) 2023 Grid-based Path Planning Competition and Contributors <https://gppc.search-conference.org/>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef GPPC_VALIDATEPATH_HPP
#define GPPC_VALIDATEPATH_HPP

#include <vector>
#include <cstddef>
#include <cassert>
#include <cstdlib>

namespace inx {

using std::size_t;

struct alignas(size_t)
Point
{
	int x;
	int y;
};

Point operator+(Point lhs, Point rhs) noexcept { return Point{lhs.x + rhs.x, lhs.y + rhs.y}; }
Point operator-(Point lhs, Point rhs) noexcept { return Point{lhs.x - rhs.x, lhs.y - rhs.y}; }

bool operator==(Point lhs, Point rhs) noexcept { return lhs.x == rhs.x && lhs.y == rhs.y; }
bool operator!=(Point lhs, Point rhs) noexcept { return lhs.x != rhs.x || lhs.y != rhs.y; }

class PathValidator
{
public:
	PathValidator(const std::vector<bool>& map, int width, int height)
		: m_map(&map), m_width(static_cast<size_t>(width)), m_height(static_cast<size_t>(height))
	{ }

	bool get(Point u) const noexcept
	{
		return get(u.x, u.y);
	}
	bool get(int x, int y) const noexcept
	{
		assert(static_cast<size_t>(x) < m_width && static_cast<size_t>(y) < m_height);
		return (*m_map)[y * m_width + x];
	}

	bool validPoint(Point u) const noexcept
	{
		return static_cast<size_t>(u.x) < m_width && static_cast<size_t>(u.y) < m_height && get(u);
	}

	bool validEdge(Point u, Point v) const noexcept
	{
		Point uv = v - u;
		Point diff{};
		bool cardinals{};
		if (uv.x == 0) {
			cardinals = true;
			if (uv.y == 0) // u = v
				return true;
			// hori line
			diff = Point{0, uv.y > 0 ? 1 : -1};
		} else if (uv.y == 0) {
			cardinals = true;
			// vert line
			diff = Point{uv.x > 0 ? 1 : -1, 0};
		} else { // non-cardinal line
			cardinals = false;
			if (std::abs(uv.x) != std::abs(uv.y))
				return false; // non-ordinal
			// ordinal line
			diff = Point{uv.x > 0 ? 1 : -1, uv.y > 0 ? 1 : -1};
		}
		// check cells are clear in grid
		if (cardinals)
			return validCardinal(u, v, diff);
		else
			return validOrdinal(u, v, diff);
	}

private:
	bool validCardinal(Point u, Point v, Point diff) const noexcept
	{
		for (Point x = u; x != v; x = x + diff) {
			if (!get(x))
				return false;
		}
		return true;
	}
	bool validOrdinal(Point u, Point v, Point diff) const noexcept
	{
		// check every 2x2 square along u-v
		for (Point x = u; x != v; x = x + diff) {
			if (!get(x) || !get(x.x + diff.x, x.y) || !get(x.x, x.y + diff.y))
				return false;
		}
		if (!get(v)) // seperate from loop to prevent 2x2 squares past v getting checked
			return false;
		return true;
	}

private:
	const std::vector<bool>* m_map;
	size_t m_width;
	size_t m_height;
};

template <typename PathContainer>
int ValidatePath(const std::vector<bool>& map, int width, int height, const PathContainer& thePath)
{
	size_t S = static_cast<size_t>(thePath.size());
	if (S == 0)
		return -1;
	if (S == 1)
		return 0;
	PathValidator validator(map, width, height);
	// check each point in path
	for (size_t i = 0; i < S; ++i) {
		Point u{static_cast<int>(thePath[i].x), static_cast<int>(thePath[i].y)};
		if (!validator.validPoint(u))
			return static_cast<int>(i);
	}
	// check each segment
	for (size_t i = 0; i < S-1; ++i) {
		Point u{static_cast<int>(thePath[i].x), static_cast<int>(thePath[i].y)};
		Point v{static_cast<int>(thePath[i+1].x), static_cast<int>(thePath[i+1].y)};
		if (!validator.validEdge(u, v))
			return static_cast<int>(i);
	}
	return -1;
}

} // namespace inx

#endif // GPPC_VALIDATEPATH_HPP
