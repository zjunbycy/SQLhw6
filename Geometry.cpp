#include "Geometry.h"
#include <cmath>
#include <gl/freeglut.h>

#define NOT_IMPLEMENT -1.0

namespace hw6 {

	/*
	 * Envelope functions
	 */
	bool Envelope::contain(double x, double y) const {
		return x >= minX && x <= maxX && y >= minY && y <= maxY;
	}

	bool Envelope::contain(const Envelope& envelope) const {
		// Task 测试Envelope是否包含关系
		// TODO
		return this->minX <= envelope.minX && this->maxX >= envelope.maxX &&
			this->minY <= envelope.minY && this->maxY >= envelope.maxY;
	}

	bool Envelope::intersect(const Envelope& envelope) const {
		// Task 测试Envelope是否相交
		// TODO
		return !(this->maxX < envelope.minX || this->minX > envelope.maxX ||
			this->maxY < envelope.minY || this->minY > envelope.maxY);
	}

	Envelope Envelope::unionEnvelope(const Envelope& envelope) const {
		// Task 合并两个Envelope生成一个新的Envelope
		// TODO
		Envelope res;
		res.minX = std::min(this->minX, envelope.minX);
		res.maxX = std::max(this->maxX, envelope.maxX);
		res.minY = std::min(this->minY, envelope.minY);
		res.maxY = std::max(this->maxY, envelope.maxY);
		return res;
	}

	void Envelope::draw() const {
		glBegin(GL_LINE_STRIP);

		glVertex2d(minX, minY);
		glVertex2d(minX, maxY);
		glVertex2d(maxX, maxY);
		glVertex2d(maxX, minY);
		glVertex2d(minX, minY);

		glEnd();
	}


	/*
	 * 封装一个辅助函数
	 */
	double pointToSegmentDistance(double px, double py, double x1, double y1, double x2, double y2) {
		double dx = x2 - x1;
		double dy = y2 - y1;
		double len2 = dx * dx + dy * dy;

		if (len2 == 0) { // Segment degenerates to a point
			return sqrt(pow(px - x1, 2) + pow(py - y1, 2));
		}

		// Projection parameter t, clamped to [0, 1]
		double t = ((px - x1) * dx + (py - y1) * dy) / len2;
		t = std::max(0.0, std::min(1.0, t));

		// Projected point coordinates
		double proj_x = x1 + t * dx;
		double proj_y = y1 + t * dy;

		// Calculate distance
		return sqrt(pow(px - proj_x, 2) + pow(py - proj_y, 2));
	}

	/*
	 * Points functions
	 */
	double Point::distance(const Point* point) const {
		return sqrt((x - point->x) * (x - point->x) +
			(y - point->y) * (y - point->y));
	}

	double Point::distance(const LineString* line) const {
		double mindist = line->getPointN(0).distance(this);
		for (size_t i = 0; i < line->numPoints() - 1; ++i) {
			double dist = 0;
			double x1 = line->getPointN(i).getX();
			double y1 = line->getPointN(i).getY();
			double x2 = line->getPointN(i + 1).getX();
			double y2 = line->getPointN(i + 1).getY();
			// Task calculate the distance between Point P(x, y) and Line [P1(x1,
			// y1), P2(x2, y2)] (less than 10 lines)
			// TODO
			double dist_seg = pointToSegmentDistance(this->x, this->y, x1, y1, x2, y2);
			if (dist_seg < mindist)
				mindist = dist_seg;
		}
		return mindist;
	}

	double Point::distance(const Polygon* polygon) const {
		LineString line = polygon->getExteriorRing();
		size_t n = line.numPoints();

		bool inPolygon = false;
		// Task whether Point P(x, y) is within Polygon (less than 15 lines)
		// TODO

		for (size_t i = 0; i < n - 1; ++i) {
			Point p1 = line.getPointN(i);
			Point p2 = line.getPointN(i + 1);

			// 检查边是否与射线（向右）相交
			if ((p1.getY() > y) != (p2.getY() > y)) {
				double x_intersect = ((y - p1.getY()) * (p2.getX() - p1.getX())) / (p2.getY() - p1.getY()) + p1.getX();
				if (x < x_intersect)
					inPolygon = !inPolygon;
			}
		}

		double mindist = 0;
		if (!inPolygon)
			mindist = this->distance(&line);
		return mindist;
	}

	bool Point::intersects(const Envelope& rect) const {
		return (x >= rect.getMinX()) && (x <= rect.getMaxX()) &&
			(y >= rect.getMinY()) && (y <= rect.getMaxY());
	}

	void Point::draw() const {
		glBegin(GL_POINTS);
		glVertex2d(x, y);
		glEnd();
	}

	/*
	 * LineString functions
	 */
	void LineString::constructEnvelope() {
		double minX, minY, maxX, maxY;
		maxX = minX = points[0].getX();
		maxY = minY = points[0].getY();
		for (size_t i = 1; i < points.size(); ++i) {
			maxX = std::max(maxX, points[i].getX());
			maxY = std::max(maxY, points[i].getY());
			minX = std::min(minX, points[i].getX());
			minY = std::min(minY, points[i].getY());
		}
		envelope = Envelope(minX, maxX, minY, maxY);
	}

	double LineString::distance(const Point* point) const {
		// 计算点到线串的距离
		if (this->numPoints() == 0) return 0.0;
		if (this->numPoints() == 1) return this->getPointN(0).distance(point);

		double mindist = this->getPointN(0).distance(point);
		for (size_t i = 0; i < this->numPoints() - 1; ++i) {
			double x1 = this->getPointN(i).getX();
			double y1 = this->getPointN(i).getY();
			double x2 = this->getPointN(i + 1).getX();
			double y2 = this->getPointN(i + 1).getY();
			double dist_seg = pointToSegmentDistance(point->getX(), point->getY(), x1, y1, x2, y2);
			if (dist_seg < mindist)
				mindist = dist_seg;
		}
		return mindist;
	}

	double LineString::distance(const LineString* line) const {
		// TODO
		// 处理空线串边界情况
		if (this->numPoints() < 2 && line->numPoints() < 2) return 0.0;
		if (this->numPoints() == 1 && line->numPoints() == 1) return this->getPointN(0).distance(&line->getPointN(0));
		if (this->numPoints() == 1) return this->getPointN(0).distance(line);
		if (line->numPoints() == 1) return line->getPointN(0).distance(this);

		// 初始化最短距离为无穷大
		double min_dist = INFINITY;

		// 遍历当前线串的每个线段
		for (size_t i = 0; i < this->numPoints() - 1; ++i) {
			double xA1 = this->getPointN(i).getX();
			double yA1 = this->getPointN(i).getY();
			double xA2 = this->getPointN(i + 1).getX();
			double yA2 = this->getPointN(i + 1).getY();

			// 遍历目标线串的每个线段
			for (size_t j = 0; j < line->numPoints() - 1; ++j) {
				double xB1 = line->getPointN(j).getX();
				double yB1 = line->getPointN(j).getY();
				double xB2 = line->getPointN(j + 1).getX();
				double yB2 = line->getPointN(j + 1).getY();

				// 计算两个线段的最短距离 (发生在端点到另一条线段上)
				double d1 = pointToSegmentDistance(xA1, yA1, xB1, yB1, xB2, yB2); // segA起点到segB的距离
				double d2 = pointToSegmentDistance(xA2, yA2, xB1, yB1, xB2, yB2); // segA终点到segB的距离
				double d3 = pointToSegmentDistance(xB1, yB1, xA1, yA1, xA2, yA2); // segB起点到segA的距离
				double d4 = pointToSegmentDistance(xB2, yB2, xA1, yA1, xA2, yA2); // segB终点到segA的距离

				double seg_min = std::min({ d1, d2, d3, d4 });

				// 更新全局最短距离
				if (seg_min < min_dist) {
					min_dist = seg_min;
				}
			}
		}

		return (min_dist == INFINITY) ? 0.0 : min_dist;
	}



	double LineString::distance(const Polygon* polygon) const {
		// TODO
		// 1. 获取线串到多边形外环的距离 (如果相交，距离为0)
		double min_dist = this->distance(&polygon->getExteriorRing());

		if (min_dist == 0.0) {
			return 0.0; // 相交，距离为0
		}

		// 2. 检查线串是否完全在多边形内部 (如果完全在内部且不相交，距离为0)
		// 遍历所有顶点，如果任一点在外部，则线串不完全在内部。
		bool all_inside = true;
		for (size_t i = 0; i < this->numPoints(); ++i) {
			// Point::distance(polygon) == 0.0 表示点在内部或边界上
			if (this->getPointN(i).distance(polygon) > 0.0) {
				all_inside = false;
				break;
			}
		}

		if (all_inside) {
			return 0.0;
		}

		// 3. 如果不相交，且不完全在内部，则最小距离就是线串到外环的最小距离
		return min_dist;
	}

	typedef int OutCode;

	const int INSIDE = 0; // 0000
	const int LEFT = 1;   // 0001
	const int RIGHT = 2;  // 0010
	const int BOTTOM = 4; // 0100
	const int TOP = 8;    // 1000

	// Compute the bit code for a point (x, y) using the clip rectangle
	// bounded diagonally by (xmin, ymin), and (xmax, ymax)
	// ASSUME THAT xmax, xmin, ymax and ymin are global constants.
	OutCode ComputeOutCode(double x, double y, double xmin, double xmax,
		double ymin, double ymax) {
		OutCode code;

		code = INSIDE; // initialised as being inside of [[clip window]]

		if (x < xmin) // to the left of clip window
			code |= LEFT;
		else if (x > xmax) // to the right of clip window
			code |= RIGHT;
		if (y < ymin) // below the clip window
			code |= BOTTOM;
		else if (y > ymax) // above the clip window
			code |= TOP;

		return code;
	}

	// Cohen–Sutherland clipping algorithm clips a line from
	// P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with
	// diagonal from (xmin, ymin) to (xmax, ymax).
	bool intersectTest(double x0, double y0, double x1, double y1, double xmin,
		double xmax, double ymin, double ymax) {
		// compute outcodes for P0, P1, and whatever point lies outside the clip
		// rectangle
		OutCode outcode0 = ComputeOutCode(x0, y0, xmin, xmax, ymin, ymax);
		OutCode outcode1 = ComputeOutCode(x1, y1, xmin, xmax, ymin, ymax);
		bool accept = false;

		while (true) {
			if (!(outcode0 | outcode1)) {
				// bitwise OR is 0: both points inside window; trivially accept and
				// exit loop
				accept = true;
				break;
			}
			else if (outcode0 & outcode1) {
				// bitwise AND is not 0: both points share an outside zone (LEFT,
				// RIGHT, TOP, or BOTTOM), so both must be outside window; exit loop
				// (accept is false)
				break;
			}
			else {
				// failed both tests, so calculate the line segment to clip
				// from an outside point to an intersection with clip edge
				double x, y;

				// At least one endpoint is outside the clip rectangle; pick it.
				OutCode outcodeOut = outcode0 ? outcode0 : outcode1;

				// Now find the intersection point;
				// use formulas:
				//   slope = (y1 - y0) / (x1 - x0)
				//   x = x0 + (1 / slope) * (ym - y0), where ym is ymin or ymax
				//   y = y0 + slope * (xm - x0), where xm is xmin or xmax
				// No need to worry about divide-by-zero because, in each case, the
				// outcode bit being tested guarantees the denominator is non-zero
				if (outcodeOut & TOP) { // point is above the clip window
					x = x0 + (x1 - x0) * (ymax - y0) / (y1 - y0);
					y = ymax;
				}
				else if (outcodeOut & BOTTOM) { // point is below the clip window
					x = x0 + (x1 - x0) * (ymin - y0) / (y1 - y0);
					y = ymin;
				}
				else if (outcodeOut &
					RIGHT) { // point is to the right of clip window
					y = y0 + (y1 - y0) * (xmax - x0) / (x1 - x0);
					x = xmax;
				}
				else if (outcodeOut &
					LEFT) { // point is to the left of clip window
					y = y0 + (y1 - y0) * (xmin - x0) / (x1 - x0);
					x = xmin;
				}

				// Now we move outside point to intersection point to clip
				// and get ready for next pass.
				if (outcodeOut == outcode0) {
					x0 = x;
					y0 = y;
					outcode0 = ComputeOutCode(x0, y0, xmin, xmax, ymin, ymax);
				}
				else {
					x1 = x;
					y1 = y;
					outcode1 = ComputeOutCode(x1, y1, xmin, xmax, ymin, ymax);
				}
			}
		}
		return accept;
	}

	bool LineString::intersects(const Envelope& rect) const {
		double xmin = rect.getMinX();
		double xmax = rect.getMaxX();
		double ymin = rect.getMinY();
		double ymax = rect.getMaxY();

		for (size_t i = 1; i < points.size(); ++i)
			if (intersectTest(points[i - 1].getX(), points[i - 1].getY(),
				points[i].getX(), points[i].getY(), xmin, xmax, ymin,
				ymax))
				return true;
		return false;
	}

	void LineString::draw() const {
		if (points.empty()) {
			std::cerr << "Warning: LineString points array is empty." << std::endl;
			return;
		}

		glBegin(GL_LINE_STRIP);
		for (size_t i = 0; i < points.size(); ++i)
			glVertex2d(points[i].getX(), points[i].getY());
		glEnd();

		GLenum err;
		while ((err = glGetError()) != GL_NO_ERROR) {
			std::cerr << "OpenGL Error in draw(): " << gluErrorString(err) << std::endl;
		}
	}

	void LineString::print() const {
		std::cout << "LineString(";
		for (size_t i = 0; i < points.size(); ++i) {
			if (i != 0)
				std::cout << ", ";
			std::cout << points[i].getX() << " " << points[i].getY();
		}
		std::cout << ")";
	}

	/*
	 * Polygon
	 */
	double Polygon::distance(const Polygon* polygon) const {
		return std::min(exteriorRing.distance(polygon),
			polygon->getExteriorRing().distance(this));
	}

	bool Polygon::intersects(const Envelope& rect) const {
		// TODO
		// 1. Polygon的MBR与rect相交
		if (!this->getEnvelope().intersect(rect)) {
			return false;
		}

		// 2. Polygon的边界（外环）与rect的边界相交
		if (exteriorRing.intersects(rect)) {
			return true;
		}

		// 3. Polygon包含rect (任一rect角点在Polygon内)
		double rx_min = rect.getMinX();
		double rx_max = rect.getMaxX();
		double ry_min = rect.getMinY();
		double ry_max = rect.getMaxY();

		Point p1(rx_min, ry_min);
		Point p2(rx_max, ry_min);
		Point p3(rx_max, ry_max);
		Point p4(rx_min, ry_max);

		// Point::distance(this) == 0.0 表示点在 Polygon 内部或边界
		if (p1.distance(this) == 0.0 || p2.distance(this) == 0.0 ||
			p3.distance(this) == 0.0 || p4.distance(this) == 0.0) {
			return true;
		}

		// 4. rect包含Polygon (任一Polygon点在rect内)
		// 检查Polygon的任一顶点是否在rect内部
		for (size_t i = 0; i < exteriorRing.numPoints(); ++i) {
			if (exteriorRing.getPointN(i).intersects(rect)) {
				return true;
			}
		}

		return false;
	}

	void Polygon::draw() const { exteriorRing.draw(); }

	void Polygon::print() const {
		std::cout << "Polygon(";
		for (size_t i = 0; i < exteriorRing.numPoints(); ++i) {
			if (i != 0)
				std::cout << ", ";
			Point p = exteriorRing.getPointN(i);
			std::cout << p.getX() << " " << p.getY();
		}
		std::cout << ")";
	}

} // namespace hw6
