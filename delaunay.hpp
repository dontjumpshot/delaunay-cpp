/* 
 * Bowyer-Watson algorithm
 * C++ implementation of http://paulbourke.net/papers/triangulate .
 **/
#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

namespace delaunay {

constexpr double eps = 1e-4;

/*创建通用的类或函数*/
template <typename T>
struct Point {
  T x, y;
  /*在定义变量时直接对其进行初始化*/
  Point() : x{0}, y{0} {}
  Point(T _x, T _y) : x{_x}, y{_y} {}

  /*用于把U类型转化为T类型*/
  template <typename U>
  Point(U _x, U _y) : x{static_cast<T>(_x)}, y{static_cast<T>(_y)}
  {
  }

  /*通过声明为friend，这个函数可以访问Point结构体的私有成员变量x和y，并直接输出它们的值。*/
  friend std::ostream& operator<<(std::ostream& os, const Point<T>& p)
  {
    os << "x=" << p.x << "  y=" << p.y;
    return os;
  }

  /*在Point结构体中定义了一个相等比较运算符的函数重载，用于比较两个Point对象是否相等*/
  bool operator==(const Point<T>& other) const
  {
    return (other.x == x && other.y == y);
  }

  /*在Point结构体中定义了一个不等比较运算符的函数重载，用于比较两个Point对象是否不相等*/
  bool operator!=(const Point<T>& other) const { return !operator==(other); }
};

/*定义了一个名为Edge的结构体，表示一条边。其中包含了一个内部类型别名Node，表示结构体Point<T>*/
template <typename T>
struct Edge {
  using Node = Point<T>;
  Node p0, p1;
  /*初始化*/
  Edge(Node const& _p0, Node const& _p1) : p0{_p0}, p1{_p1} {}
  
  friend std::ostream& operator<<(std::ostream& os, const Edge& e)
  {
    os << "p0: [" << e.p0 << " ] p1: [" << e.p1 << "]";
    return os;
  }

  /*两个边是否相等*/
  bool operator==(const Edge& other) const
  {
    return ((other.p0 == p0 && other.p1 == p1) ||
            (other.p0 == p1 && other.p1 == p0));
  }
};

/*Circle结构体有三个成员变量：x、y和radius，分别表示圆的中心点的 x 坐标、y 坐标以及半径。
构造函数使用了默认的成员初始化方式，即使用= default;进行声明。这意味着在创建Circle对象时，这些成员变量会被默认初始化。*/
template <typename T>
struct Circle {
  T x, y, radius;
  Circle() = default;
};

/*Triangle结构体使用了内部类型别名Node，表示结构体Point<T>。
Triangle结构体有六个成员变量：p0、p1、p2、e0、e1、e2和circle。其中p0、p1、p2是三角形的三个顶点，类型为Node（即Point<T>类型）；
e0、e1、e2是三角形的三条边，类型为Edge<T>；circle是三角形的外接圆，类型为Circle<T>。*/
template <typename T>
struct Triangle {
  using Node = Point<T>;
  Node p0, p1, p2;
  Edge<T> e0, e1, e2;
  Circle<T> circle;

  /*定义了一个构造函数，用于初始化Triangle结构体对象。构造函数接受三个顶点参数，并根据这些参数计算并初始化三角形的边和外接圆信息*/
  Triangle(const Node& _p0, const Node& _p1, const Node& _p2)
      : p0{_p0},
        p1{_p1},
        p2{_p2},
        e0{_p0, _p1},
        e1{_p1, _p2},
        e2{_p0, _p2},
        circle{}
  {
    const auto ax = p1.x - p0.x;
    const auto ay = p1.y - p0.y;
    const auto bx = p2.x - p0.x;
    const auto by = p2.y - p0.y;

    const auto m = p1.x * p1.x - p0.x * p0.x + p1.y * p1.y - p0.y * p0.y;
    const auto u = p2.x * p2.x - p0.x * p0.x + p2.y * p2.y - p0.y * p0.y;
    const auto s = 1. / (2. * (ax * by - ay * bx));

    circle.x = ((p2.y - p0.y) * m + (p0.y - p1.y) * u) * s;
    circle.y = ((p0.x - p2.x) * m + (p1.x - p0.x) * u) * s;

    const auto dx = p0.x - circle.x;
    const auto dy = p0.y - circle.y;
    circle.radius = dx * dx + dy * dy;
  }
};

/*Delaunay结构体包含两个成员变量：triangles和edges，分别是std::vector<Triangle<T>>和std::vector<Edge<T>>类型的向量。
triangles向量用于存储Delaunay三角剖分的三角形对象，每个三角形对象的类型为Triangle<T>，其中T是模板参数，表示顶点的数据类型。
edges向量用于存储Delaunay三角剖分的边对象，每个边对象的类型为Edge<T>。*/
template <typename T>
struct Delaunay {
  std::vector<Triangle<T>> triangles;
  std::vector<Edge<T>> edges;
};

/*用于确保模板参数T必须是浮点数类型才能使用*/
template <
    typename T,
    typename = typename std::enable_if<std::is_floating_point<T>::value>::type>

/*定义了一个名为triangulate的函数模板，用于进行Delaunay三角剖分。*/
Delaunay<T> triangulate(const std::vector<Point<T>>& points)
{
  using Node = Point<T>;
  if (points.size() < 3) {
    return Delaunay<T>{};
  }
  auto xmin = points[0].x;
  auto xmax = xmin;
  auto ymin = points[0].y;
  auto ymax = ymin;
  /*找出点集的x,y范围*/
  for (auto const& pt : points) {
    xmin = std::min(xmin, pt.x);
    xmax = std::max(xmax, pt.x);
    ymin = std::min(ymin, pt.y);
    ymax = std::max(ymax, pt.y);
  }

  const auto dx = xmax - xmin;
  const auto dy = ymax - ymin;
  const auto dmax = std::max(dx, dy);
  const auto midx = (xmin + xmax) / static_cast<T>(2.);
  const auto midy = (ymin + ymax) / static_cast<T>(2.);

  /* Init Delaunay triangulation. */
  auto d = Delaunay<T>{};  /*定义了一个 Delaunay<T> 类型的变量 d，并对其进行了值初始化*/
  
  /*通过计算坐标并创建节点，将这些节点组成的三角形添加到 Delaunay<T> 对象的三角形列表中*/
  const auto p0 = Node{midx - 20 * dmax, midy - dmax};
  const auto p1 = Node{midx, midy + 20 * dmax};
  const auto p2 = Node{midx + 20 * dmax, midy - dmax};
  d.triangles.emplace_back(Triangle<T>{p0, p1, p2});

  /*检查计算得到的距离的平方减去三角形圆的半径是否小于等于 eps，
如果满足条件，说明当前点位于三角形的外接圆内或边界上，它将该三角形的三条边添加到 edges（非德劳内三角形） 中。否则（即当前点位于三角形外接圆外），将该三角形添加到 tmps 中。
循环结束后，edges 中存储了与当前点相关的边，而 tmps 中存储了剩余的三角形。
*/
  for (auto const& pt : points) {
    std::vector<Edge<T>> edges;  \\存储边
    std::vector<Triangle<T>> tmps;  \\存储临时三角形
    for (auto const& tri : d.triangles) {
      /* 检查该点是否在三角形外接圆内. */
      const auto dist = (tri.circle.x - pt.x) * (tri.circle.x - pt.x) +
                        (tri.circle.y - pt.y) * (tri.circle.y - pt.y);
      if ((dist - tri.circle.radius) <= eps) {
        edges.push_back(tri.e0);
        edges.push_back(tri.e1);
        edges.push_back(tri.e2);
      }
      else {
        tmps.push_back(tri);
      }
    }

    /* 删除重复边 */
    std::vector<bool> remove(edges.size(), false);
    for (auto it1 = edges.begin(); it1 != edges.end(); ++it1) {
      for (auto it2 = edges.begin(); it2 != edges.end(); ++it2) {
        if (it1 == it2) {
          continue;
        }
        if (*it1 == *it2) {
          remove[std::distance(edges.begin(), it1)] = true;
          remove[std::distance(edges.begin(), it2)] = true;
        }
      }
    }

    edges.erase(
        std::remove_if(edges.begin(), edges.end(),
                       [&](auto const& e) { return remove[&e - &edges[0]]; }),
        edges.end());

    /* 更新三角剖分 */
    for (auto const& e : edges) {
      tmps.push_back({e.p0, e.p1, {pt.x, pt.y}});
    }
    d.triangles = tmps;
  }

  /* 去掉原来的超三角形 */
  d.triangles.erase(
      std::remove_if(d.triangles.begin(), d.triangles.end(),
                     [&](auto const& tri) {
                       return ((tri.p0 == p0 || tri.p1 == p0 || tri.p2 == p0) ||
                               (tri.p0 == p1 || tri.p1 == p1 || tri.p2 == p1) ||
                               (tri.p0 == p2 || tri.p1 == p2 || tri.p2 == p2));
                     }),
      d.triangles.end());

  /* 添加边 */
  for (auto const& tri : d.triangles) {
    d.edges.push_back(tri.e0);
    d.edges.push_back(tri.e1);
    d.edges.push_back(tri.e2);
  }
  return d;
}

} /* namespace delaunay */
