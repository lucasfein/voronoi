#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

#define PI 3.1415926535897932
#define EPSILON 0.0025

#define XMLNS "http://graphml.graphdrawing.org/xmlns"
#define XMLNS_XSI "http://www.w3.org/2001/XMLSchema-instance"
#define XSI_SL                                                                 \
  "http://graphml.graphdrawing.org/xmlns "                                     \
  "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd"

double deg(double rad) {
  return (rad > 0 ? rad : (2 * PI + rad)) * 360 / (2 * PI);
}

struct Point {
  double x;
  double y;
  double r;

  Point(double x = 0.0, double y = 0.0, double r = 0.0) : x(x), y(y), r(r) {}

  Point(std::vector<double> v, double r = 1.0) : x(v[0]), y(v[1]), r(r) {}

  bool operator==(const Point &point) const {
    return abs(point.x - x) < EPSILON && abs(point.y - y) < EPSILON;
  }
  bool operator<(const Point &point) const {
    return (x < point.x || (abs(point.x - x) < EPSILON && y < point.y));
  }
  double distance(const Point &point) const {
    return sqrt(pow(x - point.x, 2) + pow(y - point.y, 2));
  }
  double distance(const std::vector<double> &point) const {
    return sqrt(pow(x - point[0], 2) + pow(y - point[1], 2));
  }
  std::vector<double> vector_from(const Point &point_b) const {
    return std::vector<double>{x - point_b.x, y - point_b.y};
  }
  std::vector<double> vector_from(const std::vector<double> &point_b) const {
    return std::vector<double>{x - point_b[0], y - point_b[1]};
  }

  std::vector<double> unit_vector_from(const Point &point_b) const {
    double dis = distance(point_b);
    return std::vector<double>{(x - point_b.x) / dis, (y - point_b.y) / dis};
  }
  std::vector<double>
  unit_vector_from(const std::vector<double> &point_b) const {
    double dis = distance(point_b);
    return std::vector<double>{(x - point_b[0]) / dis, (y - point_b[1]) / dis};
  }
  std::vector<double> pos() { return std::vector<double>{x, y}; }

  double cross_product(const Point &a, const Point &b) {
    return (b.x - a.x) * (y - a.y) - (b.y - a.y) * (x - a.x);
  }
};

struct Edge {
  Point a;
  Point b;

  Edge(Point a = Point(), Point b = Point()) : a(a), b(b) {}

  bool operator==(const Edge &edge) const {
    return (edge.a == a && edge.b == b) || (edge.a == b && edge.b == a);
  }
  bool contains(const Point &point) const { return (point == a || point == b); }
  double length() const { return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2)); }
  std::vector<double> midpoint() const {
    return std::vector<double>{a.x + (b.x - a.x) * 0.5,
                               a.y + (b.y - a.y) * 0.5};
  }
  std::vector<double> weighted_midpoint() const {
    double s = (a.r - b.r) / a.distance(b);
    return std::vector<double>{a.x + (b.x - a.x) * (0.5 + s),
                               a.y + (b.y - a.y) * (0.5 + s)};
  }
  double rotation() const {
    std::vector<double> mid_point = midpoint();
    double temp_rad = atan2(a.y - mid_point[1], a.x - mid_point[0]);
    return temp_rad * 180 / PI;
  }
  bool shares_point_with(const Edge &edge) const {
    return (edge.a == a || edge.a == b || edge.b == a || edge.b == b);
  }
};

struct Triangle {
  Point a;
  Point b;
  Point c;

  Triangle(Point point_a = Point(), Point point_b = Point(),
           Point point_c = Point())
      : a(point_a), b(point_b), c(point_c) {
    sort_points_counter_clockwise();
  }

  Triangle(Edge edge = Edge(), Point point = Point())
      : a(edge.a), b(edge.b), c(point) {
    sort_points_counter_clockwise();
  }

  bool operator==(const Triangle &triangle) const {
    return (triangle.a == a && triangle.b == b && triangle.c == c) ||
           (triangle.a == a && triangle.b == c && triangle.c == b) ||
           (triangle.a == b && triangle.b == a && triangle.c == c) ||
           (triangle.a == b && triangle.b == c && triangle.c == a) ||
           (triangle.a == c && triangle.b == a && triangle.c == b) ||
           (triangle.a == c && triangle.b == b && triangle.c == a);
  }
  bool contains(const Point &point) const {
    return (a == point || b == point || c == point);
  }
  bool contains(const Edge &edge) const {
    return ((edge.contains(a) && edge.contains(b)) ||
            (edge.contains(b) && edge.contains(c)) ||
            (edge.contains(c) && edge.contains(a)));
  }

  bool triangle_contains(const Point &point) {
    double a_pab = abs(point.x * (a.y - b.y) + a.x * (b.y - point.y) +
                       b.x * (point.y - a.y)) *
                   0.5;
    double a_pbc = abs(point.x * (b.y - c.y) + b.x * (c.y - point.y) +
                       c.x * (point.y - b.y)) *
                   0.5;
    double a_pca = abs(point.x * (c.y - a.y) + c.x * (a.y - point.y) +
                       a.x * (point.y - c.y)) *
                   0.5;
    double a_abc =
        abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) * 0.5;
    return (abs((a_pab + a_pbc + a_pca) - a_abc) < EPSILON);
  }
  bool circumcircle_contains(const Point &point) const {
    double matrix[3][3];
    matrix[0][0] = a.x - point.x;
    matrix[0][1] = a.y - point.y;
    matrix[0][2] = pow(a.x - point.x, 2) + pow(a.y - point.y, 2);
    matrix[1][0] = b.x - point.x;
    matrix[1][1] = b.y - point.y;
    matrix[1][2] = pow(b.x - point.x, 2) + pow(b.y - point.y, 2);
    matrix[2][0] = c.x - point.x;
    matrix[2][1] = c.y - point.y;
    matrix[2][2] = pow(c.x - point.x, 2) + pow(c.y - point.y, 2);

    return (matrix[0][0] *
                (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
            matrix[0][1] *
                (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
            matrix[0][2] * (matrix[1][0] * matrix[2][1] -
                            matrix[1][1] * matrix[2][0])) > 0;
  }
  bool shares_vertex_with(const Triangle &triangle) const {
    return (triangle.a == a || triangle.a == b || triangle.a == c) ||
           (triangle.b == a || triangle.b == b || triangle.b == c) ||
           (triangle.c == a || triangle.c == b || triangle.c == c);
  }
  bool shares_edge_with(const Triangle &triangle) const {
    return (contains(Edge(triangle.a, triangle.b)) ||
            contains(Edge(triangle.b, triangle.c)) ||
            contains(Edge(triangle.c, triangle.a)));
  }
  Edge get_common_edge(const Triangle &triangle) const {
    if (triangle.contains(a) && triangle.contains(b)) {
      return Edge(a, b);
    } else if (triangle.contains(b) && triangle.contains(c)) {
      return Edge(b, c);
    } else {
      return Edge(c, a);
    }
  }

  std::vector<double> center() const {
    return std::vector<double>{(a.x + b.x + c.x) / 3.0,
                               (a.y + b.y + c.y) / 3.0};
  }
  std::vector<double> circumcenter() const {
    double temp;
    std::vector<double> midpoint_a = Edge(a, b).midpoint();
    std::vector<double> midpoint_b = Edge(b, c).midpoint();

    double a1 = b.y - a.y;
    double b1 = a.x - b.x;
    double c1 = -b1 * midpoint_a[0] + a1 * midpoint_a[1];
    temp = a1;
    a1 = -b1;
    b1 = temp;

    double a2 = c.y - b.y;
    double b2 = b.x - c.x;
    double c2 = -b2 * midpoint_b[0] + a2 * midpoint_b[1];
    temp = a2;
    a2 = -b2;
    b2 = temp;

    double det = a1 * b2 - a2 * b1;

    double x = (b2 * c1 - b1 * c2) / det;
    double y = (a1 * c2 - a2 * c1) / det;

    return std::vector<double>{x, y};
  }
  void sort_points_counter_clockwise() {
    std::vector<double> triangle_center = center();

    double angle_a =
        atan2(a.x - triangle_center[0], -(a.y - triangle_center[1]));
    double angle_b =
        atan2(b.x - triangle_center[0], -(b.y - triangle_center[1]));
    double angle_c =
        atan2(c.x - triangle_center[0], -(c.y - triangle_center[1]));

    std::vector<Point> points;
    if (angle_a < angle_b && angle_b < angle_c) {
      points = std::vector<Point>{a, b, c};
    } else if (angle_a < angle_c && angle_c < angle_b) {
      points = std::vector<Point>{a, c, b};
    } else if (angle_b < angle_a && angle_a < angle_c) {
      points = std::vector<Point>{b, a, c};
    } else if (angle_b < angle_c && angle_c < angle_a) {
      points = std::vector<Point>{b, c, a};
    } else if (angle_c < angle_a && angle_a < angle_b) {
      points = std::vector<Point>{c, a, b};
    } else if (angle_c < angle_b && angle_b < angle_a) {
      points = std::vector<Point>{c, b, a};
    } else {
      points = std::vector<Point>{a, b, c};
    }
    a = points[0];
    b = points[1];
    c = points[2];
  }

  std::vector<Edge> get_edges() {
    return std::vector<Edge>{Edge(a, b), Edge(b, c), Edge(c, a)};
  }

  void translate(std::vector<double> delta) {
    a.x = a.x + delta[0];
    a.y = a.y + delta[1];
    b.x = b.x + delta[0];
    b.y = b.y + delta[1];
    c.x = c.x + delta[0];
    c.y = c.y + delta[1];
  }

  void scale(std::vector<double> scale, std::vector<double> cntr) {
    a.x = cntr[0] + scale[0] * (a.x - cntr[0]);
    a.y = cntr[1] + scale[1] * (a.y - cntr[1]);
    b.x = cntr[0] + scale[0] * (b.x - cntr[0]);
    b.y = cntr[1] + scale[1] * (b.y - cntr[1]);
    c.x = cntr[0] + scale[0] * (c.x - cntr[0]);
    c.y = cntr[1] + scale[1] * (c.y - cntr[1]);
  }

  void rotate(double rad, std::vector<double> cntr) {
    double deg = rad * PI / 180;
    double a_x =
        (a.x - cntr[0]) * cos(deg) + (a.y - cntr[1]) * (-sin(deg)) + cntr[0];
    double a_y =
        (a.x - cntr[0]) * sin(deg) + (a.y - cntr[1]) * cos(deg) + cntr[1];
    double b_x =
        (b.x - cntr[0]) * cos(deg) + (b.y - cntr[1]) * (-sin(deg)) + cntr[0];
    double b_y =
        (b.x - cntr[0]) * sin(deg) + (b.y - cntr[1]) * cos(deg) + cntr[1];
    double c_x =
        (c.x - cntr[0]) * cos(deg) + (c.y - cntr[1]) * (-sin(deg)) + cntr[0];
    double c_y =
        (c.x - cntr[0]) * sin(deg) + (c.y - cntr[1]) * cos(deg) + cntr[1];

    a = Point(a_x, a_y);
    b = Point(b_x, b_y);
    c = Point(c_x, c_y);
  }
};

struct Polygon {
  std::vector<Point> vertices;

  Polygon() { vertices = std::vector<Point>{}; }

  Polygon(std::vector<Point> nodes) : vertices(nodes) {}

  Point centroid() const {
    double temp;
    double temp_area = 0.0;
    double temp_x = 0.0;
    double temp_y = 0.0;
    for (unsigned int i = 0; i < vertices.size() - 1; i++) {
      temp =
          vertices[i].x * vertices[i + 1].y - vertices[i + 1].x * vertices[i].y;
      temp_area += temp;
      temp_x += (vertices[i].x + vertices[i + 1].x) * temp;
      temp_y += (vertices[i].y + vertices[i + 1].y) * temp;
    }
    return Point(temp_x / (3 * temp_area), temp_y / (3 * temp_area));
  }

  double area() const {
    double temp = 0.0;
    for (unsigned int i = 0; i < vertices.size() - 1; i++) {
      temp +=
          vertices[i].x * vertices[i + 1].y - vertices[i + 1].x * vertices[i].y;
    }
    return abs(temp * 0.5);
  }

  bool operator<(const Polygon &polygon) const {
    return ((polygon.area() < area() ? 1 : 0));
  }
};

struct Generators {
  std::vector<Point> generators;

  Generators(std::vector<Point> generators = std::vector<Point>{})
      : generators(generators) {}

  Generators(std::string filename) {
    std::ifstream csv(filename);
    std::string x, y, r;
    while (getline(csv, x, ',') && getline(csv, y, ',') && getline(csv, r)) {
      generators.push_back(Point(std::stod(x), std::stod(y), std::stod(r)));
    }
    csv.close();
  }

  void save(const std::string &fileName) {
    std::ofstream csv;
    csv.open(fileName);
    for (unsigned long i = 0; i < generators.size(); i++) {
      csv << generators[i].x << "," << generators[i].y << "," << generators[i].r
          << "\n";
    }
    csv.close();
  }
};

struct DelaunayTriangulation {
  int dim_x, dim_y;
  Generators generators;
  std::vector<Point> vertices;
  std::vector<Edge> edges;
  std::vector<Triangle> triangles;
  std::vector<std::vector<unsigned long>> triangle_edge_map;
  std::vector<std::pair<unsigned int, unsigned int>> edge_vertex_map;
  std::vector<std::vector<bool>> adj;
  std::vector<std::vector<double>> vertex_angle_map;

  DelaunayTriangulation(Generators generators = Generators(), int dim_x = 1000,
                        int dim_y = 1000)
      : dim_x(dim_x), dim_y(dim_y), generators(generators) {
    // Bower-Watson-Algorithm
    Triangle super_triangle =
        Triangle(Point(-3 * dim_x, -3 * dim_y), Point(3 * dim_x, -3 * dim_y),
                 Point(0.5 * dim_x, 3 * dim_y));
    triangles = std::vector<Triangle>{super_triangle};

    for (std::vector<Point>::iterator generator = generators.generators.begin();
         generator != generators.generators.end(); generator++) {
      std::vector<Triangle> bad_triangles = std::vector<Triangle>{};
      for (std::vector<Triangle>::iterator triangle = triangles.begin();
           triangle != triangles.end(); triangle++) {
        if (triangle->circumcircle_contains(*generator)) {
          bad_triangles.push_back(*triangle);
        }
      }
      std::vector<Edge> polygon = std::vector<Edge>{};
      for (std::vector<Triangle>::iterator bad_triangle = bad_triangles.begin();
           bad_triangle != bad_triangles.end(); bad_triangle++) {
        std::vector<Edge> bad_triangle_edges = bad_triangle->get_edges();
        for (std::vector<Edge>::iterator edge = bad_triangle_edges.begin();
             edge != bad_triangle_edges.end(); edge++) {
          bool unique_edge = 1;
          for (std::vector<Triangle>::iterator triangle = bad_triangles.begin();
               triangle != bad_triangles.end(); triangle++) {
            if (triangle == bad_triangle) {
              continue;
            } else if (triangle->contains(*edge)) {
              unique_edge = 0;
              break;
            }
          }
          if (unique_edge) {
            polygon.push_back(*edge);
          }
        }
      }
      for (std::vector<Triangle>::iterator bad_triangle = bad_triangles.begin();
           bad_triangle != bad_triangles.end(); bad_triangle++) {
        for (std::vector<Triangle>::iterator triangle = triangles.begin();
             triangle != triangles.end(); triangle++) {
          if (*bad_triangle == *triangle) {
            triangles.erase(triangle);
            break;
          }
        }
      }
      for (std::vector<Edge>::iterator edge = polygon.begin();
           edge != polygon.end(); edge++) {
        Triangle triangle = Triangle(*edge, *generator);
        triangles.push_back(triangle);
      }
    }
    for (std::vector<Triangle>::iterator triangle = triangles.begin();
         triangle != triangles.end();) {
      if (triangle->shares_vertex_with(super_triangle)) {
        triangles.erase(triangle);
      } else {
        ++triangle;
      }
    }

    vertices = generators.generators;
    adj = std::vector<std::vector<bool>>(vertices.size(),
                                         std::vector<bool>(vertices.size(), 0));
    triangle_edge_map =
        std::vector<std::vector<unsigned long>>(triangles.size());

    for (unsigned long i = 0; i < vertices.size(); i++) {
      for (unsigned long j = 0; j < i; j++) {
        for (unsigned long t = 0; t < triangles.size(); t++) {
          if (triangles[t].contains(vertices[i]) &&
              triangles[t].contains(vertices[j])) {
            adj[i][j] = 1;
            adj[j][i] = 1;
          }
        }
      }
    }
    for (unsigned long i = 0; i < vertices.size(); i++) {
      for (unsigned long j = 0; j < i; j++) {
        if (adj[i][j]) {
          edges.push_back(Edge(vertices[i], vertices[j]));
          edge_vertex_map.push_back(std::make_pair(i, j));
          for (unsigned int t = 0; t < triangles.size(); t++) {
            if (triangles[t].contains(Edge(vertices[i], vertices[j]))) {
              triangle_edge_map[t].push_back(edges.size() - 1);
            }
          }
        }
      }
    }
    double det, dot, angle, rad;
    std::vector<double> vec_a, vec_b;
    std::vector<std::vector<std::pair<unsigned int, double>>> adj_vertices(
        vertices.size(), std::vector<std::pair<unsigned int, double>>{});
    for (unsigned int i = 0; i < vertices.size(); i++) {
      for (unsigned int j = 0; j < vertices.size(); j++) {
        if (adj[i][j]) {
          rad = -atan2(vertices[j].x - vertices[i].x,
                       -(vertices[j].y - vertices[i].y));
          rad = (rad < 0) ? rad + 2 * PI : rad;
          adj_vertices[i].push_back(std::make_pair(j, rad));
        }
      }
      std::sort(adj_vertices[i].begin(), adj_vertices[i].end(),
                [](std::pair<unsigned int, double> a,
                   std::pair<unsigned int, double> b) {
                  return a.second < b.second;
                });
      vertex_angle_map.push_back(std::vector<double>(1));
      if (adj_vertices[i].size() > 2) {
        vec_a = vertices[adj_vertices[i][0].first].vector_from(vertices[i]);
        vec_b = vertices[adj_vertices[i].back().first].vector_from(vertices[i]);
        dot = vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1];
        det = vec_a[0] * vec_b[1] - vec_b[0] * vec_a[1];
        angle = atan2(det, dot);
        angle = (angle < 0) ? angle + 2 * PI : angle;
        vertex_angle_map[i][0] = angle;
        for (int k = adj_vertices[i].size() - 1; k > 0; k--) {
          vec_a = vertices[adj_vertices[i][k].first].vector_from(vertices[i]);
          vec_b =
              vertices[adj_vertices[i][k - 1].first].vector_from(vertices[i]);
          dot = vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1];
          det = vec_a[0] * vec_b[1] - vec_b[0] * vec_a[1];
          angle = atan2(det, dot);
          angle = (angle < 0) ? angle + 2 * PI : angle;
          vertex_angle_map[i].push_back(angle);
        }
      } else if (adj_vertices[i].size() == 2) {
        vec_a = vertices[adj_vertices[i][0].first].vector_from(vertices[i]);
        vec_b = vertices[adj_vertices[i][1].first].vector_from(vertices[i]);
        dot = vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1];
        det = vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];
        angle = atan2(det, dot);
        angle = (angle < 0) ? angle + 2 * PI : angle;
        vertex_angle_map[i][0] = angle;
        vertex_angle_map[i].push_back(2 * PI - angle);
      } else if (adj_vertices[i].size() == 1) {
        vertex_angle_map[i][0] = 2 * PI;
      } else {
        continue;
      }
    }
  }

  bool triangulation_contains(std::vector<double> pos) {
    for (unsigned int i = 0; i < triangles.size(); i++) {
      if (triangles[i].triangle_contains(pos)) {
        return 1;
      }
    }
    return 0;
  }

  void save(const std::string &fileName) {
    std::ofstream graphml;
    graphml.open(fileName);
    graphml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    int indentation = 0;
    graphml << open_tag(
        "graphml",
        std::map<std::string, std::string>{{"xmlns", XMLNS},
                                           {"xmlns:xsi", XMLNS_XSI},
                                           {"xsi:schemaLocation", XSI_SL}},
        indentation, 0, 1);
    indentation++;
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "x"},
                                           {"attr.type", "double"},
                                           {"for", "node"},
                                           {"id", "d0"}},
        indentation, 1, 1);
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "y"},
                                           {"attr.type", "double"},
                                           {"for", "node"},
                                           {"id", "d1"}},
        indentation, 1, 1);

    unsigned int max_num_angles = 0;
    for (unsigned int i = 0; i < vertices.size(); i++) {
      if (vertex_angle_map[i].size() > max_num_angles) {
        max_num_angles = vertex_angle_map[i].size();
      }
    }
    for (unsigned int i = 0; i < max_num_angles; i++) {
      graphml << open_tag("key",
                          std::map<std::string, std::string>{
                              {"attr.name", "angle" + std::to_string(i)},
                              {"attr.type", "double"},
                              {"for", "node"},
                              {"id", "d2" + std::to_string(i)}},
                          indentation, 1, 1);
    }
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "length"},
                                           {"attr.type", "double"},
                                           {"for", "edge"},
                                           {"id", "d3"}},
        indentation, 1, 1);
    graphml << open_tag(
        "graph",
        std::map<std::string, std::string>{{"edgedefault", "undirected"},
                                           {"dim_x", std::to_string(dim_x)},
                                           {"dim_y", std::to_string(dim_y)}},
        indentation, 0, 1);
    indentation++;

    for (unsigned long i = 0; i < vertices.size(); i++) {
      graphml << open_tag(
          "node", std::map<std::string, std::string>{{"id", std::to_string(i)}},
          indentation, 0, 1);
      indentation++;
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d0"}},
                          indentation, 0, 0);
      graphml << vertices[i].x;
      graphml << close_tag("data", indentation, 1);
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d1"}},
                          indentation, 0, 0);
      graphml << vertices[i].y;
      graphml << close_tag("data", indentation, 1);

      for (unsigned int j = 0; j < vertex_angle_map[i].size(); j++) {
        graphml << open_tag("data",
                            std::map<std::string, std::string>{
                                {"key", "d2" + std::to_string(j)}},
                            indentation, 0, 0);
        graphml << vertex_angle_map[i][j];
        graphml << close_tag("data", indentation, 1);
      }

      indentation--;
      graphml << close_tag("node", indentation);
    }

    for (unsigned long i = 0; i < edges.size(); i++) {
      graphml << open_tag(
          "edge",
          std::map<std::string, std::string>{
              {"source", std::to_string(std::min(edge_vertex_map[i].first,
                                                 edge_vertex_map[i].second))},
              {"target", std::to_string(std::max(edge_vertex_map[i].first,
                                                 edge_vertex_map[i].second))}},
          indentation, 0, 1);
      indentation++;
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d3"}},
                          indentation, 0, 0);
      graphml << edges[i].length();
      graphml << close_tag("data", indentation, 1);
      indentation--;
      graphml << close_tag("edge", indentation);
    }

    indentation--;
    graphml << close_tag("graph", indentation);
    indentation--;
    graphml << close_tag("graphml", indentation);
    graphml.close();
  }
  std::string open_tag(std::string name,
                       std::map<std::string, std::string> attributes,
                       int indentation, bool closed = 0, bool newline = 0) {
    std::stringstream graphml;
    for (int i = 0; i < indentation; i++) {
      graphml << "\t";
    }
    graphml << "<" << name;
    for (auto const &attribute : attributes) {
      graphml << " " << attribute.first << "=\"" << attribute.second << "\"";
    }
    if (closed) {
      graphml << "/";
    }
    graphml << ">";
    if (newline) {
      graphml << "\n";
    }
    return graphml.str();
  }
  std::string close_tag(std::string name, int indentation, bool in_line = 0) {
    std::stringstream graphml;
    if (!in_line) {
      for (int i = 0; i < indentation; i++) {
        graphml << "\t";
      }
    }
    graphml << "</" << name << ">\n";
    return graphml.str();
  }
};

struct VoronoiTessellation {
  int dim_x, dim_y;
  Generators generators;
  std::vector<Point> vertices;
  std::vector<Edge> edges;
  std::vector<std::pair<unsigned long, unsigned long>> edge_vertex_map;
  std::vector<Polygon> polygons;
  std::vector<std::vector<unsigned long>> polygon_vertex_map;
  std::vector<std::vector<double>> vertex_angle_map;
  std::vector<std::vector<bool>> adj;

  void voronoi_from_delaunay(DelaunayTriangulation delaunay_triangulation) {
    std::vector<std::vector<double>> cc_dis(
        delaunay_triangulation.triangles.size(),
        std::vector<double>(delaunay_triangulation.triangles.size(), 0.0));
    std::vector<std::vector<double>> circumcenters(
        delaunay_triangulation.triangles.size(), std::vector<double>{0.0, 0.0});
    for (unsigned int i = 0; i < delaunay_triangulation.triangles.size(); i++) {
      circumcenters[i] = delaunay_triangulation.triangles[i].circumcenter();
    }
    double dis;
    for (unsigned int i = 0; i < delaunay_triangulation.triangles.size(); i++) {
      for (unsigned int j = 0; j < i; j++) {
        dis = sqrt(pow(circumcenters[j][0] - circumcenters[i][0], 2) +
                   pow(circumcenters[j][1] - circumcenters[i][1], 2));
        cc_dis[i][j] = dis;
        cc_dis[j][i] = dis;
      }
    }
    for (unsigned long i = 0; i < delaunay_triangulation.triangles.size();
         i++) {
      std::vector<Edge> edges = delaunay_triangulation.triangles[i].get_edges();
      std::vector<double> vec, partial_shift,
          total_shift = std::vector<double>{0.0, 0.0};
      std::vector<double> cc = circumcenters[i];
      double ratio;
      for (unsigned int k = 0; k < 3; k++) {
        vec = edges[k].b.vector_from(edges[k].a);
        ratio = (edges[k].a.r - edges[k].b.r) / edges[k].length();
        partial_shift = std::vector<double>{ratio * vec[0], ratio * vec[1]};
        total_shift[0] += partial_shift[0] / 3;
        total_shift[1] += partial_shift[1] / 3;
      }
      for (unsigned int j = 0; j < i; j++) {
        if (i != j && cc_dis[i][j] < sqrt(pow(total_shift[0], 2) +
                                          pow(total_shift[1], 2))) {
          total_shift[0] = 0.0;
          total_shift[1] = 0.0;
          break;
        }
      }
      vertices.push_back(Point(cc[0] + total_shift[0], cc[1] + total_shift[1]));
    }

    std::vector<double> vec;
    std::vector<std::vector<int>> periphery_vertices(
        delaunay_triangulation.triangles.size(), std::vector<int>(3, -1));
    double delta_x, delta_y, scale;
    bool vertex_i_outside, vertex_j_outside;
    for (unsigned long i = 0; i < delaunay_triangulation.triangles.size();
         i++) {
      for (unsigned int j = 0; j < i; j++) {
        if (delaunay_triangulation.triangles[j].shares_edge_with(
                delaunay_triangulation.triangles[i])) {
          vertex_i_outside = (vertices[i].x < 0 || dim_x < vertices[i].x ||
                              vertices[i].y < 0 || dim_y < vertices[i].y);
          vertex_j_outside = (vertices[j].x < 0 || dim_x < vertices[j].x ||
                              vertices[j].y < 0 || dim_y < vertices[j].y);
          if (vertex_i_outside && vertex_j_outside) {
            continue;
          } else if (vertex_i_outside && !vertex_j_outside) {
            vec = vertices[i].vector_from(vertices[j]);
            delta_x = (0 < vec[0]) ? dim_x - vertices[j].x : vertices[j].x;
            delta_y = (0 < vec[1]) ? dim_y - vertices[j].y : vertices[j].y;
            scale = (delta_x / abs(vec[0]) < delta_y / abs(vec[1]))
                        ? delta_x / abs(vec[0])
                        : delta_y / abs(vec[1]);
            Point point =
                Point(abs(std::round(vertices[j].x + vec[0] * scale)),
                      abs(std::round(vertices[j].y + vec[1] * scale)));
            edges.push_back(Edge(vertices[j], point));
            vertices.push_back(point);
            edge_vertex_map.push_back(std::make_pair(j, vertices.size() - 1));
            for (int k = 0; k < 3; k++) {
              if (delaunay_triangulation.triangles[j].contains(
                      delaunay_triangulation.triangles[i].get_edges()[k])) {
                periphery_vertices[i][k] = vertices.size() - 1;
              }
            }
          } else if (!vertex_i_outside && vertex_j_outside) {
            vec = vertices[j].vector_from(vertices[i]);
            delta_x = (0 < vec[0]) ? dim_x - vertices[i].x : vertices[i].x;
            delta_y = (0 < vec[1]) ? dim_y - vertices[i].y : vertices[i].y;
            scale = (delta_x / abs(vec[0]) < delta_y / abs(vec[1]))
                        ? delta_x / abs(vec[0])
                        : delta_y / abs(vec[1]);
            Point point =
                Point(abs(std::round(vertices[i].x + vec[0] * scale)),
                      abs(std::round(vertices[i].y + vec[1] * scale)));
            edges.push_back(Edge(vertices[i], point));
            vertices.push_back(point);
            edge_vertex_map.push_back(std::make_pair(i, vertices.size() - 1));
            for (int k = 0; k < 3; k++) {
              if (delaunay_triangulation.triangles[i].contains(
                      delaunay_triangulation.triangles[j].get_edges()[k])) {
                periphery_vertices[j][k] = vertices.size() - 1;
              }
            }
          } else {
            edges.push_back(Edge(vertices[i], vertices[j]));
            edge_vertex_map.push_back(std::make_pair(i, j));
          }
        }
      }
    }

    bool edge_shared_by_triangles = 0;
    std::vector<double> bisector, a_b;
    std::vector<Edge> triangle_edges;
    for (unsigned long i = 0; i < delaunay_triangulation.triangles.size();
         i++) {
      triangle_edges = delaunay_triangulation.triangles[i].get_edges();
      for (unsigned int k = 0; k < 3; k++) {
        edge_shared_by_triangles = 0;
        for (unsigned int j = 0; j < delaunay_triangulation.triangles.size();
             j++) {
          if (i != j &&
              delaunay_triangulation.triangles[j].contains(triangle_edges[k])) {
            edge_shared_by_triangles = 1;
            break;
          }
        }
        if (edge_shared_by_triangles) {
          continue;
        } else {
          if (vertices[i].x < 0 || dim_x < vertices[i].x || vertices[i].y < 0 ||
              dim_y < vertices[i].y) {
            continue;
          } else {
            a_b = triangle_edges[k].b.unit_vector_from(triangle_edges[k].a);
            bisector = std::vector<double>{a_b[1], -a_b[0]};
            delta_x = (0 < bisector[0]) ? dim_x - vertices[i].x : vertices[i].x;
            delta_y = (0 < bisector[1]) ? dim_y - vertices[i].y : vertices[i].y;
            scale = (delta_x / abs(bisector[0]) < delta_y / abs(bisector[1]))
                        ? delta_x / abs(bisector[0])
                        : delta_y / abs(bisector[1]);
            Point point =
                Point(abs(std::round(vertices[i].x + bisector[0] * scale)),
                      abs(std::round(vertices[i].y + bisector[1] * scale)));
            edges.push_back(Edge(vertices[i], point));
            vertices.push_back(point);
            edge_vertex_map.push_back(std::make_pair(i, vertices.size() - 1));
            periphery_vertices[i][k] = vertices.size() - 1;
          }
        }
      }
    }
    for (unsigned int i = 0; i < delaunay_triangulation.triangles.size(); i++) {
      for (unsigned int j = 0; j <= i; j++) {
        for (unsigned int k = 0; k < 3; k++) {
          for (unsigned int l = 0; l < 3; l++) {
            if ((i != j || l < k) && periphery_vertices[i][k] != -1 &&
                periphery_vertices[j][l] != -1 &&
                delaunay_triangulation.triangles[i]
                    .get_edges()[k]
                    .shares_point_with(
                        delaunay_triangulation.triangles[j].get_edges()[l])) {
              edges.push_back(Edge(vertices[periphery_vertices[i][k]],
                                   vertices[periphery_vertices[j][l]]));
              edge_vertex_map.push_back(std::make_pair(
                  periphery_vertices[i][k], periphery_vertices[j][l]));
            }
          }
        }
      }
    }
    for (unsigned int i = 0; i < vertices.size(); i++) {
      if (vertices[i].x < 0 || dim_x < vertices[i].x || vertices[i].y < 0 ||
          dim_y < vertices[i].y) {
        vertices.erase(vertices.begin() + i);
        for (unsigned int j = 0; j < edge_vertex_map.size(); j++) {
          if (edge_vertex_map[j].first > i) {
            edge_vertex_map[j].first--;
          }
          if (edge_vertex_map[j].second > i) {
            edge_vertex_map[j].second--;
          }
        }
        i--;
      }
    }
    compute_adj();
    compute_polygons();
    compute_vertex_angles();
  }

  void compute_adj() {
    adj = std::vector<std::vector<bool>>(vertices.size(),
                                         std::vector<bool>(vertices.size(), 0));
    for (unsigned int i = 0; i < edge_vertex_map.size(); i++) {
      adj[edge_vertex_map[i].first][edge_vertex_map[i].second] = 1;
      adj[edge_vertex_map[i].second][edge_vertex_map[i].first] = 1;
    }
  }

  void compute_polygons() {
    for (unsigned int i = 0; i < vertices.size(); i++) {
      find_cycle(i, i, i, std::vector<unsigned long>{});
    }

    unsigned int k = 0;
    for (unsigned long i = 0; i < polygons.size(); i++) {
      for (unsigned long j = 0; j < i; j++) {
        if (polygons[j].vertices.size() == polygons[i].vertices.size()) {
          for (k = 0; k < polygons[i].vertices.size(); k++) {
            if (std::find(
                    polygons[j].vertices.begin(), polygons[j].vertices.end(),
                    polygons[i].vertices[k]) == polygons[j].vertices.end()) {
              break;
            }
          }
          if (k == polygons[i].vertices.size()) {
            polygons.erase(polygons.begin() + j);
            polygon_vertex_map.erase(polygon_vertex_map.begin() + j);
            --i;
            --j;
          }
        }
      }
    }
    unsigned long max_index = 0;
    double temp_area;
    double max_area = 0.0;
    if (polygons.size() > 1) {
      for (unsigned long i = 0; i < polygons.size(); i++) {
        temp_area = polygons[i].area();
        if (temp_area > max_area) {
          max_index = i;
          max_area = temp_area;
        }
      }
      double sum_area = 0.0;
      for (unsigned long i = 0; i < polygons.size(); i++) {
        if (i != max_index) {
          sum_area += polygons[i].area();
        }
      }
      if (abs(sum_area - polygons[max_index].area()) < EPSILON) {
        polygons.erase(polygons.begin() + max_index);
        polygon_vertex_map.erase(polygon_vertex_map.begin() + max_index);
      }
    }
  }

  void compute_vertex_angles() {
    double det, dot, angle, rad;
    std::vector<double> vec_a, vec_b;
    std::vector<std::vector<std::pair<unsigned int, double>>> adj_vertices(
        vertices.size(), std::vector<std::pair<unsigned int, double>>{});
    for (unsigned int i = 0; i < vertices.size(); i++) {
      for (unsigned int j = 0; j < vertices.size(); j++) {
        if (adj[i][j]) {
          rad = -atan2(vertices[j].x - vertices[i].x,
                       -(vertices[j].y - vertices[i].y));
          rad = (rad < 0) ? rad + 2 * PI : rad;
          adj_vertices[i].push_back(std::make_pair(j, rad));
        }
      }
      std::sort(adj_vertices[i].begin(), adj_vertices[i].end(),
                [](std::pair<unsigned int, double> a,
                   std::pair<unsigned int, double> b) {
                  return a.second < b.second;
                });
      vertex_angle_map.push_back(std::vector<double>(1));
      if (adj_vertices[i].size() > 2) {
        vec_a = vertices[adj_vertices[i][0].first].vector_from(vertices[i]);
        vec_b = vertices[adj_vertices[i].back().first].vector_from(vertices[i]);
        dot = vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1];
        det = vec_a[0] * vec_b[1] - vec_b[0] * vec_a[1];
        angle = atan2(det, dot);
        angle = (angle < 0) ? angle + 2 * PI : angle;
        vertex_angle_map[i][0] = angle;
        for (int k = adj_vertices[i].size() - 1; k > 0; k--) {
          vec_a = vertices[adj_vertices[i][k].first].vector_from(vertices[i]);
          vec_b =
              vertices[adj_vertices[i][k - 1].first].vector_from(vertices[i]);
          dot = vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1];
          det = vec_a[0] * vec_b[1] - vec_b[0] * vec_a[1];
          angle = atan2(det, dot);
          angle = (angle < 0) ? angle + 2 * PI : angle;
          vertex_angle_map[i].push_back(angle);
        }
      } else if (adj_vertices[i].size() == 2) {
        vec_a = vertices[adj_vertices[i][0].first].vector_from(vertices[i]);
        vec_b = vertices[adj_vertices[i][1].first].vector_from(vertices[i]);
        dot = vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1];
        det = vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0];
        angle = atan2(det, dot);
        angle = (angle < 0) ? angle + 2 * PI : angle;
        vertex_angle_map[i][0] = angle;
        vertex_angle_map[i].push_back(2 * PI - angle);
      } else if (adj_vertices[i].size() == 1) {
        vertex_angle_map[i][0] = 2 * PI;
      } else {
        continue;
      }
    }
  }

  VoronoiTessellation(Generators generators = Generators(), int dim_x = 1000,
                      int dim_y = 1000)
      : dim_x(dim_x), dim_y(dim_y), generators(generators) {
    DelaunayTriangulation delaunay_triangulation(generators);
    voronoi_from_delaunay(delaunay_triangulation);
  }

  VoronoiTessellation(
      DelaunayTriangulation delaunay_triangulation = DelaunayTriangulation(),
      int dim_x = 1000, int dim_y = 1000)
      : dim_x(dim_x), dim_y(dim_y),
        generators(delaunay_triangulation.generators) {
    voronoi_from_delaunay(delaunay_triangulation);
  }

  VoronoiTessellation(const std::string &fileName) {
    boost::property_tree::ptree tree;
    boost::property_tree::read_xml(fileName, tree);
    std::string key;
    std::vector<Point> polygon_vertices;
    std::vector<unsigned long> polygon_indices;
    unsigned long target, source, vertex_id;
    double x = 0, y = 0;
    bool x_set, y_set;
    unsigned int node_id;
    std::map<unsigned int, unsigned int> xml_node_map;
    std::map<std::string, std::string> key_map;
    std::string id, value;

    const boost::property_tree::ptree &nodes = tree.get_child("graphml");
    for (const boost::property_tree::ptree::value_type &node : nodes) {
      if (node.first == "key") {
        const boost::property_tree::ptree &key_attributes =
            node.second.get_child("<xmlattr>");
        for (const boost::property_tree::ptree::value_type &key_attribute :
             key_attributes) {
          if (key_attribute.first == "attr.name") {
            value = key_attribute.second.data();
          } else if (key_attribute.first == "id") {
            id = key_attribute.second.data();
          }
        }
        key_map[id] = value;
      }
    }
    dim_x = stoi(tree.get("graphml.graph.<xmlattr>.dim_x", "1000"));
    dim_y = stoi(tree.get("graphml.graph.<xmlattr>.dim_y", "1000"));
    const boost::property_tree::ptree &graph_nodes =
        tree.get_child("graphml.graph");
    for (const boost::property_tree::ptree::value_type &graph_node :
         graph_nodes) {
      if (graph_node.first == "node") {
        node_id = std::stoi(graph_node.second.get("<xmlattr>.id", ""));
        const boost::property_tree::ptree &node_entries = graph_node.second;
        x_set = 0;
        y_set = 0;
        for (const boost::property_tree::ptree::value_type &node_entry :
             node_entries) {
          if (node_entry.first == "data" &&
              key_map[node_entry.second.get("<xmlattr>.key", "")] == "x") {
            x = std::stod(node_entry.second.data());
            x_set = 1;
          }
          if (node_entry.first == "data" &&
              key_map[node_entry.second.get("<xmlattr>.key", "")] == "y") {
            y = std::stod(node_entry.second.data());
            y_set = 1;
          }
        }
        if (x_set && y_set) {
          vertices.push_back(Point(x, y));
          xml_node_map[node_id] = vertices.size() - 1;
        }
      }
    }
    for (const boost::property_tree::ptree::value_type &graph_node :
         graph_nodes) {
      if (graph_node.first == "edge") {
        source = static_cast<unsigned long>(
            stoi(graph_node.second.get("<xmlattr>.source", "")));
        target = static_cast<unsigned long>(
            stoi(graph_node.second.get("<xmlattr>.target", "")));
        edges.push_back(Edge(vertices[xml_node_map[source]],
                             vertices[xml_node_map[target]]));
        edge_vertex_map.push_back(std::make_pair(source, target));

      } else if (graph_node.first == "face") {
        polygon_vertices = std::vector<Point>{};
        polygon_indices = std::vector<unsigned long>{};
        const boost::property_tree::ptree &face_entries = graph_node.second;
        for (const boost::property_tree::ptree::value_type &face_entry :
             face_entries) {
          if (face_entry.first == "node") {
            vertex_id = static_cast<unsigned long>(
                stoi(face_entry.second.get("<xmlattr>.node", "")));
            polygon_vertices.push_back(vertices[xml_node_map[vertex_id]]);
            polygon_indices.push_back(xml_node_map[vertex_id]);
          }
        }
        polygon_vertices.push_back(vertices[polygon_indices[0]]);
        polygons.push_back(Polygon(polygon_vertices));
        polygon_indices.push_back(polygon_indices[0]);
        polygon_vertex_map.push_back(polygon_indices);
      }
    }
    compute_adj();
    compute_polygons();
    compute_vertex_angles();
  }

  Generators get_centroids() {
    std::vector<Point> centroids = std::vector<Point>{};
    for (unsigned int i = 0; i < polygons.size(); i++) {
      Point c = polygons[i].centroid();
      centroids.push_back(Point(c.x, c.y, 1.0));
    }
    return Generators(centroids);
  }

  void find_cycle(unsigned long target_vertex, unsigned long current_vertex,
                  unsigned long previous_vertex,
                  std::vector<unsigned long> current_polygon) {
    // adapted from "Finding faces in a planar embedding of a graph", Schneider
    // et al. (2015)
    current_polygon.push_back(current_vertex);
    if (current_vertex == target_vertex && current_polygon.size() > 1) {
      Polygon polygon = Polygon();
      std::vector<unsigned long> temp_vertex_map = std::vector<unsigned long>{};
      for (unsigned int j = 0; j < current_polygon.size(); j++) {
        polygon.vertices.push_back(vertices[current_polygon[j]]);
        temp_vertex_map.push_back(current_polygon[j]);
      }
      polygons.push_back(polygon);
      polygon_vertex_map.push_back(temp_vertex_map);
    } else {
      if (current_polygon.size() == 1) {
        for (unsigned int i = 0; i < vertices.size(); i++) {
          if (adj[current_vertex][i]) {
            find_cycle(target_vertex, i, current_vertex, current_polygon);
          }
        }
      } else {
        double min_rad = 2 * PI;
        unsigned long next_vertex = current_vertex;
        std::vector<double> vec_prev_curr =
            vertices[previous_vertex].vector_from(vertices[current_vertex]);
        for (unsigned int i = 0; i < vertices.size(); i++) {
          if (adj[current_vertex][i] && i != previous_vertex &&
              std::find(std::next(current_polygon.begin()),
                        current_polygon.end(), i) == current_polygon.end()) {
            std::vector<double> vec_curr_i =
                vertices[i].vector_from(vertices[current_vertex]);
            double dot = vec_prev_curr[0] * vec_curr_i[0] +
                         vec_prev_curr[1] * vec_curr_i[1];
            double det = vec_prev_curr[0] * vec_curr_i[1] -
                         vec_prev_curr[1] * vec_curr_i[0];
            double temp_rad = atan2(det, dot);
            temp_rad = (temp_rad < 0) ? temp_rad + 2 * PI : temp_rad;
            if (temp_rad < min_rad) {
              min_rad = temp_rad;
              next_vertex = i;
            }
          }
        }
        if (next_vertex != current_vertex) {
          find_cycle(target_vertex, next_vertex, current_vertex,
                     current_polygon);
        }
      }
    }
  }

  std::map<int, int> map_closest_vertices(VoronoiTessellation voronoi) {
    unsigned int i_dim = voronoi.vertices.size();
    unsigned int j_dim = vertices.size();

    double min_i = 0, min_j = 0;
    std::vector<std::vector<double>> dis(i_dim, std::vector<double>(j_dim));
    std::map<int, int> closest_vertices;
    for (unsigned int i = 0; i < i_dim; i++) {
      for (unsigned int j = 0; j < j_dim; j++) {
        dis[i][j] = voronoi.vertices[i].distance(vertices[j]);
      }
    }
    while (i_dim > 0 && j_dim > 0) {
      double min = std::numeric_limits<double>::max();
      for (unsigned int i = 0; i < voronoi.vertices.size(); i++) {
        for (unsigned int j = 0; j < vertices.size(); j++) {
          if (dis[i][j] < min) {
            min = dis[i][j];
            min_i = i;
            min_j = j;
          }
        }
      }
      closest_vertices[min_i] = min_j;
      for (unsigned int i = 0; i < voronoi.vertices.size(); i++) {
        dis[i][min_j] = std::numeric_limits<double>::max();
      }
      for (unsigned int j = 0; j < vertices.size(); j++) {
        dis[min_i][j] = std::numeric_limits<double>::max();
      }

      i_dim--;
      j_dim--;
    }
    return closest_vertices;
  }

  double mean_error(VoronoiTessellation voronoi,
                    std::map<int, int> closest_vertices) {
    double me = 0.0;
    for (const auto &pair : closest_vertices) {
      me += voronoi.vertices[pair.first].distance(vertices[pair.second]);
    }
    return me / closest_vertices.size();
  }

  Generators approximate_generators(const std::string &fileNameBase,
                                    double epsilon) {
    std::vector<double> errors;
    std::vector<std::pair<int, double>> accepted_errors;

    std::vector<double> vertex_distance, generator_distance;
    Generators approximation_generators = get_centroids();
    Generators new_approximation_generators;
    VoronoiTessellation approximation_voronoi =
        VoronoiTessellation(approximation_generators);
    VoronoiTessellation new_approximation_voronoi =
        VoronoiTessellation(approximation_voronoi);
    std::map<int, int> vertex_map = map_closest_vertices(approximation_voronoi),
                       new_vertex_map;
    double me = mean_error(approximation_voronoi, vertex_map), new_me;

    errors.push_back(me);
    accepted_errors.push_back(std::make_pair(0, me));

    std::vector<unsigned int> generator_polygon_map(
        approximation_generators.generators.size());
    double min_distance, temp_distance;
    unsigned int generator_index = 0;
    for (unsigned int j = 0; j < approximation_voronoi.polygons.size(); j++) {
      min_distance = std::numeric_limits<double>::max();
      for (unsigned int i = 0; i < approximation_generators.generators.size();
           i++) {
        temp_distance = Point(approximation_voronoi.polygons[j].centroid())
                            .distance(approximation_generators.generators[i]);
        if (temp_distance < min_distance) {
          generator_index = i;
          min_distance = temp_distance;
        }
      }
      generator_polygon_map[generator_index] = j;
    }

    std::vector<double> translation;
    std::vector<double> partial_translation;
    int num_translations, ratio;
    bool regression = 1;
    int t = 0;
    while (regression) {
      regression = 0;
      for (unsigned int i = 0; i < approximation_generators.generators.size();
           i++) {
        t++;
        new_approximation_generators = Generators(approximation_generators);
        ratio = polygons[i].area() /
                approximation_voronoi.polygons[generator_polygon_map[i]].area();
        new_approximation_generators.generators[i].r = std::min(
            new_approximation_generators.generators[i].r * ratio,
            0.1 * sqrt(approximation_voronoi.polygons[generator_polygon_map[i]]
                           .area() /
                       PI));
        new_approximation_generators.generators[i].r =
            std::max(new_approximation_generators.generators[i].r, 1.0);

        new_approximation_voronoi =
            VoronoiTessellation(new_approximation_generators);

        if (new_approximation_voronoi.polygon_vertex_map.size() <=
            generator_polygon_map[i]) {
          continue;
        }

        translation = std::vector<double>{0.0, 0.0};
        num_translations = 0;
        for (unsigned int j = 0;
             j < new_approximation_voronoi
                         .polygon_vertex_map[generator_polygon_map[i]]
                         .size() -
                     1;
             j++) {
          partial_translation =
              vertices[vertex_map[new_approximation_voronoi.polygon_vertex_map
                                      [generator_polygon_map[i]][j]]]
                  .vector_from(
                      new_approximation_voronoi
                          .vertices[new_approximation_voronoi.polygon_vertex_map
                                        [generator_polygon_map[i]][j]]);
          translation[0] += partial_translation[0];
          translation[1] += partial_translation[1];
          num_translations++;
        }
        translation[0] /= num_translations;
        translation[1] /= num_translations;

        new_approximation_generators.generators[i].x += translation[0];
        new_approximation_generators.generators[i].y += translation[1];
        new_approximation_voronoi =
            VoronoiTessellation(new_approximation_generators);

        new_vertex_map = map_closest_vertices(new_approximation_voronoi);
        new_me = mean_error(new_approximation_voronoi, new_vertex_map);
        errors.push_back(new_me);

        if (new_me < me - epsilon) {
          approximation_voronoi = new_approximation_voronoi;
          approximation_generators = new_approximation_generators;
          vertex_map = new_vertex_map;
          me = new_me;
          regression = 1;
          accepted_errors.push_back(std::make_pair(t, me));
        }
      }
    }
    std::ofstream errors_file;
    errors_file.open(fileNameBase + "_errors.csv");
    ;
    for (unsigned int i = 0; i < errors.size(); i++) {
      errors_file << errors[i] << "\n";
    }

    std::ofstream accepted_errors_file;
    accepted_errors_file.open(fileNameBase + "_accepted_errors.csv");
    for (unsigned int i = 0; i < accepted_errors.size(); i++) {
      accepted_errors_file << accepted_errors[i].first << ","
                           << accepted_errors[i].second << "\n";
    }
    errors_file.close();
    accepted_errors_file.close();

    return approximation_generators;
  }
  void save(const std::string &fileName) {
    std::ofstream graphml;
    graphml.open(fileName);
    graphml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    int indentation = 0;
    Point centroid;
    graphml << open_tag(
        "graphml",
        std::map<std::string, std::string>{{"xmlns", XMLNS},
                                           {"xmlns:xsi", XMLNS_XSI},
                                           {"xsi:schemaLocation", XSI_SL}},
        indentation, 0, 1);
    indentation++;
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "x"},
                                           {"attr.type", "double"},
                                           {"for", "node"},
                                           {"id", "d0"}},
        indentation, 1, 1);
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "y"},
                                           {"attr.type", "double"},
                                           {"for", "node"},
                                           {"id", "d1"}},
        indentation, 1, 1);

    unsigned int max_num_angles = 0;
    for (unsigned int i = 0; i < vertices.size(); i++) {
      if (vertex_angle_map[i].size() > max_num_angles) {
        max_num_angles = vertex_angle_map[i].size();
      }
    }
    for (unsigned int i = 0; i < max_num_angles; i++) {
      graphml << open_tag("key",
                          std::map<std::string, std::string>{
                              {"attr.name", "angle" + std::to_string(i)},
                              {"attr.type", "double"},
                              {"for", "node"},
                              {"id", "d2" + std::to_string(i)}},
                          indentation, 1, 1);
    }

    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "length"},
                                           {"attr.type", "double"},
                                           {"for", "edge"},
                                           {"id", "d3"}},
        indentation, 1, 1);
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "area"},
                                           {"attr.type", "double"},
                                           {"for", "face"},
                                           {"id", "d4"}},
        indentation, 1, 1);
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "centroid_x"},
                                           {"attr.type", "double"},
                                           {"for", "face"},
                                           {"id", "d5"}},
        indentation, 1, 1);
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "centroid_y"},
                                           {"attr.type", "double"},
                                           {"for", "face"},
                                           {"id", "d6"}},
        indentation, 1, 1);
    graphml << open_tag(
        "graph",
        std::map<std::string, std::string>{{"edgedefault", "undirected"},
                                           {"dim_x", std::to_string(dim_x)},
                                           {"dim_y", std::to_string(dim_y)}},
        indentation, 0, 1);
    indentation++;

    for (unsigned long i = 0; i < vertices.size(); i++) {
      graphml << open_tag(
          "node", std::map<std::string, std::string>{{"id", std::to_string(i)}},
          indentation, 0, 1);
      indentation++;
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d0"}},
                          indentation, 0, 0);
      graphml << vertices[i].x;
      graphml << close_tag("data", indentation, 1);
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d1"}},
                          indentation, 0, 0);
      graphml << vertices[i].y;
      graphml << close_tag("data", indentation, 1);

      for (unsigned int j = 0; j < vertex_angle_map[i].size(); j++) {
        graphml << open_tag("data",
                            std::map<std::string, std::string>{
                                {"key", "d2" + std::to_string(j)}},
                            indentation, 0, 0);
        graphml << vertex_angle_map[i][j];
        graphml << close_tag("data", indentation, 1);
      }

      indentation--;
      graphml << close_tag("node", indentation);
    }
    for (unsigned long i = 0; i < edges.size(); i++) {
      graphml << open_tag(
          "edge",
          std::map<std::string, std::string>{
              {"source", std::to_string(std::min(edge_vertex_map[i].first,
                                                 edge_vertex_map[i].second))},
              {"target", std::to_string(std::max(edge_vertex_map[i].first,
                                                 edge_vertex_map[i].second))}},
          indentation, 0, 1);
      indentation++;
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d3"}},
                          indentation, 0, 0);
      graphml << edges[i].length();
      graphml << close_tag("data", indentation, 1);
      indentation--;
      graphml << close_tag("edge", indentation);
    }
    for (unsigned long i = 0; i < polygons.size(); i++) {
      centroid = polygons[i].centroid();
      graphml << open_tag(
          "face", std::map<std::string, std::string>{{"id", std::to_string(i)}},
          indentation, 0, 1);
      indentation++;
      for (unsigned long j = 0; j < polygon_vertex_map[i].size() - 1; j++) {
        graphml << open_tag(
            "node",
            std::map<std::string, std::string>{
                {"node", std::to_string(polygon_vertex_map[i][j])}},
            indentation, 1, 1);
      }
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d4"}},
                          indentation, 0, 0);
      graphml << polygons[i].area();
      graphml << close_tag("data", indentation, 1);
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d5"}},
                          indentation, 0, 0);
      graphml << centroid.x;
      graphml << close_tag("data", indentation, 1);
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d6"}},
                          indentation, 0, 0);
      graphml << centroid.y;
      graphml << close_tag("data", indentation, 1);
      indentation--;
      graphml << close_tag("face", indentation);
    }
    indentation--;
    graphml << close_tag("graph", indentation);
    indentation--;
    graphml << close_tag("graphml", indentation);
    graphml.close();
  }
  std::string open_tag(std::string name,
                       std::map<std::string, std::string> attributes,
                       int indentation, bool closed = 0, bool newline = 0) {
    std::stringstream graphml;
    for (int i = 0; i < indentation; i++) {
      graphml << "\t";
    }
    graphml << "<" << name;
    for (auto const &attribute : attributes) {
      graphml << " " << attribute.first << "=\"" << attribute.second << "\"";
    }
    if (closed) {
      graphml << "/";
    }
    graphml << ">";
    if (newline) {
      graphml << "\n";
    }
    return graphml.str();
  }
  std::string close_tag(std::string name, int indentation, bool in_line = 0) {
    std::stringstream graphml;
    if (!in_line) {
      for (int i = 0; i < indentation; i++) {
        graphml << "\t";
      }
    }
    graphml << "</" << name << ">\n";
    return graphml.str();
  }
};

struct Circle {
  std::vector<double> c;
  double r;

  Circle(std::vector<double> c = std::vector<double>{}, double r = 0.0)
      : c(c), r(r) {}

  std::vector<std::vector<double>> intersection(const Circle &circle) const {
    double distance = Point(c).distance(circle.c);
    if (distance <= r + circle.r && distance >= abs(r - circle.r)) {
      std::vector<double> v = Point(circle.c).unit_vector_from(c);
      double x =
          (pow(r, 2) - pow(circle.r, 2) + pow(distance, 2)) / (2 * distance);
      double y = sqrt(pow(r, 2) - pow(x, 2));
      std::vector<double> pos1 = std::vector<double>{
          c[0] + x * v[0] - y * v[1], c[1] + x * v[1] + y * v[0]};
      std::vector<double> pos2 = std::vector<double>{
          c[0] + x * v[0] + y * v[1], c[1] + x * v[1] - y * v[0]};
      if (abs(pos1[0] - pos2[0]) < EPSILON &&
          abs(pos1[1] - pos2[1]) < EPSILON) {
        return std::vector<std::vector<double>>{pos1};
      } else {
        return std::vector<std::vector<double>>{pos1, pos2};
      }
    } else {
      return std::vector<std::vector<double>>{};
    }
  }

  std::vector<Edge> draw_between(std::vector<double> pos_a,
                                 std::vector<double> pos_b,
                                 bool smallest_angle = 0) {
    double a_rad = atan2(pos_a[1] - c[1], pos_a[0] - c[0]);
    double b_rad = atan2(pos_b[1] - c[1], pos_b[0] - c[0]);
    double a_deg = deg(a_rad);
    double b_deg = deg(b_rad);
    double rad, rad_inc;
    int a_deg_int, b_deg_int;

    if (abs(a_deg - b_deg) < 1) {
      return std::vector<Edge>{Edge(Point(c[0] + cos(a_deg * PI / 180) * r,
                                          c[1] + sin(a_deg * PI / 180) * r),
                                    Point(c[0] + cos(b_deg * PI / 180) * r,
                                          c[1] + sin(b_deg * PI / 180) * r))};
    }
    if (smallest_angle &&
        ((a_deg > b_deg && b_deg + 360 - a_deg > a_deg - b_deg) ||
         (a_deg < b_deg && a_deg + 360 - b_deg < b_deg - a_deg))) {
      return draw_between(pos_b, pos_a);
    }
    std::vector<Edge> curve;
    if (a_deg < b_deg) {
      a_deg_int = int(a_deg + 1.0);
      b_deg_int = int(b_deg);
      curve.push_back(Edge(
          Point(c[0] + cos(a_deg * PI / 180) * r,
                c[1] + sin(a_deg * PI / 180) * r),
          Point(c[0] + cos(static_cast<double>(a_deg_int) * PI / 180) * r,
                c[1] + sin(static_cast<double>(a_deg_int) * PI / 180) * r)));
      for (int deg = a_deg_int; deg < b_deg_int; deg++) {
        rad = static_cast<double>(deg) * PI / 180;
        rad_inc = static_cast<double>(deg + 1) * PI / 180;
        curve.push_back(
            Edge(Point(c[0] + cos(rad) * r, c[1] + sin(rad) * r),
                 Point(c[0] + cos(rad_inc) * r, c[1] + sin(rad_inc) * r)));
      }
      curve.push_back(
          Edge(Point(c[0] + cos(static_cast<double>(b_deg_int) * PI / 180) * r,
                     c[1] + sin(static_cast<double>(b_deg_int) * PI / 180) * r),
               Point(c[0] + cos(b_deg * PI / 180) * r,
                     c[1] + sin(b_deg * PI / 180) * r)));
    } else {
      a_deg_int = int(a_deg + 1.0);
      b_deg_int = int(b_deg);
      curve.push_back(Edge(
          Point(c[0] + cos(a_deg * PI / 180) * r,
                c[1] + sin(a_deg * PI / 180) * r),
          Point(c[0] + cos(static_cast<double>(a_deg_int) * PI / 180) * r,
                c[1] + sin(static_cast<double>(a_deg_int) * PI / 180) * r)));
      for (int deg = a_deg_int; deg < 360; deg++) {
        rad = static_cast<double>(deg) * PI / 180;
        rad_inc = static_cast<double>(deg + 1) * PI / 180;
        curve.push_back(
            Edge(Point(c[0] + cos(rad) * r, c[1] + sin(rad) * r),
                 Point(c[0] + cos(rad_inc) * r, c[1] + sin(rad_inc) * r)));
      }
      for (int deg = 0; deg < b_deg_int; deg++) {
        rad = static_cast<double>(deg) * PI / 180;
        rad_inc = static_cast<double>(deg + 1) * PI / 180;
        curve.push_back(
            Edge(Point(c[0] + cos(rad) * r, c[1] + sin(rad) * r),
                 Point(c[0] + cos(rad_inc) * r, c[1] + sin(rad_inc) * r)));
      }
      curve.push_back(
          Edge(Point(c[0] + cos(static_cast<double>(b_deg_int) * PI / 180) * r,
                     c[1] + sin(static_cast<double>(b_deg_int) * PI / 180) * r),
               Point(c[0] + cos(b_deg * PI / 180) * r,
                     c[1] + sin(b_deg * PI / 180) * r)));
    }

    return curve;
  }
};

struct CircularVoronoiTessellation {
  // adapted from "Generalized Voronoi Tessellation as a Model of
  // Two-dimensional Cell Tissue Dynamics", Bock et. al (2009)
  std::vector<Point> vertices;
  std::vector<Edge> edges;
  int dim_x, dim_y;
  std::vector<std::pair<unsigned long, unsigned long>> edge_vertex_map;

  CircularVoronoiTessellation(Generators generators = Generators(),
                              int dim_x = 1000, int dim_y = 1000)
      : dim_x(dim_x), dim_y(dim_y) {
    DelaunayTriangulation delaunay_triangulation(generators);
    voronoi_from_delaunay(delaunay_triangulation);
  }

  CircularVoronoiTessellation(
      DelaunayTriangulation delaunay_triangulation = DelaunayTriangulation(),
      int dim_x = 1000, int dim_y = 1000)
      : dim_x(dim_x), dim_y(dim_y) {
    voronoi_from_delaunay(delaunay_triangulation);
  }

  void voronoi_from_delaunay(DelaunayTriangulation delaunay_triangulation) {
    std::vector<Circle> circles(delaunay_triangulation.edges.size());
    Edge edge;
    for (unsigned long i = 0; i < delaunay_triangulation.triangles.size();
         i++) {
      for (unsigned int k = 0; k < 3; k++) {
        if (delaunay_triangulation
                .edges[delaunay_triangulation.triangle_edge_map[i][k]]
                .b.r >
            delaunay_triangulation
                .edges[delaunay_triangulation.triangle_edge_map[i][k]]
                .a.r) {
          edge = Edge(delaunay_triangulation
                          .edges[delaunay_triangulation.triangle_edge_map[i][k]]
                          .b,
                      delaunay_triangulation
                          .edges[delaunay_triangulation.triangle_edge_map[i][k]]
                          .a);
        } else {
          edge = delaunay_triangulation
                     .edges[delaunay_triangulation.triangle_edge_map[i][k]];
        }
        std::vector<double> ccmp = edge.midpoint();
        double scale = -((pow(edge.a.r, 2) + pow(edge.b.r, 2)) /
                         (pow(edge.a.r, 2) - pow(edge.b.r, 2)));
        std::vector<double> a_ccmp = edge.a.vector_from(ccmp);
        std::vector<double> center = std::vector<double>{
            ccmp[0] + a_ccmp[0] * scale, ccmp[1] + a_ccmp[1] * scale};
        double radius =
            ((edge.a.r * edge.b.r) / (pow(edge.a.r, 2) - pow(edge.b.r, 2))) *
            edge.a.distance(edge.b);
        circles[delaunay_triangulation.triangle_edge_map[i][k]] =
            Circle(center, radius);
      }
    }
    std::vector<std::vector<double>> intersections(
        delaunay_triangulation.triangles.size());
    std::vector<std::vector<double>> temp_01, temp_12, temp;
    for (unsigned long i = 0; i < delaunay_triangulation.triangles.size();
         i++) {
      temp_01 =
          circles[delaunay_triangulation.triangle_edge_map[i][0]].intersection(
              circles[delaunay_triangulation.triangle_edge_map[i][1]]);
      temp_12 =
          circles[delaunay_triangulation.triangle_edge_map[i][1]].intersection(
              circles[delaunay_triangulation.triangle_edge_map[i][2]]);

      temp = std::vector<std::vector<double>>{};
      for (unsigned int j = 0; j < temp_01.size(); j++) {
        for (unsigned int k = 0; k < temp_12.size(); k++) {
          if (abs(temp_01[j][0] - temp_12[k][0]) < EPSILON &&
              abs(temp_01[j][1] - temp_12[k][1]) < EPSILON) {
            temp.push_back(temp_01[j]);
          }
        }
      }
      if (temp.empty()) {
        intersections[i] = std::vector<double>{};
      } else if (temp.size() == 1 ||
                 (temp.size() == 2 &&
                  Point(delaunay_triangulation.triangles[i].circumcenter())
                          .distance(temp[0]) <
                      Point(delaunay_triangulation.triangles[i].circumcenter())
                          .distance(temp[1]))) {
        intersections[i] = temp[0];
      } else {
        intersections[i] = temp[1];
      }
    }
    for (unsigned long i = 0; i < delaunay_triangulation.triangles.size();
         i++) {
      if (intersections[i].empty()) {
        continue;
      }
      for (unsigned int k = 0; k < 3; k++) {
        for (unsigned long j = 0; j < i; j++) {
          if (j != i &&
              delaunay_triangulation.triangles[j].contains(
                  delaunay_triangulation
                      .edges[delaunay_triangulation.triangle_edge_map[i][k]])) {
            if (!intersections[j].empty()) {
              convert_curve(
                  circles[delaunay_triangulation.triangle_edge_map[i][k]]
                      .draw_between(intersections[i], intersections[j], 1));
            }
          }
        }
      }
    }
    double boundary_size = 0.0;
    for (unsigned int i = 0; i < delaunay_triangulation.triangles.size(); i++) {
      if (intersections[i].empty()) {
        continue;
      }
      for (unsigned int k = 0; k < 3; k++) {
        if (!intersections[i].empty()) {
          edge = delaunay_triangulation
                     .edges[delaunay_triangulation.triangle_edge_map[i][k]];
          if (edge.a.distance(intersections[i]) > edge.a.r * boundary_size) {
            boundary_size = edge.a.distance(intersections[i]) /
                            delaunay_triangulation.edges[i].a.r;
          }
          if (edge.b.distance(intersections[i]) > edge.b.r * boundary_size) {
            boundary_size = edge.b.distance(intersections[i]) /
                            delaunay_triangulation.edges[i].b.r;
          }
        }
      }
    }
    std::vector<unsigned int> edges_at_periphery;
    std::vector<unsigned int> triangles_at_periphery;
    bool edge_at_periphery;
    for (unsigned long i = 0; i < delaunay_triangulation.triangles.size();
         i++) {
      for (unsigned int k = 0; k < 3; k++) {
        edge_at_periphery = 1;
        for (unsigned long j = 0; j < delaunay_triangulation.triangles.size();
             j++) {
          if (j != i &&
              delaunay_triangulation.triangles[j].contains(
                  delaunay_triangulation
                      .edges[delaunay_triangulation.triangle_edge_map[i][k]])) {
            edge_at_periphery = 0;
            break;
          }
        }
        if (edge_at_periphery) {
          edges_at_periphery.push_back(
              delaunay_triangulation.triangle_edge_map[i][k]);
          triangles_at_periphery.push_back(i);
        }
      }
    }
    Point a, b;
    std::map<unsigned int, std::vector<std::vector<double>>> outer_segments_map;
    std::vector<std::vector<double>> outer_intersections;
    std::vector<double> outer_intersection;
    for (unsigned int i = 0; i < edges_at_periphery.size(); i++) {
      a = delaunay_triangulation.edges[edges_at_periphery[i]].a;
      b = delaunay_triangulation.edges[edges_at_periphery[i]].b;
      outer_intersections =
          Circle(a.pos(), a.r * boundary_size)
              .intersection(Circle(b.pos(), b.r * boundary_size));

      if (outer_intersections.empty() ||
          intersections[triangles_at_periphery[i]].empty()) {
        continue;
      }
      if (delaunay_triangulation.triangulation_contains(
              outer_intersections[0]) &&
          !delaunay_triangulation.triangulation_contains(
              outer_intersections[1])) {
        outer_intersection = outer_intersections[1];
      } else if (delaunay_triangulation.triangulation_contains(
                     outer_intersections[1]) &&
                 !delaunay_triangulation.triangulation_contains(
                     outer_intersections[0])) {
        outer_intersection = outer_intersections[0];
      } else {
        if (Point(outer_intersections[1])
                .distance(
                    delaunay_triangulation.triangles[triangles_at_periphery[i]]
                        .center()) <
            Point(outer_intersections[0])
                .distance(
                    delaunay_triangulation.triangles[triangles_at_periphery[i]]
                        .center())) {
          outer_intersection = outer_intersections[0];
        } else {
          outer_intersection = outer_intersections[1];
        }
      }
      convert_curve(circles[edges_at_periphery[i]].draw_between(
          outer_intersection, intersections[triangles_at_periphery[i]], 1));
      outer_segments_map
          [delaunay_triangulation.edge_vertex_map[edges_at_periphery[i]].first]
              .push_back(outer_intersection);
      outer_segments_map
          [delaunay_triangulation.edge_vertex_map[edges_at_periphery[i]].second]
              .push_back(outer_intersection);
    }
    Point boundary_vertex;
    unsigned int boundary_index;
    Circle boundary_circle;
    double angle_i, angle_a, angle_b;
    for (unsigned int i = 0; i < edges_at_periphery.size(); i++) {
      for (unsigned int j = 0; j < i; j++) {
        if (delaunay_triangulation.edges[edges_at_periphery[i]]
                .shares_point_with(
                    delaunay_triangulation.edges[edges_at_periphery[j]])) {
          if (delaunay_triangulation.edge_vertex_map[edges_at_periphery[i]]
                      .first ==
                  delaunay_triangulation.edge_vertex_map[edges_at_periphery[j]]
                      .first ||
              delaunay_triangulation.edge_vertex_map[edges_at_periphery[i]]
                      .first ==
                  delaunay_triangulation.edge_vertex_map[edges_at_periphery[j]]
                      .second) {
            boundary_index =
                delaunay_triangulation.edge_vertex_map[edges_at_periphery[i]]
                    .first;
            boundary_vertex = delaunay_triangulation.vertices[boundary_index];
            angle_i = deg(
                atan2(delaunay_triangulation.edges[edges_at_periphery[i]].b.y -
                          boundary_vertex.y,
                      delaunay_triangulation.edges[edges_at_periphery[i]].b.x -
                          boundary_vertex.x));
          } else {
            boundary_index =
                delaunay_triangulation.edge_vertex_map[edges_at_periphery[i]]
                    .second;
            boundary_vertex = delaunay_triangulation.vertices[boundary_index];
            angle_i = deg(
                atan2(delaunay_triangulation.edges[edges_at_periphery[i]].a.y -
                          boundary_vertex.y,
                      delaunay_triangulation.edges[edges_at_periphery[i]].a.x -
                          boundary_vertex.x));
            ;
          }
          if (outer_segments_map[boundary_index].size() < 2) {
            continue;
          }
          angle_a = deg(atan2(
              outer_segments_map[boundary_index][0][1] - boundary_vertex.y,
              outer_segments_map[boundary_index][0][0] - boundary_vertex.x));
          angle_b = deg(atan2(
              outer_segments_map[boundary_index][1][1] - boundary_vertex.y,
              outer_segments_map[boundary_index][1][0] - boundary_vertex.x));
          boundary_circle =
              Circle(boundary_vertex.pos(), boundary_size * boundary_vertex.r);
          if ((angle_a < angle_b && (angle_a < angle_i && angle_i < angle_b)) ||
              (angle_a > angle_b &&
               ((angle_a > angle_i && angle_i < angle_b) ||
                (angle_a < angle_i && angle_i > angle_b)))) {
            convert_curve(boundary_circle.draw_between(
                outer_segments_map[boundary_index][1],
                outer_segments_map[boundary_index][0]));
          } else {
            convert_curve(boundary_circle.draw_between(
                outer_segments_map[boundary_index][0],
                outer_segments_map[boundary_index][1]));
          }
        }
      }
    }
  }

  void convert_curve(std::vector<Edge> curve) {
    std::vector<Point>::iterator pos_first;
    std::vector<Point>::iterator pos_last;
    if (!curve.empty()) {
      if (std::find(vertices.begin(), vertices.end(), curve[0].a) ==
          vertices.end()) {
        vertices.push_back(curve[0].a);
        vertices.push_back(curve[0].b);
        edges.push_back(
            Edge(vertices[vertices.size() - 2], vertices[vertices.size() - 1]));
        edge_vertex_map.push_back(
            std::make_pair(vertices.size() - 2, vertices.size() - 1));
      } else {
        vertices.push_back(curve[0].b);
        pos_first = std::find(vertices.begin(), vertices.end(), curve[0].a);
        edges.push_back(Edge(vertices[pos_first - vertices.begin()],
                             vertices[vertices.size() - 1]));
        edge_vertex_map.push_back(
            std::make_pair(pos_first - vertices.begin(), vertices.size() - 1));
      }
      for (unsigned int i = 1; i < curve.size() - 1; i++) {
        vertices.push_back(curve[i].b);
        edges.push_back(
            Edge(vertices[vertices.size() - 2], vertices[vertices.size() - 1]));
        edge_vertex_map.push_back(
            std::make_pair(vertices.size() - 2, vertices.size() - 1));
      }
      if (std::find(vertices.begin(), vertices.end(),
                    curve[curve.size() - 1].b) == vertices.end()) {
        vertices.push_back(curve[curve.size() - 1].a);
        vertices.push_back(curve[curve.size() - 1].b);
        edges.push_back(
            Edge(vertices[vertices.size() - 2], vertices[vertices.size() - 1]));
        edge_vertex_map.push_back(
            std::make_pair(vertices.size() - 2, vertices.size() - 1));
      } else {
        vertices.push_back(curve[curve.size() - 1].a);
        pos_last = std::find(vertices.begin(), vertices.end(),
                             curve[curve.size() - 1].b);
        edges.push_back(Edge(vertices[vertices.size() - 1],
                             vertices[pos_last - vertices.begin()]));
        edge_vertex_map.push_back(
            std::make_pair(vertices.size() - 1, pos_last - vertices.begin()));
      }
    }
  }
  void save(const std::string &fileName) {
    std::ofstream graphml;
    graphml.open(fileName);
    graphml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    int indentation = 0;
    graphml << open_tag(
        "graphml",
        std::map<std::string, std::string>{{"xmlns", XMLNS},
                                           {"xmlns:xsi", XMLNS_XSI},
                                           {"xsi:schemaLocation", XSI_SL}},
        indentation, 0, 1);
    indentation++;
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "x"},
                                           {"attr.type", "double"},
                                           {"for", "node"},
                                           {"id", "d0"}},
        indentation, 1, 1);
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "y"},
                                           {"attr.type", "double"},
                                           {"for", "node"},
                                           {"id", "d1"}},
        indentation, 1, 1);
    graphml << open_tag(
        "key",
        std::map<std::string, std::string>{{"attr.name", "area"},
                                           {"attr.type", "double"},
                                           {"for", "face"},
                                           {"id", "d2"}},
        indentation, 1, 1);
    graphml << open_tag(
        "graph",
        std::map<std::string, std::string>{{"edgedefault", "undirected"},
                                           {"dim_x", std::to_string(dim_x)},
                                           {"dim_y", std::to_string(dim_y)}},
        indentation, 0, 1);

    indentation++;
    for (unsigned long i = 0; i < vertices.size(); i++) {
      graphml << open_tag(
          "node", std::map<std::string, std::string>{{"id", std::to_string(i)}},
          indentation, 0, 1);
      indentation++;
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d0"}},
                          indentation, 0, 0);
      graphml << vertices[i].x;
      graphml << close_tag("data", indentation, 1);
      graphml << open_tag("data",
                          std::map<std::string, std::string>{{"key", "d1"}},
                          indentation, 0, 0);
      graphml << vertices[i].y;
      graphml << close_tag("data", indentation, 1);
      indentation--;
      graphml << close_tag("node", indentation);
    }
    for (unsigned long i = 0; i < edges.size(); i++) {
      graphml << open_tag(
          "edge",
          std::map<std::string, std::string>{
              {"source", std::to_string(std::min(edge_vertex_map[i].first,
                                                 edge_vertex_map[i].second))},
              {"target", std::to_string(std::max(edge_vertex_map[i].first,
                                                 edge_vertex_map[i].second))}},
          indentation, 1, 1);
    }
    indentation--;
    graphml << close_tag("graph", indentation);
    indentation--;
    graphml << close_tag("graphml", indentation);
    graphml.close();
  }
  std::string open_tag(std::string name,
                       std::map<std::string, std::string> attributes,
                       int indentation, bool closed = 0, bool newline = 0) {
    std::stringstream graphml;
    for (int i = 0; i < indentation; i++) {
      graphml << "\t";
    }
    graphml << "<" << name;
    for (auto const &attribute : attributes) {
      graphml << " " << attribute.first << "=\"" << attribute.second << "\"";
    }
    if (closed) {
      graphml << "/";
    }
    graphml << ">";
    if (newline) {
      graphml << "\n";
    }
    return graphml.str();
  }
  std::string close_tag(std::string name, int indentation, bool in_line = 0) {
    std::stringstream graphml;
    if (!in_line) {
      for (int i = 0; i < indentation; i++) {
        graphml << "\t";
      }
    }
    graphml << "</" << name << ">\n";
    return graphml.str();
  }
};

int main(int argc, char *argv[]) {
  std::string fileName, fileNameBase, fileExtension;
  int dim_x = 0, dim_y = 0;
  double approx_eps = 0.0;
  if (argc > 2) {
    fileName = argv[2];
    boost::filesystem::path pathObj(fileName);
    fileNameBase = pathObj.stem().string();
    fileExtension = pathObj.extension().string();
  }
  if (argc == 5) {
    dim_x = std::stoi(argv[3]);
    dim_y = std::stoi(argv[4]);
  }
  if (argc == 6) {
    approx_eps = std::stod(argv[5]);
  }

  if (fileExtension == ".csv") {
    if (strcmp(argv[1], "--voronoi-tessellation") == 0) {
      VoronoiTessellation(Generators(fileName), dim_x, dim_y)
          .save(fileNameBase + ".graphml");
    } else if (strcmp(argv[1], "--circular-voronoi-tessellation") == 0) {
      CircularVoronoiTessellation(Generators(fileName), dim_x, dim_y)
          .save(fileNameBase + ".graphml");
    } else if (strcmp(argv[1], "--delaunay-triangulation") == 0) {
      DelaunayTriangulation(Generators(fileName), dim_x, dim_y)
          .save(fileNameBase + ".graphml");
    }
  } else if (fileExtension == ".graphml") {
    if (strcmp(argv[1], "--generators") == 0) {
      VoronoiTessellation voronoi_tessellation{fileName};
      voronoi_tessellation.approximate_generators(fileNameBase, approx_eps)
          .save(fileNameBase + ".csv");
      voronoi_tessellation.save(fileNameBase + ".graphml");
    } else if (strcmp(argv[1], "--centroids") == 0) {
      VoronoiTessellation voronoi_tessellation{fileName};
      voronoi_tessellation.get_centroids().save(fileNameBase + ".csv");
      voronoi_tessellation.save(fileNameBase + ".graphml");
    }
  }
  return 0;
}