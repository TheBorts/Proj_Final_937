#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <unordered_map>
#include <sstream>
#include <iomanip>
#include <random>
#include <cmath>

using namespace std;
using namespace glm;

struct Material {
    string name;
    vec3 ambient = vec3(0.0f);
    vec3 diffuse = vec3(0.0f);
    vec3 specular = vec3(0.0f);
    float shineness = 1.0f;
};

struct Vertex {
    int index;
    vec3 position;
    vec3 normal;
    vec3 color = vec3(0.0f);
    Material material;
    // ADS calculado por vértice -> Gourard Shading
};

struct Triangle {
    Vertex v0, v1, v2;
    vec3 color = vec3(0.0f); // RGB
    Material material;
};


vector<Vertex> vertices;
vector<Triangle> triangles;
vec3 light_pos = vec3(0.0f, 5.0f, 5.0f);
vec3 view_pos = vec3(0.0f);
vec3 light_ambient = vec3(1.0);
vec3 light_diffuse = vec3(1.0);
vec3 light_specular = vec3(1.0);

Material getMTL(string path, string material) {
    ifstream file(path);
    string line;
    Material m;

    while (getline(file, line)) {
        if (line.substr(0, 5) == "newmtl" && line.substr(7) == material) {
            m.name = line.substr(7);
            getline(file, line);
            vec3 ret;
            ret.r = stof(line.substr(3, 5));   // R
            ret.g = stof(line.substr(7, 9));   // G
            ret.b = stof(line.substr(11, 13)); // B
            if (line.substr(0, 1) == "Ka") { // cor ambiente
                m.ambient = ret;
            } else if (line.substr(0, 1) == "Kd") {
                m.diffuse = ret;
            } else if (line.substr(0, 1) == "Ks") {
                m.specular = ret;
            }
        }
    }
    return m;
}

Triangle setTriangle(Vertex v0, Vertex v1, Vertex v2, Material m) {
    Triangle t;
    t.material = m;

    t.v0 = v0;
    v0.material = m;

    t.v1 = v1;
    v1.material = m;

    t.v2 = v2;
    v2.material = m;

    return t;
}

void loadOBJ(const string& path) {
    ifstream file(path);
    string line;
    vector<vec3> temp_normais;
    Material current_material;
    string mtl;
    int index = 1;
    while (getline(file, line)) {
        if (line.substr(0, 5) == "mtllib") {
            mtl = line.substr(7);
        } else if (line[0] == 'v') {
            vec3 v;
            v.x = (int) line[2];
            v.y = (int) line[4];
            v.z = (int) line[6];
            Vertex ver;
            ver.index = index;
            ver.position = v;
            ver.normal = normalize(v);
            vertices.push_back(ver);
            ++index;
        } else if (line.substr(0, 5) == "usemtl") {
            current_material = getMTL(mtl, line.substr(7));
        }
        else if (line[0] == 'f') {
            int a = (int) line[2];
            int b = (int) line[4];
            int c = (int) line[6];

            Triangle t = setTriangle(vertices.at(a-1), vertices.at(b-1),
            vertices.at(c-1), current_material);

            triangles.push_back(t);
        }
    }
}

vector<vec3> generateSphericalRays(float theta_begin, float theta_end, float phi_begin, float phi_end, int num_rays) {
    vector<vec3> rays;
    float phi_i = (phi_begin + phi_end) / num_rays;
    for (int i = phi_begin; i < phi_end; i += phi_i) {
        float theta_i = (theta_begin + theta_end) / num_rays;
        for (int j = theta_begin; j < theta_end; j += theta_i) {
            float x = sin(i) * cos(j);
            float y = sin(i) * sin(j);
            float z = cos(i);
            vec3 direction = vec3(x, y, z);
            direction = direction / normalize(direction);
            rays.push_back(direction);
        }
    }

    return rays;
}

bool ray_triangle_intersection(vec3 orig, vec3 dir, Triangle t, vec3& intersection) {
    vec3 e1 = t.v1.position - t.v0.position;
    vec3 e2 = t.v2.position - t.v0.position;
    
    vec3 p = cross(dir, e2);
    float det = dot(e1, p);

    if (abs(det) < epsilon<float>())
        return false;

    vec3 s = orig - t.v0.position;
    float inv_det = 1.0 / det;

    float u = inv_det * dot(s, p);
    if (u < 0.0 || u > 1.0)
        return false;

    vec3 q = cross(s, e1);

    float v = inv_det * dot(dir, q);
    if (v < 0.0 || u + v > 1.0)
        return false;

    float td = inv_det * dot(e2, q);
    if (td >= epsilon<float>()) {
        intersection = orig + td * dir;
        return true;
    }

    return false;
}

void castRaysAndDetectHits(vector<vec3>& ray_directions, vector<vec3>& intersections, vector<Triangle>& hit_faces) {
    for (auto& rd : ray_directions) {
        for (int i = 0; i < triangles.size(); ++i) {
            Triangle t = triangles[i];
            vec3 v0 = t.v0.position;
            vec3 v1 = t.v1.position;           
            vec3 v2 = t.v2.position;

            // descarta faces opostas
            vec3 normal = normalize(cross(v1 - v0, v2 - v0));
            if (dot(normal, rd) >= 0)
                continue;

            // se tiver intersecção com a face, adiciona ao vetor
            vec3 intersection;
            if (ray_triangle_intersection(light_pos, rd, t, intersection)) {
                intersections.push_back(intersection);
                hit_faces.push_back(t);
            }
        }
    }
}

vec3 calculateADS(Vertex v) {
    vec3 N = v.normal;
    vec3 L = normalize(light_pos - v.position);
    vec3 V = normalize(view_pos - v.position);
    vec3 R = reflect(-L, N);
    vec3 a = light_ambient * v.material.ambient;
    vec3 d = light_diffuse * v.material.diffuse * glm::max(dot(N, L), 0.0f);
    vec3 s = light_specular * v.material.specular * pow(glm::max(dot(R, V), 0.0f), v.material.shineness);
    return clamp(a + d + s, 0.0f, 1.0f);
}

void update_colors(vector<Triangle>& hit_faces) {
    // atualizamos as cores de todos os vértices pertencentes à triângulos atingidos
    for (int i = 0; i < hit_faces.size(); ++i) {
        Triangle& t = hit_faces[i];
        t.v0.color = calculateADS(t.v0);
        t.v1.color = calculateADS(t.v1);
        t.v2.color = calculateADS(t.v1);
        t.color = (t.v0.color + t.v1.color + t.v2.color) / 3.0f;
    }
}

void writeOBJandMTL(const string& objPath, const string& mtlPath) {
    ofstream obj(objPath);
    ofstream mtl(mtlPath);
    // referenciamos no .obj onde está o .mtl
    obj << "mtllib " << mtlPath << "\n";
    // escrevemos todos os vértices e suas posições, separados por espaço (ex. "v 0 0 1\n")
    for (auto& v : vertices)
        obj << "v " << v.position.x << " " << v.position.y << " " << v.position.z << "\n";

    // escrevemos todas as normais dos vértices e suas posições, separados por espaço
    for (auto& v : vertices)
        obj << "vn " << v.normal.x << " " << v.normal.y << " " << v.normal.z << "\n";

    // escrevemos sobre os materiais e os triângulos
    for (int i = 0; i < triangles.size(); ++i) {
        auto& t = triangles[i];
        mtl << "newmtl " << t.material.name << "\n";
        mtl << "Ka " << t.material.ambient.r << " " << t.material.ambient.g << " " << t.material.ambient.b << "\n";
        mtl << "Kd " << t.material.diffuse.r << " " << t.material.diffuse.g << " " << t.material.diffuse.b << "\n";
        mtl << "Ks " << t.material.specular.r << " " << t.material.specular.g << " " << t.material.specular.b << "\n\n";

        obj << "usemtl " << t.material.name << "\n";
        obj << "f " << t.v0.index << " " << t.v1.index << " " << t.v2.index << "\n";
    }
}


int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Uso: ./Lab3 mesh.obj\n";
        return 1;
    }

    loadOBJ(argv[1]);

    vector<vec3> ray_directions = generateSphericalRays(0, 2 * pi<float>(), 0, pi<float>(), 50);
    vector<vec3> intersections;
    vector<Triangle> hit_faces;
    castRaysAndDetectHits(ray_directions, intersections, hit_faces);
    update_colors(hit_faces);

    writeOBJandMTL("output.obj", "output.mtl");

    return 0;
}
