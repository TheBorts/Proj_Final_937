#include <math.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "clas.h"

class AABB {
    public:
        glm::vec3 min;
        glm::vec3 max;
        long long num_faces = 0; // Number of faces contained within this AABB
        vector<face> faces; // Faces contained within this AABB
        mesh* m = nullptr; // Pointer to the object this AABB represents
        AABB* children[2] = {nullptr}; // Pointers to child AABBs for octree structure
        bool isLeaf = true; // Indicates if this AABB is a leaf node in the octree

        AABB() : m(nullptr), min(glm::vec3(0.0f)), max(glm::vec3(0.0f)) {
            for (int i = 0; i < 8; ++i) {
                children[i] = nullptr;
            }
        }

        AABB(const std::string& filename) {
            m = new mesh(filename);
            calculate_bounds();
            generate_children();
        }



    private:

        void set_faces(){
            for (int i = 0; i < m->mod->num_faces; i++) {
                face f = m->mod->faces[i];
                faces.push_back(f);
            }
        }

        bool contains(const glm::vec3& point) const {
            return (point.x >= min.x && point.x <= max.x) &&
                (point.y >= min.y && point.y <= max.y) &&
                (point.z >= min.z && point.z <= max.z);
        }

        void calculate_bounds() {
            if (!m || m->data.empty()) return;

            min = glm::vec3(INFINITY);
            max = glm::vec3(-INFINITY);

            for (const auto& vertex : m->data) {
                min = glm::min(min, vertex.position);
                max = glm::max(max, vertex.position);
            }
        }

        void find_faces(vector<face>& faces_pai, long long& num_faces_pai) {
            if (!m || m->data.empty()) return;
            
            for (int i = 0; i < num_faces_pai; i++) {
                face f = faces_pai[i];
                bool faceContained = false;
                for (int j = 0; j < f.num_vertex; j++) {
                    glm::vec3 vertexPos = glm::vec3(m->mod->vertices[f.vertices[j] - 1].x,
                                                    m->mod->vertices[f.vertices[j] - 1].y,
                                                    m->mod->vertices[f.vertices[j] - 1].z);
                    if (contains(vertexPos)) {
                        faceContained = true;
                        break;
                    }
                }
                if (faceContained) {
                    faces.push_back(f);
                    num_faces++;
                }
            }
           
        }

        void generate_children() {
            if (!m || m->data.empty()) return;

            if (num_faces <= 10) {
                // If there's only one face or less, no need to split further
                isLeaf = true;
                return;
            }

            glm::vec3 center = (min + max) * 0.5f;

            // Create 2 children AABBs on the dominant axis
            glm::vec3 size = max - min;
            int dominantAxis = (size.x > size.y && size.x > size.z) ? 0 :
                               (size.y > size.x && size.y > size.z) ? 1 :
                               2;
            glm::vec3 childMin1 = min;
            glm::vec3 childMax1 = max;
            glm::vec3 childMin2 = min;
            glm::vec3 childMax2 = max;
            if (dominantAxis == 0) { // Split along x-axis
                childMax1.x = center.x;
                childMin2.x = center.x;
            } else if (dominantAxis == 1) { // Split along y-axis
                childMax1.y = center.y;
                childMin2.y = center.y;
            } else { // Split along z-axis
                childMax1.z = center.z;
                childMin2.z = center.z;
            }
            children[0] = new AABB();
            children[0]->min = childMin1;
            children[0]->max = childMax1;
            children[0]->m = m; // Assign the same mesh to the first child
            children[0]->find_faces(faces, num_faces);

            children[1] = new AABB();
            children[1]->min = childMin2;
            children[1]->max = childMax2;
            children[1]->m = m; // Assign the same mesh to the second child
            children[1]->find_faces(faces, num_faces);

            faces.clear(); // Clear the faces in the parent AABB after distributing them to children
            num_faces = 0; // Reset the face count in the parent AABB
            isLeaf = false; // This AABB is no longer a leaf node
            children[0]->generate_children(); // Recursively generate children for the first child
            children[1]->generate_children(); // Recursively generate children for the second child
        }
};

class object{
    public:
        AABB* tree = nullptr;
        glm::vec3 position;
        glm::vec3 rotation;
        glm::vec3 scale;
        object() : tree(nullptr), position(glm::vec3(0.0f)), rotation(glm::vec3(0.0f)), scale(glm::vec3(1.0f)) {}
        object(const std::string& filename) {
            tree = new AABB(filename);
            position = glm::vec3(0.0f);
            rotation = glm::vec3(0.0f);
            scale = glm::vec3(1.0f);
        }

        void 



};       

class world {
    public:
        std::vector<mesh> meshes;
        std::vector<lightSource> lights;
        glm::vec3 ambientLight;

        world() : ambientLight(glm::vec3(0.1f, 0.1f, 0.1f)) {}

        void add_mesh(const mesh& m) {
            meshes.push_back(m);
        }

        void add_light(const lightSource& light) {
            lights.push_back(light);
        }

        void set_ambient_light(const glm::vec3& color) {
            ambientLight = color;
        }

        void render(Shader& shader) {
            shader.use();
            shader.setVec3("ambientLight", ambientLight);

            for (const auto& light : lights) {
                shader.setVec3("light.position", light.position[0], light.position[1], light.position[2]);
                shader.setVec3("light.ambient", light.ambient[0], light.ambient[1], light.ambient[2]);
                shader.setVec3("light.diffuse", light.diffuse[0], light.diffuse[1], light.diffuse[2]);
                shader.setVec3("light.specular", light.specular[0], light.specular[1], light.specular[2]);
            }

            for (auto& m : meshes) {
                m.draw_mesh(shader);
            }
        }
};