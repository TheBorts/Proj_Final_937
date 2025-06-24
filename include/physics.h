#ifndef PHYSICS_H
#define PHYSICS_H

#include <math.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "clas.h"
#include "shaders.h"

float deltaTime = 0.0f;
float lastFrame = 0.0f;
glm::vec3 gravity(0.0f, -0.5f, 0.0f); // Gravity std::vector
class AABB {
    public:
        glm::vec3 min;
        glm::vec3 max;
        long long num_faces = 0; // Number of faces contained within this AABB
        std::vector<face> faces; // Faces contained within this AABB
        mesh* m = nullptr; // Pointer to the object this AABB represents
        AABB* children[2] = {nullptr}; // Pointers to child AABBs for octree structure
        bool isLeaf = true; // Indicates if this AABB is a leaf node in the octree

        AABB() : m(nullptr), min(glm::vec3(0.0f)), max(glm::vec3(0.0f)) {
            for (int i = 0; i < 8; ++i) {
                children[i] = nullptr;
            }
            num_faces = 0;
            set_faces();
            isLeaf = true;
        }

        AABB(const std::string& filename) {
            m = new mesh(filename);
            calculate_bounds();
            generate_children();
        }

        ~AABB() {
            for (int i = 0; i < 2; ++i) {
                if (children[i]) {
                    delete children[i]; // Clean up child AABBs
                }
            }
            
            delete m; // Clean up the mesh object               
        }

        void detect_collision(AABB* other, std::vector<face>& collision_faces) {
            if (!other || !m || !other->m) return;

            // Check for overlap in AABBs
            if (min.x > other->max.x || max.x < other->min.x ||
                min.y > other->max.y || max.y < other->min.y ||
                min.z > other->max.z || max.z < other->min.z) {
                return; // No collision
            }

            // If this AABB is a leaf node, check its faces against the other's faces
            if (isLeaf && other->isLeaf) {
                for (const auto& f : faces) {
                    for (const auto& of : other->faces) {
                        
                        collision_faces.push_back(f);
                    }
                }
            } else {
                // If not leaf, check children recursively
                for (int i = 0; i < 2; ++i) {
                    if (children[i]) {
                        other->detect_collision(children[i], collision_faces);
                    }
                }
            }
        }

    bool triangleTriangleIntersectionComplete(const face& tri1, const face& tri2) {

        // Compute plane of triangle 1
        glm::vec3 tr1[3] = {
            glm::vec3(m->mod->vertices[tri1.vertices[0]-1].x, 
                       m->mod->vertices[tri1.vertices[0]-1].y, 
                       m->mod->vertices[tri1.vertices[0]-1].z),
            glm::vec3(m->mod->vertices[tri1.vertices[1]-1].x,
                          m->mod->vertices[tri1.vertices[1]-1].y, 
                          m->mod->vertices[tri1.vertices[1]-1].z),
            glm::vec3(m->mod->vertices[tri1.vertices[2]-1].x,
                          m->mod->vertices[tri1.vertices[2]-1].y, 
                          m->mod->vertices[tri1.vertices[2]-1].z)
        };
        glm::vec3 tr2[3] = {
            glm::vec3(m->mod->vertices[tri2.vertices[0]-1].x, 
                       m->mod->vertices[tri2.vertices[0]-1].y, 
                       m->mod->vertices[tri2.vertices[0]-1].z),
            glm::vec3(m->mod->vertices[tri2.vertices[1]-1].x,
                          m->mod->vertices[tri2.vertices[1]-1].y, 
                          m->mod->vertices[tri2.vertices[1]-1].z),
            glm::vec3(m->mod->vertices[tri2.vertices[2]-1].x,
                          m->mod->vertices[tri2.vertices[2]-1].y, 
                          m->mod->vertices[tri2.vertices[2]-1].z)
        };
        glm::vec3 normal1 = glm::cross(tr1[1] - tr1[0], tr1[2] - tr1[0]);

        float d1 = -glm::dot(normal1, tr1[0]);

        // Distances of triangle 2 vertices to plane 1
        float dist2[3];
        dist2[0] = glm::dot(normal1, tr2[0]) + d1;
        dist2[1] = glm::dot(normal1, tr2[1]) + d1;
        dist2[2] = glm::dot(normal1, tr2[2]) + d1;

        // Check if all distances are same sign
        if ((dist2[0] * dist2[1] > 0.0f) && (dist2[1] * dist2[2] > 0.0f)) {
            return false;
        }

        // Compute plane of triangle 2
        glm::vec3 normal2 = glm::cross(tr2[1] - tr2[0], tr2[2] - tr2[0]);
        float d2 = -glm::dot(normal2, tr2[0]);

        // Distances of triangle 1 vertices to plane 2
        float dist1[3];
        dist1[0] = glm::dot(normal2, tr1[0]) + d2;
        dist1[1] = glm::dot(normal2, tr1[1]) + d2;
        dist1[2] = glm::dot(normal2, tr1[2]) + d2;

        // Check if all distances are same sign
        if ((dist1[0] * dist1[1] > 0.0f) && (dist1[1] * dist1[2] > 0.0f)) {
            return false;
        }

        // Compute direction of intersection line
        glm::vec3 direction = glm::cross(normal1, normal2);

        // Compute and compare intervals
        float max1 = -FLT_MAX, min1 = FLT_MAX;
        float max2 = -FLT_MAX, min2 = FLT_MAX;

        for (int i = 0; i < 3; ++i) {
            float projection = glm::dot(direction, tr1[i]);
            max1 = std::max(max1, projection);
            min1 = std::min(min1, projection);
            
            projection = glm::dot(direction, tr2[i]);
            max2 = std::max(max2, projection);
            min2 = std::min(min2, projection);
        }

        if (max1 < min2 || max2 < min1) {
            return false;
        }

        // More detailed tests for edge intersections
        glm::vec3 edges1[3] = {tr1[1] - tr1[0], tr1[2] - tr1[1], tr1[0] - tr1[2]};
        glm::vec3 edges2[3] = {tr2[1] - tr2[0], tr2[2] - tr2[1], tr2[0] - tr2[2]};

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                glm::vec3 axis = glm::cross(edges1[i], edges2[j]);
                if (glm::length(axis) < 1e-6f) continue; // Skip parallel edges

                // Project both triangles onto the axis
                float minProj1 = FLT_MAX, maxProj1 = -FLT_MAX;
                float minProj2 = FLT_MAX, maxProj2 = -FLT_MAX;

                for (const auto& v : tr1) {
                    float proj = glm::dot(axis, v);
                    minProj1 = std::min(minProj1, proj);
                    maxProj1 = std::max(maxProj1, proj);
                }

                for (const auto& v : tr2) {
                    float proj = glm::dot(axis, v);
                    minProj2 = std::min(minProj2, proj);
                    maxProj2 = std::max(maxProj2, proj);
                }

                if (maxProj1 < minProj2 || maxProj2 < minProj1) {
                    return false;
                }
            }
    }

    return true;
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

        void find_faces(std::vector<face>& faces_pai, long long& num_faces_pai) {
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
        
        std::vector<glm::vec3> position;
        std::vector<glm::vec3> velocity;

        bool isStatic = true; // Indicates if the object is static or dynamic
        std::vector<bool> isVertStatic; // Indicates if the object is static or dynamic

        float mass = 1.0f; // Default mass for the object
        float springRestitution = 0.5f; // Default restitution for the object
        float friction = 0.5f; // Default friction for the object
        float linearDamping = 0.005f; // Default linear damping for the object
        float angularDamping = 0.1f; // Default angular damping for the object
        float dragCoefficient = 0.47f; // Default drag coefficient for the object
        float restitution = 0.5f; // Default restitution for the object


        object(const std::string& filename) {
            tree = new AABB(filename);
            for (long long i = 0; i < tree->m->mod->num_vertices; i++) {
                position.push_back(glm::vec3(tree->m->mod->vertices[i].x, 
                                             tree->m->mod->vertices[i].y, 
                                             tree->m->mod->vertices[i].z));
            }
            velocity = std::vector<glm::vec3>(position.size(), glm::vec3(0.0f));
            isVertStatic = std::vector<bool>(position.size(), true); // Initialize all vertices as static
            
        }

        ~object() {
            if (tree) {
                delete tree; // Clean up the AABB tree
            }
            
        }

        void set_positions(){
            for (long long i = 0; i < tree->m->mod->num_vertices; i++) {
                position[i] = glm::vec3(tree->m->mod->vertices[i].x, 
                                        tree->m->mod->vertices[i].y, 
                                        tree->m->mod->vertices[i].z);
            }
        }

        void calculate_forces() {
            if (isStatic) return; // Skip force calculations for static objects
            // using every edge as a spring
            for (size_t i = 0; i < tree->m->mod->num_vertices; i++) {
                long long v1_index = i;
                for (size_t j = i+1; j < tree->m->mod->num_vertices; j++) {
                    long long v2_index = j;

                    glm::vec3 v1 = position[v1_index];
                    glm::vec3 v2 = position[v2_index];

                    glm::vec3 initial_v1 = glm::vec3(tree->m->mod->vertices[v1_index].x, 
                                                      tree->m->mod->vertices[v1_index].y, 
                                                      tree->m->mod->vertices[v1_index].z);
                    glm::vec3 initial_v2 = glm::vec3(tree->m->mod->vertices[v2_index].x, 
                                                      tree->m->mod->vertices[v2_index].y, 
                                                      tree->m->mod->vertices[v2_index].z);
                    float restLength = glm::length(initial_v1 - initial_v2);
                    glm::vec3 currentLength = v1 - v2;
                    glm::vec3 force = - springRestitution * (glm::length(currentLength) - restLength) * glm::normalize(currentLength);

                    // Apply the force to both vertices
                    velocity[v1_index] += force * deltaTime / mass;
                    velocity[v2_index] -= force * deltaTime / mass; // Opposite direction for the second vertex
                    }                       
                // Apply gravity
                velocity[v1_index] += gravity * deltaTime;
                                
                // Apply damping
                velocity[v1_index] *= (1.0f - linearDamping);

                if (position[v1_index].y < -1.0f) { // Simple ground collision
                    position[v1_index].y = -1.0f; // Reset position to ground level
                    velocity[v1_index].y = -velocity[v1_index].y * restitution; // Reverse velocity with restitution
                }
            }
        }

        void update_positions() {
            if (isStatic) return; // Skip position updates for static objects
            for (size_t i = 0; i < position.size(); i++) {
                if (isVertStatic[i]) continue; // Skip static vertices
                position[i] += velocity[i] * deltaTime;
            }
        }

        void draw(Shader shader) {
            
            tree->m->update_data(position);
            tree->m->draw_mesh(shader);
        }
};       

class world {
    public:
        object** objs;
        int num_obj = 0;
        world(){
            objs = nullptr;
        }

        ~world() {
            for (int i = 0; i < num_obj; i++) {
                object* obj = objs[i];
                if (obj) {
                    delete obj; // Clean up each object
                }
            }
            delete[] objs; // Clean up the array of objects
        }

        void add_object(const std::string& filename) {
            object** new_obj = new object*[num_obj + 1];
            for (int i = 0; i < num_obj; i++) {
                new_obj[i] = objs[i];
            }
            new_obj[num_obj] = new object(filename);
            delete[] objs; // Free the old array
            objs = new_obj;
            num_obj++;
        }

        void start_world() {
            for (int i = 0; i < num_obj; i++) {
                object* obj = objs[i];
                obj->tree->m->setup_mesh(i);
            }
        }

        void draw_world(Shader shader) {
            update_world(); // Update positions and forces before drawing

            for (int i = 0; i < num_obj; i++) {
                objs[i]->draw(shader);
                
            }

        }

        void update_world() {
            for (int i = 0; i < num_obj; i++) {
                objs[i]->calculate_forces();
                objs[i]->update_positions();
            }
            //detect_collisions(); // Detect collisions after updating positions
        }

        void detect_collisions() {
            for (int i = 0; i < num_obj; i++) {
                object* obj1 = objs[i];
                for (int j = i + 1; j < num_obj; j++) {
                    object* obj2 = objs[j];
                    std::vector<face> collision_faces = std::vector<face>();
                    obj1->tree->detect_collision(obj2->tree, collision_faces);
                    if (!collision_faces.empty()) {
                        
                        std::cout << "Collision detected between object " << i << " and object " << j << std::endl;
                    }
                }
            }
        }



};

class writer {
    public:
        std::string objPath;
        std::string mtlPath;
        std::ofstream fileMTL;
        std::ofstream fileOBJ;
        
        writer(const std::string& filenameOBJ){
            objPath = filenameOBJ;
            fileOBJ.open(filenameOBJ);
            if (!fileOBJ.is_open()) {
                std::cerr << "Error opening file for writing: " << filenameOBJ << std::endl;
                return;
            }
        }

        writer(const std::string& filenameOBJ, const std::string& filenameMTL) {
            objPath = filenameOBJ;
            mtlPath = filenameMTL;
        }

        ~writer() {
            if (fileOBJ.is_open()) {
                fileOBJ.close();
            }
            if (fileMTL.is_open()) {
                fileMTL.close();
            }
        }

        void setPaths(const std::string& filenameOBJ, const std::string& filenameMTL) {
            objPath = filenameOBJ;
            mtlPath = filenameMTL;
        }

        void writeMTL(const world& myWorld) {
            fileMTL.open(mtlPath);
            if (!fileMTL.is_open()) {
                std::cerr << "File not open for writing: " << mtlPath << std::endl;
                return;
            }

            
            fileMTL << "# MTL fileMTL generated by writer\n";
            for (int i = 0; i < myWorld.num_obj; i++) {
                object* obj = myWorld.objs[i];
                fileMTL << "newmtl " << obj->tree->m->mod->mat.name << "\n";
                fileMTL << "Ka " << obj->tree->m->mod->mat.ambient[0] << " "
                << obj->tree->m->mod->mat.ambient[1] << " "
                << obj->tree->m->mod->mat.ambient[2] << "\n";
                fileMTL << "Kd " << obj->tree->m->mod->mat.diffuse[0] << " "
                << obj->tree->m->mod->mat.diffuse[1] << " "
                << obj->tree->m->mod->mat.diffuse[2] << "\n";
                fileMTL << "Ks " << obj->tree->m->mod->mat.specular[0] << " "
                << obj->tree->m->mod->mat.specular[1] << " "
                << obj->tree->m->mod->mat.specular[2] << "\n";
                fileMTL << "Ns " << obj->tree->m->mod->mat.shininess << "\n";
                if (!obj->tree->m->mod->path_texture.empty()) {
                    fileMTL << "map_Kd " << obj->tree->m->mod->path_texture << "\n";
                }
                fileMTL << "\n"; // Add a newline for readability
            }
            fileMTL.flush();
            if (fileMTL.fail()) {
                std::cerr << "Error writing to file: " << mtlPath << std::endl;
            }

            std::cout << "Writing MTL file: " << mtlPath << std::endl;          
            
            fileMTL.close();

            if (std::filesystem::file_size(mtlPath) == 0) {
                std::cerr << "Warning: File was created but is empty!" << std::endl;
            } else {
                std::cout << "Successfully wrote MTL file with " 
                        << std::filesystem::file_size(mtlPath) 
                        << " bytes to " << mtlPath << std::endl;
            }

        }

        void writeOBJ(const world& myWorld) {
            fileOBJ.open(objPath);
            if (!fileOBJ.is_open()) {
                std::cerr << "FileOBJ not open for writing: " << objPath << std::endl;
                return;
            }

            fileOBJ << "mtllib " << mtlPath << "\n";
            for (int i = 0; i < myWorld.num_obj; i++) {
                object* obj = myWorld.objs[i];
                fileOBJ << "o " << obj->tree->m->mod->mod_name << "\n";
                fileOBJ << "usemtl " << obj->tree->m->mod->mat.name << "\n";
                for (size_t j = 0; j < obj->position.size(); j++) {
                    fileOBJ << "v " << obj->position[j].x << " " << obj->position[j].y << " " << obj->position[j].z << "\n";
                }
                for (size_t j = 0; j < obj->tree->m->mod->num_vertices; j++) {
                    fileOBJ << "vn " << obj->tree->m->proper_normals[j].x << " "
                         << obj->tree->m->proper_normals[j].y << " "
                         << obj->tree->m->proper_normals[j].z << "\n";
                }
                for (int j = 0; j < obj->tree->m->mod->num_textures; j++) {
                    fileOBJ << "vt " << obj->tree->m->mod->textures[j].u << " "
                         << obj->tree->m->mod->textures[j].v << "\n";
                }
                for (long long j = 0; j < obj->tree->m->mod->num_faces; j++) {
                    face& f = obj->tree->m->mod->faces[j];
                    fileOBJ << "f";
                    for (long long k = 0; k < f.num_vertex; k++) {
                        fileOBJ << " ";

                        if (obj->tree->m->mod->num_textures == 0) {
                            fileOBJ << " " << f.vertices[k] << "//" << f.normals[k];
                        } else {
                            fileOBJ << " " << f.vertices[k] << "/" << f.textures[k] << "/" << f.normals[k];
                        }
                    }
                    fileOBJ << "\n";
                }
            }
            fileOBJ.close();
        }
};

#endif // PHYSICS_H