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

// Class that represents a 3D object in the physics simulation
class object{
    public:
        // Pointer to the mesh with the 3D model data and the connection to openGl
        mesh* m = nullptr;
        
        // Vectors to store the positions and velocities of the vertices
        std::vector<glm::vec3> position;
        std::vector<glm::vec3> velocity;

        // Flags to indicate if the object is static or dynamic
        bool isStatic = true; 
        std::vector<bool> isVertStatic;

        float mass = 1.0f; // Default mass for the object
        float springRestitution = 5.0f; // Default restitution for the object
        float friction = 0.5f; // Default friction for the object
        float linearDamping = 0.005f; // Default linear damping for the object
        float angularDamping = 0.1f; // Default angular damping for the object
        float dragCoefficient = 0.47f; // Default drag coefficient for the object
        float restitution = 0.5f; // Default restitution for the object

        // Constructor to initialize the object with a mesh from a file
        object(const std::string& filename) {
            m = new mesh(filename);
            for (long long i = 0; i < m->mod->num_vertices; i++) {
                position.push_back(glm::vec3(m->mod->vertices[i].x, 
                                             m->mod->vertices[i].y, 
                                             m->mod->vertices[i].z));
            }
            velocity = std::vector<glm::vec3>(position.size(), glm::vec3(0.0f));
            isVertStatic = std::vector<bool>(position.size(), true); // Initialize all vertices as static
            
        }

        // Default destructor
        ~object() {
            if (m) {
                delete m; // Clean up
            }
        }

        // Function to set the positions of the vertices based on the mesh data
        void set_positions(){
            for (long long i = 0; i < m->mod->num_vertices; i++) {
                position[i] = glm::vec3(m->mod->vertices[i].x, 
                                        m->mod->vertices[i].y, 
                                        m->mod->vertices[i].z);
            }
        }

        // Function to update the velocities of the vertices based on the physics simulation
        void calculate_forces() {
            if (isStatic) return; // Skip force calculations for static objects
            // using every edge as a spring
            for (size_t i = 0; i < m->mod->num_vertices; i++) {
                long long v1_index = i;
                for (size_t j = i+1; j < m->mod->num_vertices; j++) {
                    long long v2_index = j;

                    glm::vec3 v1 = position[v1_index];
                    glm::vec3 v2 = position[v2_index];

                    glm::vec3 initial_v1 = glm::vec3(m->mod->vertices[v1_index].x, 
                                                      m->mod->vertices[v1_index].y, 
                                                      m->mod->vertices[v1_index].z);
                    glm::vec3 initial_v2 = glm::vec3(m->mod->vertices[v2_index].x, 
                                                      m->mod->vertices[v2_index].y, 
                                                      m->mod->vertices[v2_index].z);
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

                // Cheat to prevent objects from falling through the ground
                // This is a very simple collision detection and response
                if (position[v1_index].y < -0.9f) { // Simple ground collision
                    position[v1_index].y = -0.9f; // Reset position to ground level
                    velocity[v1_index].y = -velocity[v1_index].y * restitution; // Reverse velocity with restitution
                }
            }
        }

        // Function to update the positions of the vertices based on the velocities
        // This function also skips updates for static objects and static vertices
        void update_positions() {
            if (isStatic) return; // Skip position updates for static objects
            for (size_t i = 0; i < position.size(); i++) {
                if (isVertStatic[i]) continue; // Skip static vertices
                position[i] += velocity[i] * deltaTime;
            }
        }

        // Function to update the mesh data and draw the object
        void draw(Shader shader) {
            m->update_data(position);
            m->draw_mesh(shader);
        }
};

// Class that represents an Axis-Aligned Bounding Box (AABB) for collision detection
class AABB {
    public:

        // Minimum and maximum coordinates of the AABB
        glm::vec3 min;
        glm::vec3 max;

        // Number of faces contained within this AABB and a vector to store them
        long long num_faces = 0;
        std::vector<face*> faces;

        // Pointers to the object this AABB represents and its children
        object* obj = nullptr;
        AABB* children[2] = {nullptr};

        // Flag to indicate if this AABB is a leaf node in the octree
        bool isLeaf = true;

        // Default constructor initializes an empty AABB
        AABB() : obj(nullptr), min(glm::vec3(0.0f)), max(glm::vec3(0.0f)) {
            
            for (int i = 0; i < 2; ++i) {
                children[i] = nullptr;
            }
            num_faces = 0;
            isLeaf = true;
            faces =  std::vector<face*>();
        }

        // Constructor that initializes the AABB with an object from a file
        AABB(const std::string& filename) {
            obj = new object(filename);
            isLeaf = true; // Initially set to leaf
            set_faces();
            calculate_bounds();
            generate_children();
        }

        // Destructor to clean up the AABB and its children
        ~AABB() {
            for (int i = 0; i < 2; ++i) {
                if (children[i]) {
                    delete children[i]; // Clean up child AABBs
                }
            }
        }

        // Function to set the faces of the AABB based on the object's mesh
        void set_positions() {
            if (!obj) return;
            obj->set_positions();
        }

        // Function that detects collisions between this AABB and another AABB
        void detect_collision(AABB* other, std::vector<face*>& collision_faces) {
            if (!other || !obj || !other->obj) return;


            //std::cout << min.x << ' ' << min.y << ' ' << min.z << ' ' << max.x << ' ' << max.y << ' ' << max.z << std::endl;
            //std::cout << other->min.x << ' ' << other->min.y << ' ' << other->min.z << ' ' << other->max.x << ' ' << other->max.y << ' ' << other->max.z << std::endl;
            
            // Check for overlap in AABBs
            if (min.x > other->max.x || max.x < other->min.x ||
                min.y > other->max.y || max.y < other->min.y ||
                min.z > other->max.z || max.z < other->min.z) {
                return; // No collision
            }

            // If this AABB is a leaf node, check its faces against the other's faces
            //std::cout << "AABB 1 is leaf? " << (num_faces < 20 ? "Yes" : "No") << std::endl;
            //std::cout << "AABB 1 faces count: " << num_faces << std::endl;
            //std::cout << "AABB 2 is leaf? " << (other->num_faces < 20 ? "Yes" : "No") << std::endl;
            //std::cout << "AABB 2 faces count: " << other->num_faces << std::endl;
            
            // If both AABBs are leaf nodes, check their faces for collision
            if (num_faces <= 20 and other->num_faces <= 20) {
                //std::cout << "Both AABBs are leaf nodes, checking faces for collision." << std::endl;
                for (const auto& f : faces) {
                    for (const auto& of : other->faces) {
                        glm::vec3 tr1[3] = { obj->position[f->vertices[0]-1],
                                                obj->position[f->vertices[1]-1],
                                                obj->position[f->vertices[2]-1] };
                        glm::vec3 tr2[3] = {other->obj->position[of->vertices[0]-1],
                                                other->obj->position[of->vertices[1]-1],
                                                other->obj->position[of->vertices[2]-1] };
                        if(triangleTriangleIntersectionComplete(tr1, tr2)){ 
                            collision_faces.push_back(f);
                            collision_faces.push_back(of);
                        }
                    }
                }
            } else {
                // If not leaf, check children recursively
                
                //Generate children now to spare useless recursion
                if (num_faces > 20 && num_faces >= other->num_faces){
                    generate_children();
                    other->detect_collision(children[0], collision_faces);
                    other->detect_collision(children[1], collision_faces);
                }else if (other->num_faces > 20){
                    other->generate_children();
                    detect_collision(other->children[0], collision_faces);
                    detect_collision(other->children[1], collision_faces);
                }
            }
        }

        // Function to return the number of faces in the leaf nodes of the AABB
        int calc_faces(){
            if (isLeaf){
                return num_faces;
            }

            return children[0]->calc_faces() + children[1]->calc_faces();
        }

        // Function to identify if two triangles intersect
        bool triangleTriangleIntersectionComplete(const glm::vec3 tr1[3], const glm::vec3 tr2[3]) {

            glm::vec3 normal1 = glm::cross(tr1[1] - tr1[0], tr1[2] - tr1[0]);

            float d1 = -glm::dot(normal1, tr1[0]);
            //std::cout << "D1: " << d1 << std::endl;

            // Distances of triangle 2 vertices to plane 1
            float dist2[3] =  {0.0f, 0.0f, 0.0f};
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

                    for (int i = 0; i < 3; i++) {
                        float proj = glm::dot(axis, tr1[i]);
                        minProj1 = std::min(minProj1, proj);
                        maxProj1 = std::max(maxProj1, proj);
                    }

                    for (int i = 0; i < 3; i++) {
                        float proj = glm::dot(axis, tr2[i]);
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

        // Function to calculate the bounds of the AABB based on the faces it contains
        void calculate_bounds() {
            if (!obj->m) return;

            min = glm::vec3(INFINITY);
            max = glm::vec3(-INFINITY);

            for (const auto& face : faces) {
                for (int i  = 0; i < face->num_vertex; i++) {
                    const glm::vec3& vert = obj->position[face->vertices[i] - 1];
                    
                    // Update min and max bounds based on vertex positions
                    min = glm::min(min, vert);
                    max = glm::max(max,  vert);
                }
            }
        }

    private:

        // Function to set the faces of the AABB based on the object's mesh
        void set_faces(){
            for (int i = 0; i < obj->m->mod->num_faces; i++) {
                face* f = &obj->m->mod->faces[i];
                faces.push_back(f);
                num_faces++;
            }
        }

        // Function to check if a point is contained within the AABB
        bool contains(const glm::vec3& point) const {
            return (point.x >= min.x && point.x <= max.x) &&
                (point.y >= min.y && point.y <= max.y) &&
                (point.z >= min.z && point.z <= max.z);
        }

        // Funtion to find the faces that are inherited from the parent AABB
        void find_faces(std::vector<face*>& faces_pai, long long& num_faces_pai) {
            if (!obj->m) return;
            
            for (int i = 0; i < num_faces_pai; i++) {
                face* f = faces_pai[i];
                bool faceContained = false;
                for (int j = 0; j < f->num_vertex; j++) {
                    glm::vec3 vertexPos = glm::vec3(obj->m->mod->vertices[f->vertices[j] - 1].x,
                                                    obj->m->mod->vertices[f->vertices[j] - 1].y,
                                                    obj->m->mod->vertices[f->vertices[j] - 1].z);
                    //std::cout << "Checking vertex: " << min.x << " <= " << vertexPos.x << max.x << " ";
                    //std::cout << "Checking vertex: " << min.y << " <= " << vertexPos.y << max.y << " ";
                    //std::cout << "Checking vertex: " << min.z << " <= " << vertexPos.z << max.z << std::endl;
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

        // Function to generate child AABBs by splitting the current AABB
        void generate_children() {
            if (!obj->m) return;

            if (num_faces <= 20) {
                isLeaf = true;
                return;
            }
            calculate_bounds(); // Ensure bounds are calculated before splitting

            // Calculate split plane
            glm::vec3 center = (min + max) * 0.5f;

            //std::cout << center.x << ' ' << center.y << ' ' << center.z << std::endl;

            glm::vec3 size = max - min;
            
            //std::cout << size.x << ' ' << size.y << ' ' << size.z << std::endl;
            
            // Determine dominant axis for splitting
            int dominantAxis = 0;
            if (size.y >= size.x && size.y >= size.z) dominantAxis = 1;
            else if (size.z >= size.x && size.z >= size.y) dominantAxis = 2;
            
            // Create child bounds
            glm::vec3 childMin1 = min;
            glm::vec3 childMax1 = max;
            glm::vec3 childMin2 = min;
            glm::vec3 childMax2 = max;
            //std::cout << "Child 1 min: (" << childMin1.x << ", " << childMin1.y << ", " << childMin1.z << ")" << " ";
            //std::cout << "Child 1 max: (" << childMax1.x << ", " << childMax1.y << ", " << childMax1.z << ")" << std::endl;
            //std::cout << "Child 2 min: (" << childMin2.x << ", " << childMin2.y << ", " << childMin2.z << ")" << " ";
            //std::cout << "Child 2 max: (" << childMax2.x << ", " << childMax2.y << ", " << childMax2.z << ")" << std::endl;

            // Split along dominant axis
            switch (dominantAxis) {
                case 0: // X-axis
                    childMax1.x = center.x;
                    childMin2.x = center.x;
                    break;
                case 1: // Y-axis
                    childMax1.y = center.y;
                    childMin2.y = center.y;
                    break;
                case 2: // Z-axis
                    childMax1.z = center.z;
                    childMin2.z = center.z;
                    break;
            }

            children[0] = new AABB();
            children[0]->min = childMin1;
            children[0]->max = childMax1;
            children[0]->obj = obj; // Assign the same mesh to the first child
            children[0]->find_faces(faces, num_faces);



            children[1] = new AABB();
            children[1]->min = childMin2;
            children[1]->max = childMax2;
            children[1]->obj = obj; // Assign the same mesh to the second child
            children[1]->find_faces(faces, num_faces);            



            //faces.clear(); // Clear the faces in the parent AABB after distributing them to children
            //num_faces = 0; // Reset the face count in the parent AABB
            isLeaf = false; // This AABB is no longer a leaf node
            //std::cout << "AABB split into two children along axis " << dominantAxis << std::endl;
            //std::cout << "Child 1 faces: " << std::to_string(children[0]->num_faces) << std::endl;
            //std::cout << "Child 2 faces: " << std::to_string(children[1]->num_faces) << std::endl;
            //children[0]->generate_children(); // Recursively generate children for the first child
            //children[1]->generate_children(); // Recursively generate children for the second child
        }
};

// Class that represents the world containing multiple AABBs
class world {
    public:
        AABB** trees;
        int num_tr = 0;

        // Default constructor initializes the world with no trees
        world(){
            trees = nullptr;
        }

        // Destructor to clean up the world and its trees
        ~world() {
            for (int i = 0; i < num_tr; i++) {
                AABB* tr = trees[i];
                if (tr) {
                    delete tr->obj; // Clean up the object associated with the AABB
                    delete tr; // Clean up each object
                }
            }
            delete[] trees; // Clean up the array of objects
        }

        // Function to add a new tree (AABB) to the world
        void add_tree(const std::string& filename) {
            AABB** new_tr = new AABB*[num_tr + 1];
            for (int i = 0; i < num_tr; i++) {
                new_tr[i] = trees[i];
            }
            new_tr[num_tr] = new AABB(filename);
            delete[] trees; // Free the old array
            trees = new_tr;
            num_tr++;
        }

        // Function to initialize the world with a given number of trees
        void start_world() {
            for (int i = 0; i < num_tr; i++) {
                AABB* tr = trees[i];
                tr->obj->m->setup_mesh(i);
            }
        }

        // Function to draw the world using a given shader
        void draw_world(Shader shader) {
            update_world(); // Update positions and forces before drawing

            for (int i = 0; i < num_tr; i++) {
                // We pass the shader to the object's draw function so that we can pass uniform variables
                trees[i]->obj->draw(shader);
            }

        }

        // Function to update the world by calculating forces and updating positions of all trees
        // It also detects collisions after updating positions
        void update_world() {
            for (int i = 0; i < num_tr; i++) {
                trees[i]->obj->calculate_forces();
                trees[i]->obj->update_positions();
                trees[i]->calculate_bounds();
            }
            detect_collisions(); // Detect collisions after updating positions
        }

        // Function to detect collisions between all trees in the world
        // It checks for collisions only between non-static objects
        // If a collision is detected, both objects involved are marked as static
        void detect_collisions() {
            for (int i = 0; i < num_tr; i++) {
                AABB* tr1 = trees[i];
                for (int j = i + 1; j < num_tr; j++) {
                    if (tr1->obj->isStatic && trees[j]->obj->isStatic) {
                        continue; // Skip collision detection if both objects are static
                    }
                    AABB* tr2 = trees[j];
                    std::vector<face*> collision_faces = std::vector<face*>();
                    tr1->detect_collision(tr2, collision_faces);
                    if (!collision_faces.empty()) {
                        tr1->obj->isStatic = true; 
                        tr2->obj->isStatic = true;
                    }
                }
            }
        }
};


#endif // PHYSICS_H