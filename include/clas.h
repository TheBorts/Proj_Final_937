#ifndef CLASSES_H
#define CLASSES_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <fstream>

#include "shaders.h"

class material {
    public:
        std::string name;
        float ambient[3];
        float diffuse[3];
        float specular[3];
        float shininess;
        
        material() : name("Unnamed"),ambient{0.2f, 0.2f, 0.2f}, diffuse{0.8f, 0.8f, 0.8f}, specular{1.0f, 1.0f, 1.0f}, shininess(32.0f) {}
        
        material(std::string nam, float ar, float ag, float ab, float dr, float dg, float db, float sr, float sg, float sb, float sh)
            : name(nam),ambient{ar, ag, ab}, diffuse{dr, dg, db}, specular{sr, sg, sb}, shininess(sh) {}
};

class vertice{
    public:
        float x, y, z;
        vertice() : x(0.0f), y(0.0f), z(0.0f) {}
        vertice(float x, float y, float z) : x(x), y(y), z(z) {}
};

class normal{
    public:
        float x, y, z;
        normal() : x(0.0f), y(0.0f), z(0.0f) {}
        normal(float x, float y, float z) : x(x), y(y), z(z) {}
    };

class texture{
    public:
        float u, v, w;
        texture() : u(0.0f), v(0.0f), w(0.0f) {}
        texture(float u) : u(u), v(0.0f), w(0.0f) {}
        texture(float u, float v) : u(u), v(v), w(0.0f) {}
        texture(float u, float v, float w) : u(u), v(v), w(w) {}
};

class face{
    public:

        GLuint* vertices; // Indices of vertices;
        GLuint* normals; // Indices of normals
        GLuint* textures; // Indices of textures

        long long num_vertex;
        long long max_num_vertex;

        face(){
            vertices = new GLuint[5];
            normals = new GLuint[5];
            textures = new GLuint[5];
            num_vertex = 0;
            max_num_vertex = 5;
        }

        void add_vertex(int v_index, int n_index, int t_index) {
            if (num_vertex >= max_num_vertex) {
                max_num_vertex *= 2;
                GLuint* new_vertices = new GLuint[max_num_vertex];
                GLuint* new_normals = new GLuint[max_num_vertex];
                GLuint* new_textures = new GLuint[max_num_vertex];
                for (long long i = 0; i < num_vertex; i++) {
                    new_vertices[i] = vertices[i];
                    new_normals[i] = normals[i];
                    new_textures[i] = textures[i];
                }
                delete[] vertices;
                delete[] normals;
                delete[] textures;
                vertices = new_vertices;
                normals = new_normals;
                textures = new_textures;
            }
            vertices[num_vertex] = v_index;
            normals[num_vertex] = n_index;
            textures[num_vertex] = t_index;
            num_vertex++;
        }
};

class vertex{
    public:
        glm::vec3 position;
        glm::vec3 normals;
        vertex() : position(glm::vec3(0.0f)), normals(glm::vec3(0.0f)) {}
        vertex(vertice v, normal n) 
            : position(glm::vec3(v.x, v.y, v.z)), normals(glm::vec3(n.x, n.y, n.z)) {}
        void print() const {
            std::cout << "Vertex Position: (" << position.x << ", " << position.y << ", " << position.z << ") "
                      << "Normal: (" << normals.x << ", " << normals.y << ", " << normals.z << ")" << std::endl;
        }
};

class object3D {
    public:
        vertice* vertices;
        normal* normals;
        texture* textures;
        face* faces;
        face* init_faces;
        GLuint* indices;
        long long num_vertices;
        long long num_textures;
        long long num_normals;
        long long num_faces;
        long long init_num_faces;
        long long max_num_faces;
        vertice center = vertice();
        std::string obj_name = "";
        std::string path_material = "";
        material mat = material();

        std::ifstream myFile;
        object3D() : vertices(nullptr), normals(nullptr), faces(nullptr), num_vertices(0), num_normals(0), num_faces(0) {}
        object3D(const std::string& filename) : vertices(nullptr), textures(nullptr), normals(nullptr), faces(nullptr), indices(nullptr), num_vertices(0), num_normals(0), num_faces(0), max_num_faces(0) {
            myFile.open(filename);
            if (!myFile.is_open()) {
                std::cerr << "Error opening file: " << filename << std::endl;
                return;
            }
            std::string line;
            while (getline(myFile, line)) {
                if (line.substr(0, 2) == "v ") {
                    num_vertices++;
                } else if (line.substr(0, 3) == "vn ") {
                    num_normals++;
                } else if (line.substr(0, 3) == "vt ") {
                    num_textures++;
                } else if (line.substr(0, 2) == "f ") {
                    num_faces++;
                }
            }
            
            init_num_faces = num_faces;

            vertices = new vertice[num_vertices];
            normals = new normal[num_normals];
            textures = new texture[num_textures];
            faces = new face[num_faces];
            init_faces = new face[num_faces];
            indices = new GLuint[3 * num_faces];
            
            max_num_faces = num_faces;

            myFile.clear();
            myFile.seekg(0, std::ios::beg);
            
            long long v_index = 0, n_index = 0, f_index = 0, t_index = 0;

            while (getline(myFile, line)) {
                if (line.substr(0, 2) == "v ") {
                    get_vertices(line, vertices, v_index++);
                } else if (line.substr(0, 3) == "vn ") {
                    get_normals(line, normals, n_index++);
                } else if (line.substr(0, 3) == "vt ") {
                    get_textures(line, textures, t_index++);
                } else if (line.substr(0, 2) == "f ") {
                    get_faces(line, faces, f_index++);
                } else if (line.substr(0, 2) == "o ") {
                    get_name(line);
                } else if (line.substr(0, 7) == "mtllib ") {
                    get_material_path(line);
                } else if (line.substr(0, 7) == "usemtl ") {
                    get_material(line);
                }
            }

            myFile.close();

            for(long long i = 0; i < num_vertices; i++) {
                center.x += vertices[i].x;
                center.y += vertices[i].y;
                center.z += vertices[i].z;
            }
            center.x /= num_vertices;
            center.y /= num_vertices;
            center.z /= num_vertices;

            for (long long i = 0; i < num_faces; i++) {
                init_faces[i] = faces[i]; // Initialize the initial faces
            }

            // Triangulate faces if necessary
            triangulate_faces();

            for (long long i = 0; i < num_faces; i++) {
                indices[3 * i] = faces[i].vertices[0] - 1;
                indices[3 * i + 1] = faces[i].vertices[1] - 1;
                indices[3 * i + 2] = faces[i].vertices[2] - 1;
            }
        }

        void translate(float dx, float dy, float dz) {
            for (long long i = 0; i < num_vertices; i++) {
                vertices[i].x += dx;
                vertices[i].y += dy;
                vertices[i].z += dz;
            }
            center.x += dx;
            center.y += dy;
            center.z += dz;
        }

        void rotate(float theta, float phi, float psi) {
            glm::mat3 R = glm::mat3(glm::rotate(glm::mat4(1.0f), glm::radians(theta), glm::vec3(1.0f, 0.0f, 0.0f)) *
                                    glm::rotate(glm::mat4(1.0f), glm::radians(phi), glm::vec3(0.0f, 1.0f, 0.0f)) *
                                    glm::rotate(glm::mat4(1.0f), glm::radians(psi), glm::vec3(0.0f, 0.0f, 1.0f)));
            for (long long i = 0; i < num_vertices; i++) {
                glm::vec3 v = glm::vec3(vertices[i].x, vertices[i].y, vertices[i].z);
                v = R * v;
                vertices[i].x = v.x;
                vertices[i].y = v.y;
                vertices[i].z = v.z;
            }
        }      

    private:

        void get_name(std::string line) {
            line = line.substr(2); // Remove "o "
            obj_name = line;
        }

        void get_material_path(std::string line) {
            line = line.substr(7); // Remove "mtllib "
            path_material = line; // Store the material file path
        }

        void get_material(std::string line){
            line = line.substr(7); // Remove "usemtl "
            mat.name = line; // Store the material name or path
            std::ifstream matFile(path_material);
            if (!matFile.is_open()) {
                std::cerr << "Error opening material file: " << path_material << std::endl;
                return;
            }
            std::string matLine;
            while (getline(matFile, matLine)) {
                if (matLine.substr(0, 2) == "Ka") { // Ambient color
                    sscanf(matLine.c_str(), "Ka %f %f %f", &mat.ambient[0], &mat.ambient[1], &mat.ambient[2]);
                } else if (matLine.substr(0, 2) == "Kd") { // Diffuse color
                    sscanf(matLine.c_str(), "Kd %f %f %f", &mat.diffuse[0], &mat.diffuse[1], &mat.diffuse[2]);
                } else if (matLine.substr(0, 2) == "Ks") { // Specular color
                    sscanf(matLine.c_str(), "Ks %f %f %f", &mat.specular[0], &mat.specular[1], &mat.specular[2]);
                } else if (matLine.substr(0, 2) == "Ns") { // Shininess
                    sscanf(matLine.c_str(), "Ns %f", &mat.shininess);
                }
            }
            matFile.close();
        }

        void get_vertices(std::string line, vertice* vertices, int index_vertices){
            std::string mynum;

            line = line.substr(2);

            for(int i = 0; i < 3; i++){
                int j = 0;
                while(line[j] != ' ' && j < line.size()){
                    mynum += line[j];
                    j++;
                }
                *((float*)(vertices + index_vertices) + i) = stof(mynum);
                mynum = "";
                if (j < line.size()) line = line.substr(j+1);
            }
        }

        void get_textures(std::string line, texture* textures, int index_textures){
            std::string mynum;

            line = line.substr(3);

            for(int i = 0; i < 3; i++){
                int j = 0;
                while(line[j] != ' ' && j < line.size()){
                    mynum += line[j];
                    j++;
                }
                *((float*)(textures + index_textures) + i) = stof(mynum);
                mynum = "";
                if (j < line.size()) line = line.substr(j+1);
                if (line == "") break; // If there are no more texture coordinates, break
            }
        }
        
        void get_normals(std::string line, normal* normals, int index_normals){
            std::string mynum;

            line = line.substr(3);

            for(int i = 0; i < 3; i++){
                int j = 0;
                while(line[j] != ' ' && j < line.size()){
                    mynum += line[j];
                    j++;
                }
                *((float*)(normals + index_normals) + i) = stof(mynum);
                mynum = "";
                if (j < line.size()) line = line.substr(j+1);
            }
        }
        
        void get_faces(std::string line, face* faces, int index_faces){
            std::string mynum= "0";

            line = line.substr(2);
            line += " "; // Add a space at the end to ensure the last vertex is processed
            for(int i = 0; i < 15; i++){
                int j = 0;
                int index[3] = {0, 0, 0};
                for (int k = 0; k < 3; k++) {
                    while(line[j] != ' ' && j < line.size() && line[j] != '/'){
                        mynum += line[j];
                        j++;
                    }
                    index[k] = stoi(mynum);
                    mynum = "0";
                    line = line.substr(j+1);
                    j = 0; // Reset j for the next index
                }
                faces[index_faces].add_vertex(index[0], index[1], index[2]);
                if (line.empty()) return; // If there are no more vertices in the face, break)
            }
        }

        void add_face(face f){
            if (num_faces >= max_num_faces) {
                max_num_faces *= 2;
                face* new_faces = new face[max_num_faces];
                for (long long i = 0; i < num_faces; i++) {
                    new_faces[i] = faces[i];
                }
                delete[] faces;
                faces = new_faces;
            }
            faces[num_faces++] = f;
        }

        void triangulate_faces() {
            int num_faces_initial = num_faces;
            for (long long i = 0; i < num_faces_initial; i++) {
                if (faces[i].num_vertex > 3) {
                    // Triangulate the face
                    GLuint first_vertex = faces[i].vertices[0];
                    face last_face = face();
                    
                    last_face.add_vertex(first_vertex, faces[i].normals[0], faces[i].textures[0]);
                    last_face.add_vertex(faces[i].vertices[1], faces[i].normals[1], faces[i].textures[1]);
                    last_face.add_vertex(faces[i].vertices[2], faces[i].normals[2], faces[i].textures[2]);

                    for (long long j = 2; j < faces[i].num_vertex - 1; j++) {
                        face n_face = face();
                        n_face.add_vertex(first_vertex, faces[i].normals[0], faces[i].textures[0]);
                        n_face.add_vertex(faces[i].vertices[j], faces[i].normals[j], faces[i].textures[j]);
                        n_face.add_vertex(faces[i].vertices[j + 1], faces[i].normals[j + 1], faces[i].textures[j + 1]);
                        add_face(n_face);
                    }

                    faces[i] = last_face; // Replace the original face with the first triangle
                }
            }
        }
};

class mesh{
    public:
        object3D* obj;
        normal* calculated_normals;
        vertex* data;
        GLuint VAO, VBO, EBO;
        
        mesh() : obj(nullptr) {}
        mesh(const std::string& filename) {
            obj = new object3D(filename);
        }

        ~mesh() {
            delete obj;
        }

        void set_data() {
            data = new vertex[obj->num_vertices];
            for (long long i = 0; i < obj->num_vertices; i++) {
                data[i] = vertex(obj->vertices[i], calculated_normals[i]);
                data[i].print(); // Print vertex data for debugging
            }
        }

        void calculate_normals() {
            calculated_normals = new normal[obj->num_vertices];
            // Initialize normals to zero
            for (long long i = 0; i < obj->num_vertices; i++) {
                calculated_normals[i] = normal(0.0f, 0.0f, 0.0f);
            }
            // Calculate normals for each face
            for (long long i = 0; i < obj->init_num_faces; i++) {
                face& f = obj->init_faces[i];
                if (f.num_vertex < 3) continue; // Skip faces with less than 3 vertices
                glm::vec3 edge1 = glm::vec3(obj->vertices[f.vertices[1] - 1].x - obj->vertices[f.vertices[0] - 1].x,
                                                obj->vertices[f.vertices[1] - 1].y - obj->vertices[f.vertices[0] - 1].y,
                                                obj->vertices[f.vertices[1] - 1].z - obj->vertices[f.vertices[0] - 1].z);
                glm::vec3 edge2 = glm::vec3(obj->vertices[f.vertices[1] - 1].x - obj->vertices[f.vertices[2] - 1].x,
                                                obj->vertices[f.vertices[1] - 1].y - obj->vertices[f.vertices[2] - 1].y,
                                                obj->vertices[f.vertices[1] - 1].z - obj->vertices[f.vertices[2] - 1].z);
                
                glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));

                for (int j = 0; j < f.num_vertex; j++) {
                    long long vertex_index = f.vertices[j] - 1; // Convert to zero-based index
                    calculated_normals[vertex_index].x += normal.x;
                    calculated_normals[vertex_index].y += normal.y;
                    calculated_normals[vertex_index].z += normal.z;
                }
                
            }
            // Normalize the normals
            for (long long i = 0; i < obj->num_vertices; i++) {
                float length = glm::length(glm::vec3(calculated_normals[i].x, calculated_normals[i].y, calculated_normals[i].z));
                if (length > 0.0f) {                
                    calculated_normals[i].x /= length;
                    calculated_normals[i].y /= length;
                    calculated_normals[i].z /= length;
                } else {
                    calculated_normals[i] = normal(0.0f, 0.0f, 0.0f); // Set to zero if length is zero
                }
            }
        }

        void setup_mesh(){
            calculate_normals(); // Calculate normals before setting up the mesh
            set_data(); // Set the vertex data

            glGenVertexArrays(1, &VAO);
            glGenBuffers(1, &VBO);
            glGenBuffers(1, &EBO);

            glBindVertexArray(VAO); // Make it active
            glBindBuffer(GL_ARRAY_BUFFER, VBO); // Make it active
            
            glBufferData(GL_ARRAY_BUFFER, 6*obj->num_vertices * sizeof(float), data, GL_STATIC_DRAW);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * obj->num_faces * sizeof(unsigned int), obj->indices, GL_STATIC_DRAW);

            //Defines how vertex data is structured
            glVertexAttribPointer(0,// Attribute location in the shader
                          3,//Number of components per vertex
                          GL_FLOAT,//Data type of each component
                          GL_FALSE,//Normalize integer values to [0,1]
                          6*sizeof(float),//Byte offset between consecutive attributes
                          nullptr);//Offset from the start of VBO where data begins;
            glEnableVertexAttribArray(0); //Enables the attribute (position)

            //Defines how normal data is structured
            glVertexAttribPointer(1, // Attribute location in the shader
                          3, // Number of components per normal
                          GL_FLOAT, // Data type of each component
                          GL_FALSE, // Normalize integer values to [0,1]
                          6*sizeof(float), // Byte offset between consecutive attributes
                          (void*)(3*sizeof(float))); // Offset from the start of VBO where normals begin
            glEnableVertexAttribArray(1); // Enables the attribute (normal)
        }

        void draw_mesh(Shader shader) {

            shader.setVec3("material.ambient", obj->mat.ambient[0], obj->mat.ambient[1], obj->mat.ambient[2]);
            shader.setVec3("material.diffuse", obj->mat.diffuse[0], obj->mat.diffuse[1], obj->mat.diffuse[2]);
            shader.setVec3("material.specular", obj->mat.specular[0], obj->mat.specular[1], obj->mat.specular[2]);
            shader.setFloat("material.shininess", obj->mat.shininess);

            glBindVertexArray(VAO); // Bind the VAO
            //glDrawArrays(GL_TRIANGLES, 0, 6*obj->num_vertices); // Draw the mesh using vertex array
            glDrawElements(GL_TRIANGLES, 3 * obj->num_faces, GL_UNSIGNED_INT, 0); // Draw the mesh
            glBindVertexArray(0); // Unbind the VAO
        }

        void translate(float dx, float dy, float dz) {
            if (obj) {
                obj->translate(dx, dy, dz);
            }
        }

        void rotate(float theta, float phi, float psi) {
            if (obj) {
                obj->rotate(theta, phi, psi);
            }
        }

        object3D* get_object() const {
            return obj;
        }
};

class lightSource {
    public:
        float position[3];
        float ambient[3];
        float diffuse[3];
        float specular[3];
        lightSource() : position{0.0f, 0.0f, 0.0f}, ambient{0.2f, 0.2f, 0.2f}, diffuse{0.8f, 0.8f, 0.8f}, specular{1.0f, 1.0f, 1.0f} {}
        lightSource(float x, float y, float z, float ar, float ag, float ab, float dr, float dg, float db, float sr, float sg, float sb) 
            : position{x, y, z}, ambient{ar, ag, ab}, diffuse{dr, dg, db}, specular{sr, sg, sb} {}
};

class scene {
    public:
        mesh** meshes;
        int num_meshes;
        scene() : meshes(nullptr), num_meshes(0) {}
        
        void add_mesh(const std::string& filename) {
            mesh** new_meshes = new mesh*[num_meshes + 1];
            for (int i = 0; i < num_meshes; i++) {
                new_meshes[i] = meshes[i];
            }
            new_meshes[num_meshes] = new mesh(filename);
            delete[] meshes;
            meshes = new_meshes;
            num_meshes++;
        }

        void start_scene() {
            for (int i = 0; i < num_meshes; i++) {
                meshes[i]->setup_mesh();    
            }
        }

        void draw_scene(Shader shader) {
            for (int i = 0; i < num_meshes; i++) {
                meshes[i]->draw_mesh(shader);
            }
        }
};

#endif