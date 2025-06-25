#ifndef CLASSES_H
#define CLASSES_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <fstream>
#include <filesystem>
#include <vector>
#include <unordered_map>

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
        float x = 0, y = 0, z = 0;
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

        face& operator=(const face& other) {
            num_vertex = other.num_vertex;
            max_num_vertex = other.max_num_vertex;
            vertices = new GLuint[max_num_vertex];
            normals = new GLuint[max_num_vertex];
            textures = new GLuint[max_num_vertex];
            for (long long i = 0; i < num_vertex; i++) {
                vertices[i] = other.vertices[i];
                normals[i] = other.normals[i];
                textures[i] = other.textures[i];
            }
            return *this;
        }

        ~face() {
            delete[] vertices;
            delete[] normals;
            delete[] textures;
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
        glm::vec2 texture_coords;
        vertex() : position(glm::vec3(0.0f)), normals(glm::vec3(0.0f)), texture_coords(glm::vec2(0.0f)) {}
        vertex(vertice v, normal n)
            : position(glm::vec3(v.x, v.y, v.z)), normals(glm::vec3(n.x, n.y, n.z)), texture_coords(glm::vec2(0.0f)) {}
        vertex(vertice v, normal n, texture t)
            : position(glm::vec3(v.x, v.y, v.z)), normals(glm::vec3(n.x, n.y, n.z)), texture_coords(glm::vec2(t.u, t.v)) {}
        vertex(glm::vec3 pos, glm::vec3 norm, glm::vec2 tex)
            : position(pos), normals(norm), texture_coords(tex) {}
        
        void print() const {
            std::cout << "Vertex Position: (" << position.x << ", " << position.y << ", " << position.z << ") "
                      << "Normal: (" << normals.x << ", " << normals.y << ", " << normals.z << ") "
                      << "Texture: (" << texture_coords.x << ", " << texture_coords.y << ")" << std::endl;
        }
        void print_textures() const {
            std::cout << "Texture Coordinates: (" << texture_coords.x << ", " << texture_coords.y << ")" << std::endl;
        }
};

class model3D {
    public:
        vertice* vertices = nullptr;
        
        texture* textures = nullptr;
        normal* normals = nullptr;
        face* faces = nullptr;

        face* init_faces = nullptr;
        GLuint* indices = nullptr;
        
        long long num_vertices = 0;
        long long num_textures = 0;
        long long num_normals = 0;
        long long num_faces = 0;
        long long init_num_faces = 0;
        long long max_num_faces = 0;
        
        vertice center = vertice();
        
        std::string mod_name = "";

        std::string path_material = "";
        material mat = material();
        
        std::string path_texture = "";

        std::ifstream myFile = std::ifstream();

        model3D() {}
        model3D(const std::string& filename){
            
            myFile.open(filename);
            if (myFile.fail()) {
                std::cerr << "Error opening file: " << filename << std::endl;
                return;
            }

            std::string line = "";
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
            
            max_num_faces = num_faces;

            if (num_textures == 0) {
                delete[] textures; // If no textures are defined, delete the texture array
                textures = new texture[1] ; // Create a default texture
                textures[0] = texture(0.0f, 0.0f); // Default texture coordinates
            }

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

            indices = new GLuint[3 * num_faces];
            
            for (long long i = 0; i < num_faces; i++) {
                indices[3 * i] = faces[i].vertices[0] - 1;
                indices[3 * i + 1] = faces[i].vertices[1] - 1;
                indices[3 * i + 2] = faces[i].vertices[2] - 1;
            }
        }

        ~model3D() {
            delete[] vertices;
            delete[] normals;
            delete[] textures;
            delete[] indices;
            delete[] faces;
            delete[] init_faces;
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

        void trim(std::string& str) {
            // Remove leading and trailing whitespace
            str.erase(0, str.find_first_not_of(" \t\n\r\f\v"));
            str.erase(str.find_last_not_of(" \t\n\r\f\v") + 1);
        }

        void get_name(std::string line) {
            line = line.substr(2); // Remove "o "
            mod_name = line;
        }

        void get_material_path(std::string line) {
            line = line.substr(7); // Remove "mtllib "
            std::string dir = "../models/";
            path_material = dir + line; // Store the material file path
            trim(path_material); // Trim whitespace
            //path_material.pop_back(); // Remove the newline character at the end
            //std::cout << "Material file path: " << path_material << std::endl;
            if (!std::filesystem::exists(path_material)) {
                std::cout << path_material << std::endl;
                std::cerr << "Please check the path or the file name. You're currently on " << std::filesystem::current_path() << std::endl;
                return;
            }
        }

        void get_material(std::string line){
            line = line.substr(7); // Remove "usemtl "
            mat.name = line; // Store the material name or path
            std::ifstream matFile(path_material);
            if (matFile.fail()) {
                std::cerr << "Error opening material file: " << path_material << std::endl;
                return;
            }
            std::string matLine = "";
            while (getline(matFile, matLine)) {
                if (matLine.substr(0, 2) == "Ka") { // Ambient color
                    sscanf(matLine.c_str(), "Ka %f %f %f", &mat.ambient[0], &mat.ambient[1], &mat.ambient[2]);
                } else if (matLine.substr(0, 2) == "Kd") { // Diffuse color
                    sscanf(matLine.c_str(), "Kd %f %f %f", &mat.diffuse[0], &mat.diffuse[1], &mat.diffuse[2]);
                } else if (matLine.substr(0, 2) == "Ks") { // Specular color
                    sscanf(matLine.c_str(), "Ks %f %f %f", &mat.specular[0], &mat.specular[1], &mat.specular[2]);
                } else if (matLine.substr(0, 2) == "Ns") { // Shininess
                    sscanf(matLine.c_str(), "Ns %f", &mat.shininess);
                } else if (matLine.substr(0, 6) == "map_Kd") { // Texture map
                    path_texture = matLine.substr(7); // Remove "map_Kd "
                    path_texture = "../models/" + path_texture; // Prepend the directory
                    trim(path_texture); // Remove the newline character at the end
                }
            }
            matFile.close();
        }

        void get_vertices(std::string line, vertice* vertices, int index_vertices){
            std::string mynum = "";

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
            std::string mynum = "";

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
            std::string mynum = "";

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
            std::string mynum= "";

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

                    if (!mynum.empty()){
                        index[k] = stoi(mynum);
                    }else{
                        index[k] = 1; // Default to 1 if no index is found
                    }
                    mynum = "";
                    line = line.substr(j+1);
                    j = 0; // Reset j for the next index
                }
                faces[index_faces].add_vertex(index[0], index[2], index[1]);
                if (line.empty()) return; // If there are no more vertices in the face, break)
            }
        }

        void add_face(face& f){
            if (num_faces >= max_num_faces) {
                max_num_faces *= 2;
                face* new_faces = new face[max_num_faces];
                for (long long i = 0; i < num_faces; i++) {
                    new_faces[i] = faces[i];
                }
                delete[] faces; // Free the old array
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

                    delete [] faces[i].vertices;
                    delete [] faces[i].normals;
                    delete [] faces[i].textures;
                    faces[i] = last_face; // Replace the original face with the first triangle
                }
            }
        }
};

struct VertexKey {
    unsigned int v, vt, vn;
    bool operator==(const VertexKey& other) const {
        return v == other.v && vt == other.vt && vn == other.vn;
    }
    // Hash function for std::unordered_map
    struct Hash {
        size_t operator()(const VertexKey& k) const {
            return ((k.v << 16) ^ (k.vt << 8) ^ k.vn);
        }
    };
};

class mesh{
    public:
        model3D* mod = nullptr;
        vertice* proper_vertices = nullptr;
        normal* proper_normals = nullptr;
        texture* proper_textures = nullptr;
        std::unordered_map<VertexKey, GLuint, VertexKey::Hash> vertex_map = std::unordered_map<VertexKey, GLuint, VertexKey::Hash>();
        std::vector<GLuint> indices = std::vector<GLuint>();
        std::vector<vertex> data = std::vector<vertex>();
        GLuint VAO, VBO, EBO;
        unsigned int textureID;
        
        mesh() : mod(nullptr) {}
        mesh(const std::string& filename) {
            mod = new model3D(filename);
        }

        ~mesh() {
           
            if (mod) delete mod; // Clean up the model3D object
            delete[] proper_vertices; // Clean up vertices
            delete[] proper_normals; // Clean up normals
            delete[] proper_textures; // Clean up textures

            glDeleteVertexArrays(1, &VAO); // Delete the Vertex Array Object
            glDeleteBuffers(1, &VBO); // Delete the Vertex Buffer Object
            glDeleteBuffers(1, &EBO); // Delete the Element Buffer Object
            if (textureID) {
                glDeleteTextures(1, &textureID); // Delete the texture if it exists
            }
        }

        void set_data() {
            for (long long i = 0; i < mod->num_faces; i++) {
                for (long long j = 0; j < mod->faces[i].num_vertex; j++) {                    
                    VertexKey key{mod->faces[i].vertices[j], 
                        mod->faces[i].textures[j], 
                        mod->faces[i].normals[j]};
                    if (vertex_map.find(key) == vertex_map.end()) {
                        // If the vertex is not in the map, add it
                        vertex v(mod->vertices[mod->faces[i].vertices[j] - 1], 
                            proper_normals[mod->faces[i].vertices[j] - 1], 
                            mod->textures[mod->faces[i].textures[j] - 1]);
                            data.push_back(v);
                            vertex_map[key] = data.size() - 1; // Store the index of the new vertex
                    }
                    indices.push_back(vertex_map[key]);
                }
            }
        }

        void calculate_normals(){
            std::vector<glm::vec3> verts =  std::vector<glm::vec3>(mod->num_vertices);
            for (long long i = 0; i < mod->num_vertices; i++) {
                verts[i] = glm::vec3(mod->vertices[i].x, mod->vertices[i].y, mod->vertices[i].z);
            }
            calculate_normals(verts); // Calculate normals using the vertex positions
        }

        void calculate_normals(std::vector<glm::vec3> verts) {
            proper_normals = new normal[mod->num_vertices];
            // Initialize normals to zero
            for (long long i = 0; i < mod->num_vertices; i++) {
                proper_normals[i] = normal(0.0f, 0.0f, 0.0f);
            }
            // Calculate normals for each face
            for (long long i = 0; i < mod->init_num_faces; i++) {
                face& f = mod->init_faces[i];
                if (f.num_vertex < 3) continue; // Skip faces with less than 3 vertices
                //std::cout << "Calculating normals for face " << i << " with " << f.num_vertex << " vertices." << std::endl;
                glm::vec3 edge1 = glm::vec3(verts[f.vertices[1] - 1].x - verts[f.vertices[0] - 1].x,
                                                verts[f.vertices[1] - 1].y - verts[f.vertices[0] - 1].y,
                                                verts[f.vertices[1] - 1].z - verts[f.vertices[0] - 1].z);
                glm::vec3 edge2 = glm::vec3(verts[f.vertices[1] - 1].x - verts[f.vertices[2] - 1].x,
                                                verts[f.vertices[1] - 1].y - verts[f.vertices[2] - 1].y,
                                                verts[f.vertices[1] - 1].z - verts[f.vertices[2] - 1].z);
                
                glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));

                for (int j = 0; j < f.num_vertex; j++) {
                    long long vertex_index = f.vertices[j] - 1; // Convert to zero-based index
                    proper_normals[vertex_index].x += normal.x;
                    proper_normals[vertex_index].y += normal.y;
                    proper_normals[vertex_index].z += normal.z;
                }
            }

            // Normalize the normals
            for (long long i = 0; i < mod->num_vertices; i++) {
                float length = glm::length(glm::vec3(proper_normals[i].x, proper_normals[i].y, proper_normals[i].z));
                if (length > 0.0f) {                
                    proper_normals[i].x /= length;
                    proper_normals[i].y /= length;
                    proper_normals[i].z /= length;
                } else {
                    proper_normals[i] = normal(0.0f, 0.0f, 0.0f); // Set to zero if length is zero
                }
            }
        }

        void set_texture(unsigned int i){
            if(mod->path_texture != "") {

                //float borderColor[] = {1.0f, 0.0f, 1.0f, 1.0f}; // Yellow border color

                glGenTextures(1, &textureID);
                glActiveTexture(GL_TEXTURE0 + i); // Activate texture unit i
                glBindTexture(GL_TEXTURE_2D, textureID); // Bind the texture
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
                //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
                //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
                //glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
                glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                
                int width, height, nrChannels;
                std::string dir = "../models/";
                stbi_set_flip_vertically_on_load(true);
                unsigned char* textureData = stbi_load(dir.append(mod->path_texture).c_str(), &width, &height, &nrChannels, 0);
                if (textureData) {
                    GLenum format = (nrChannels == 4) ? GL_RGBA : GL_RGB;
                    glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, textureData);
                    std::cout << "Loaded texture: " << mod->path_texture << " with size: " << width << "x" << height << std::endl;
                    glGenerateMipmap(GL_TEXTURE_2D);
                } else {
                    std::cerr << "Failed to load texture: " << mod->path_texture << std::endl;
                }
                stbi_image_free(textureData); // Free the texture data
            } else {
                std::cout << "No texture provided for mesh: " << mod->path_texture << std::endl;
                textureID = 0; // Set textureID to 0 if no texture is provided
            }
        }

        void setup_mesh(unsigned int i){
            calculate_normals(); // Calculate normals before setting up the mesh
            set_data(); // Set the vertex data
            set_texture(i); // Set the texture if available
            glGenVertexArrays(1, &VAO);
            glGenBuffers(1, &VBO);
            glGenBuffers(1, &EBO);

            glBindVertexArray(VAO); // Make it active
            glBindBuffer(GL_ARRAY_BUFFER, VBO); // Make it active

            glBufferData(GL_ARRAY_BUFFER, sizeof(vertex)*data.size(), data.data(), GL_STATIC_DRAW);

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int)*indices.size(), indices.data(), GL_STATIC_DRAW);

            //Defines how vertex data is structured
            glVertexAttribPointer(0,// Attribute location in the shader
                          3,//Number of components per vertex
                          GL_FLOAT,//Data type of each component
                          GL_FALSE,//Normalize integer values to [0,1]
                          8*sizeof(float),//Byte offset between consecutive attributes
                          nullptr);//Offset from the start of VBO where data begins;
            glEnableVertexAttribArray(0); //Enables the attribute (position)

            //Defines how normal data is structured
            glVertexAttribPointer(1, // Attribute location in the shader
                          3, // Number of components per normal
                          GL_FLOAT, // Data type of each component
                          GL_FALSE, // Normalize integer values to [0,1]
                          8*sizeof(float), // Byte offset between consecutive attributes
                          (void*)(3*sizeof(float))); // Offset from the start of VBO where normals begin
            glEnableVertexAttribArray(1); // Enables the attribute (normal)

            // Defines how texture coordinate data is structured
            glVertexAttribPointer(2, // Attribute location in the shader
                            2, // Number of components per texture coordinate
                            GL_FLOAT, // Data type of each component
                            GL_FALSE, // Normalize integer values to [0,1]
                            8*sizeof(float), // Byte offset between consecutive attributes
                            (void*)(6*sizeof(float))); // Offset from the start of VBO where texture coordinates
            glEnableVertexAttribArray(2); // Enables the attribute (texture coordinates)

            glBindBuffer(GL_ARRAY_BUFFER, 0); // Unbind the VBO
            glBindVertexArray(0); // Unbind the VAO
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); // Unbind the EBO                        
                
        }

        void update_data(std::vector<glm::vec3> verts) {
            if (proper_normals) {
                delete[] proper_normals; // Free existing normals if they exist
                proper_normals = nullptr; // Set to nullptr to avoid dangling pointer
            }
            calculate_normals(verts); // Calculate normals based on the vertices
            data.clear(); // Clear existing data
            vertex_map.clear(); // Clear the vertex map
            data = std::vector<vertex>(); // Reinitialize the data vector
            vertex_map = std::unordered_map<VertexKey, GLuint, VertexKey::Hash>(); // Reinitialize the vertex map
            for (long long i = 0; i < mod->num_faces; i++) {
                for (long long j = 0; j < mod->faces[i].num_vertex; j++) {
                    VertexKey key{mod->faces[i].vertices[j], 
                        mod->faces[i].textures[j], 
                        mod->faces[i].normals[j]};
                    if (vertex_map.find(key) == vertex_map.end()) {
                        // If the vertex is not in the map, add it
                        vertex v(verts[mod->faces[i].vertices[j] - 1],
                            glm::vec3(proper_normals[mod->faces[i].vertices[j] - 1].x, 
                                      proper_normals[mod->faces[i].vertices[j] - 1].y, 
                                      proper_normals[mod->faces[i].vertices[j] - 1].z),
                            glm::vec2(mod->textures[mod->faces[i].textures[j] - 1].u, mod->textures[mod->faces[i].textures[j] - 1].v));
                        data.push_back(v);
                        vertex_map[key] = data.size() - 1; // Store the index of the new vertex
                    }
                }
            }
        }

        void draw_mesh(Shader shader) {

            shader.setVec3("material.ambient", mod->mat.ambient[0], mod->mat.ambient[1], mod->mat.ambient[2]);
            shader.setVec3("material.diffuse", mod->mat.diffuse[0], mod->mat.diffuse[1], mod->mat.diffuse[2]);
            shader.setVec3("material.specular", mod->mat.specular[0], mod->mat.specular[1], mod->mat.specular[2]);
            shader.setFloat("material.shininess", mod->mat.shininess);

            if (textureID != 0) {
                glActiveTexture(GL_TEXTURE0); // Activate texture unit 0
                glBindTexture(GL_TEXTURE_2D, textureID); // Bind the texture
                shader.setInt("material.hasTexture", 1); // Set the texture uniform in the shader
                shader.setInt("OurTexture", 0); // Set the texture unit to 0
            } else {
                shader.setInt("material.hasTexture", 0); // No texture
            }
            
            glBindVertexArray(VAO); // Bind the VAO
            if (proper_normals) {
                glBindBuffer(GL_ARRAY_BUFFER, VBO); // Bind the VBO
                glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * data.size(), data.data(), GL_STATIC_DRAW); // Update the VBO with new data
            }
            //glDrawArrays(GL_TRIANGLES, 0, 6*mod->num_vertices); // Draw the mesh using vertex array
            glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0); // Draw the mesh
            glBindVertexArray(0); // Unbind the VAO
        }

        void translate(float dx, float dy, float dz) {
            if (mod) {
                mod->translate(dx, dy, dz);
            }
        }

        void rotate(float theta, float phi, float psi) {
            if (mod) {
                mod->rotate(theta, phi, psi);
            }
        }

        model3D* get_model() const {
            return mod;
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

#endif