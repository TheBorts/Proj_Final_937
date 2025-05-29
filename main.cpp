#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <filesystem>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

class object3D {
    public:
        GLfloat* vertices;
        GLfloat* normals;
        GLuint* faces;
        GLuint* indices;
        float* center;
        long long num_vertices;
        long long num_normals;
        long long num_faces;

        ifstream myFile;
        object3D() : vertices(nullptr), normals(nullptr), faces(nullptr), num_vertices(0), num_normals(0), num_faces(0) {}
        object3D(const string& filename) : vertices(nullptr), normals(nullptr), faces(nullptr), num_vertices(0), num_normals(0), num_faces(0) {
            myFile.open(filename);
            if (!myFile.is_open()) {
                cerr << "Error opening file: " << filename << endl;
                return;
            }
            string line;
            while (getline(myFile, line)) {
                if (line[0] == 'v' && line[1] == ' ') {
                    num_vertices++;
                } else if (line[0] == 'v' && line[1] == 'n') {
                    num_normals++;
                } else if (line[0] == 'f') {
                    num_faces++;
                }
            }
            vertices = new GLfloat[3 * num_vertices];
            normals = new GLfloat[3 * num_normals];
            faces = new GLuint[9 * num_faces];
            center = new float[3]{0.0f, 0.0f, 0.0f};
            myFile.clear();
            myFile.seekg(0, ios::beg);
            long long v_index = 0, n_index = 0, f_index = 0;
            while (getline(myFile, line)) {
                if (line[0] == 'v' && line[1] == ' ') {
                    get_vertices(line, vertices, v_index++);
                } else if (line[0] == 'v' && line[1] == 'n') {
                    get_normals(line, normals, n_index++);
                } else if (line[0] == 'f') {
                    get_faces(line, faces, f_index++);
                }
            }

            myFile.close();

            for(long long i = 0; i < num_vertices; i++) {
                center[0] += vertices[3 * i];
                center[1] += vertices[3 * i + 1];
                center[2] += vertices[3 * i + 2];
            }
            center[0] /= num_vertices;
            center[1] /= num_vertices;
            center[2] /= num_vertices;

            indices = new GLuint[3 * num_faces];
            for (long long i = 0; i < num_faces; i++) {
                indices[3 * i] = faces[9 * i] - 1;
                indices[3 * i + 1] = faces[9 * i + 3] - 1;
                indices[3 * i + 2] = faces[9 * i + 6] - 1;
            }
        }

        void translate(float dx, float dy, float dz) {
            for (long long i = 0; i < num_vertices; i++) {
                vertices[3 * i] += dx;
                vertices[3 * i + 1] += dy;
                vertices[3 * i + 2] += dz;
            }
            center[0] += dx;
            center[1] += dy;
            center[2] += dz;
        }

        void rotate(float theta, float phi, float psi) {
            glm::mat3 R = glm::mat3(glm::rotate(glm::mat4(1.0f), glm::radians(theta), glm::vec3(1.0f, 0.0f, 0.0f)) *
                                    glm::rotate(glm::mat4(1.0f), glm::radians(phi), glm::vec3(0.0f, 1.0f, 0.0f)) *
                                    glm::rotate(glm::mat4(1.0f), glm::radians(psi), glm::vec3(0.0f, 0.0f, 1.0f)));
            for (long long i = 0; i < num_vertices; i++) {
                glm::vec3 v = glm::vec3(vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
                v = R * v;
                vertices[3 * i] = v.x;
                vertices[3 * i + 1] = v.y;
                vertices[3 * i + 2] = v.z;
            }
        }
        
};

class scene {
    public:
        object3D** objects;
        int num_objects;
        scene() : objects(nullptr), num_objects(0) {}
        
        void add_object(const string& filename) {
            object3D** new_objects = new object3D*[num_objects + 1];
            for (int i = 0; i < num_objects; i++) {
                new_objects[i] = objects[i];
            }
            new_objects[num_objects] = &object3D(filename);
            delete[] objects;
            objects = new_objects;
            num_objects++;
        }

        void build_scene() {
            for (int i = 0; i < num_objects; i++) {

            }
        }
};

void get_vertices(string line, GLfloat* vertices, int num_vertices){
    string mynum;


    line = line.substr(2);

    for(int i = 0; i < 3; i++){
        int j = 0;
        while(line[j] != ' ' && j < line.size()){
            mynum += line[j];
            j++;
        }
        vertices[3*num_vertices + i] = stof(mynum);
        mynum = "";
        if (j < line.size()) line = line.substr(j+1);
    }
}


void get_normals(string line, GLfloat* normals, int num_normals){
    string mynum;

    line = line.substr(3);

    for(int i = 0; i < 3; i++){
        int j = 0;
        while(line[j] != ' ' && j < line.size()){
            mynum += line[j];
            j++;
        }
        normals[3*num_normals + i] = stof(mynum);
        mynum = "";
        if (j<line.size()) line = line.substr(j+1);
    }
}


void get_faces(string line , GLuint* faces, int num_faces){
    string mynum;

    line = line.substr(2);
    for(int i = 0; i < 9; i++){
        int j = 0;
        while(line[j] != ' ' && j < line.size() && line[j] != '/'){
            mynum += line[j];
            j++;
        }
        faces[9 * num_faces + i] = stoi(mynum);
        mynum = "0";
        if(j < line.size()) line = line.substr(j+1);
    }
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}  


int main(int argc, char** argv)
{
    
    
    GLFWwindow* window;
    
    /* Initialize the library */
    if (!glfwInit())
    {
        return -1;
    }
    
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    
    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(640, 480, "Scene", nullptr, nullptr);
    
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    
    /* Make the window's context current */
    glfwMakeContextCurrent(window);
    
    glewExperimental = GL_TRUE;
    
    // Initialize GLEW
    if (glewInit() != GLEW_OK)
    {
        std::cerr << "ERROR: GLEW Initialization Failed\n";
        return -1;
    }
    
    // OpenGL begins here

    

    long long num_vertices = 0;
    long long num_normals = 0;
    long long num_faces = 0;
    
    // Getting data from the file
    GLfloat* vertices;
    GLfloat* normals;
    GLuint* faces;
    
    string line;
    ifstream myFile(argv[1]);

    while(getline(myFile, line)){
        if(line[0] == 'v' && line[1] == ' '){
            num_vertices++;
        }
        else if(line[0] == 'v' && line[1] == 'n'){
            num_normals++;
        }
        else if(line[0] == 'f'){
            num_faces++;
        }
    }
    vertices = new GLfloat[3*num_vertices];
    normals = new GLfloat[3*num_normals];
    faces = new GLuint[9*num_faces];
    myFile.close();
    
    myFile.open(argv[1]);
    
    num_vertices = 0;
    num_normals = 0;
    num_faces = 0;
    
    while(getline(myFile, line)){
        if(line[0] == 'v' && line[1] == ' '){
            get_vertices(line, vertices, num_vertices);
            num_vertices++;
        }
        else if(line[0] == 'v' && line[1] == 'n'){
            get_normals(line, normals, num_normals);
            num_normals++;
        }
        else if(line[0] == 'f'){
            get_faces(line, faces, num_faces);
            num_faces++;
        }
    }
    
    myFile.close();

    GLuint* indices = new GLuint[3*num_faces];
    for(int i = 0; i < num_faces; i++){
        indices[3*i] = faces[9*i]-1;
        indices[3*i + 1] = faces[9*i + 3]-1;
        indices[3*i + 2] = faces[9*i + 6]-1;
    }   

    //Rotating vertices matrix arg[2], arg[3] and arg[4] degrees
    float theta = stof(argv[2]);
    float phi = stof(argv[3]);
    float psi = stof(argv[4]);

    glm::mat3 R = glm::mat3(glm::rotate(glm::mat4(1.0f), glm::radians(theta), glm::vec3(1.0f, 0.0f, 0.0f)) *
                           glm::rotate(glm::mat4(1.0f), glm::radians(phi), glm::vec3(0.0f, 1.0f, 0.0f)) *
                           glm::rotate(glm::mat4(1.0f), glm::radians(psi), glm::vec3(0.0f, 0.0f, 1.0f)));

    for(int i = 0; i < num_vertices; i++){
        glm::vec3 v = glm::vec3(vertices[3*i], vertices[3*i + 1], vertices[3*i + 2]);
        v = R*v;
        vertices[3*i] = 4* v.x;
        vertices[3*i + 1] = 4* v.y;
        vertices[3*i + 2] = 4* v.z;
    }

    // Create VAO, VBO and EBO
    GLuint VAO; //Stores the configuration of how vertex data is read from the VBO.
    GLuint VBO; //Stores vertex data (positions, colors, texture coordinates, etc.) in GPU memory.
    GLuint EBO; //Stores indices in GPU memory.
    
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &EBO);
    glGenBuffers(1, &VBO);
    
    glBindVertexArray(VAO); //make it active
    
    glBindBuffer(GL_ARRAY_BUFFER, VBO); //make it active
    glBufferData(GL_ARRAY_BUFFER, 3 * num_vertices * sizeof(vertices[0]), vertices, GL_STATIC_DRAW);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*num_faces*sizeof(indices[0]), indices, GL_STATIC_DRAW);

    //Defines how vertex data is structured
    glVertexAttribPointer(0,// Attribute location in the shader
                          3,//Number of components per vertex
                          GL_FLOAT,//Data type of each component
                          GL_FALSE,//Normalize integer values to [0,1]
                          0,//Byte offset between consecutive attributes
                          nullptr);//Offset from the start of VBO where data begins;
    glEnableVertexAttribArray(0); //Enables the attribute (position)

        

    // Shaders
    const char * vertex_shader =
        "#version 330 core\n"
        "in vec3 vp;"
        "varying float zDepth;"
        "void main(){"
        "gl_Position = vec4(vp, 1.0);"
        "zDepth = (2.0 - (1.0 + gl_Position.z))/2.0;"
        "}";

    const char * fragment_shader =
        "#version 330 core\n"
        "out vec4 frag_color;"
        "varying highp float zDepth;"
        "void main(){"
        "frag_color = vec4(zDepth, 0.0, 0.0, 1.0);"
        "}";


    // Compile and link the shaders
    int  success;
    char infoLog[512];
        
    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &vertex_shader, nullptr);
    glCompileShader(vs);
    glGetShaderiv(vs, GL_COMPILE_STATUS, &success);

    if(!success){
        glGetShaderInfoLog(vs, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &fragment_shader, nullptr);
    glCompileShader(fs);
    glGetShaderiv(vs, GL_COMPILE_STATUS, &success);

    if(!success){
        glGetShaderInfoLog(fs, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, fs);
    glAttachShader(shaderProgram, vs);
    glLinkProgram(shaderProgram);
    
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if(!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::LINKING_FAILED\n" << infoLog << std::endl;
    }
    
    glDeleteShader(vs);
    glDeleteShader(fs);
    
    // Set up the viewport
    glViewport(0, 0, 800, 600);
    
    
    // Wireframe mode
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glEnable(GL_DEPTH_TEST);


    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
        /* Render here */
        //glClearColor(0.0f, 0.0f, 0.5f, 0.4f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);
        glBufferData(GL_ARRAY_BUFFER, 3 * num_vertices * sizeof(vertices[0]), vertices, GL_STATIC_DRAW);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * num_faces*sizeof(indices[0]), indices, GL_STATIC_DRAW);
        glDrawElements(GL_TRIANGLES, 3*num_faces, GL_UNSIGNED_INT, 0);

        glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
        
        /* Swap front and back buffers */
        glfwSwapBuffers(window);
        
        /* Poll for and process events */
        glfwPollEvents();

        
        float dely = 0;
        float delx = 0;
        float delz = 0;
        float speed = 0.1;
        
        if(glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS){
            delx = -speed;
        }else if(glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS){
            delx = speed;
        }
        if(glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS){
            delz = speed;
        }else if(glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS){
            delz = -speed;
        } 
        if(glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS){
            dely = speed;
        }else if(glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS){
            dely = -speed;
        }

        char rotx = 'c';
        char roty = 'c';
        if(glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS){
            rotx = 'a';
        }else if(glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS){
            rotx = 'd';
        }
        if(glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS){
            roty = 'w';
        }else if(glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS){
            roty = 's';
        }

        theta = 0;
        phi = 0;
        psi = 0;
        
        if (rotx == 'a'){
            phi = 15;
        }else if(rotx == 'd'){
        phi = -15;
    }
        if(roty == 'w'){
        theta = 15;
    }else if(roty == 's'){
        theta = -15;
    }
    
    
        //making the camera rotate

        

        R = glm::mat3(glm::rotate(glm::mat4(1.0f), glm::radians(theta), glm::vec3(1.0f, 0.0f, 0.0f)) *
        glm::rotate(glm::mat4(1.0f), glm::radians(phi), glm::vec3(0.0f, 1.0f, 0.0f)) *
        glm::rotate(glm::mat4(1.0f), glm::radians(psi), glm::vec3(0.0f, 0.0f, 1.0f)));
        
        for(int i = 0; i < num_vertices; i++){
            glm::vec3 v = glm::vec3(vertices[3*i], vertices[3*i + 1], vertices[3*i + 2]);
            v = R*v;
            vertices[3*i] = v.x;
            vertices[3*i + 1] = v.y;
            vertices[3*i + 2] = v.z;
        }

           
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
           
           
           //Maknig the camera walk
        for(int i = 0; i < num_vertices; i++){
            glm::vec3 v = glm::vec3(vertices[3*i], vertices[3*i + 1], vertices[3*i + 2]);
            v.x += delx;
            v.y += dely;
            v.z += delz;
            vertices[3*i] = v.x;
            vertices[3*i + 1] = v.y;
            vertices[3*i + 2] = v.z;
        }
        
    }

    glfwTerminate();
    
    return 0;
}