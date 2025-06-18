#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <filesystem>
#include <string>
#include <fstream>
#include <iostream>

#include "include/stb_image.h"
#include "include/shaders.h"
#include "include/camera.h"
#include "include/clas.h"

using namespace std;

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

camera cam(glm::vec3(400.0f, 0.0f, 0.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}  

void mouse_callback(GLFWwindow* window, double xposIn, double yposIn){
    
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    cam.ProcessMouseMove(xoffset, yoffset);
}

void processInput(GLFWwindow *window){
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cam.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cam.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cam.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cam.ProcessKeyboard(RIGHT, deltaTime);
}

int main(int argc, char** argv)
{
    GLFWwindow* window;
    
    /* Initialize the library */
    if (!glfwInit())
    {
        return -1;
    }
    glEnable(GL_CULL_FACE);  // Enable face culling
    glCullFace(GL_BACK);     // Cull back faces (default)
    glFrontFace(GL_CCW);     // Set counter-clockwise winding order as front-facing
    
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
    
    glfwMakeContextCurrent(window);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);        

    glewExperimental = GL_TRUE;
    
    // Initialize GLEW
    if (glewInit() != GLEW_OK)
    {
        std::cerr << "ERROR: GLEW Initialization Failed\n";
        return -1;
    }
    
    // OpenGL begins here

    scene myScene = scene();

    for (int i = 1; i < argc; i++) {
        myScene.add_mesh(argv[i]);
    }

    lightSource light = lightSource(400.0f, 0.0f, 0.0f, 0.3f, 0.2f, 0.7f, 0.8f, 0.5f, 0.1f, 1.0f, 1.0f, 1.0f);

    Shader meuShader("../shaders/myshader.vs", "../shaders/myshader.fs");
    //Shader meuShader("../shaders/light.vs", "../shaders/light.fs");
    
    myScene.start_scene();

    // Set up the viewport
    glViewport(0, 0, 800, 600);
    
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);  // Fragments closer to the camera overwrite farther ones

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {

        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);

        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


        meuShader.use();
        
        glm::mat4 view = cam.getViewMatrix();
        glm::mat4 projection = glm::perspective(glm::radians(cam.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 800.0f);
        glm::mat4 model = glm::mat4(1.0f);

        meuShader.setMat4("view", view);
        meuShader.setMat4("projection", projection);
        meuShader.setMat4("model", model);
        meuShader.setVec3("viewPos", cam.Position);
        
        meuShader.setVec3("lightPos", light.position[0], light.position[1], light.position[2]);
        meuShader.setVec3("light.ambient", light.ambient[0], light.ambient[1], light.ambient[2]);
        meuShader.setVec3("light.diffuse", light.diffuse[0], light.diffuse[1], light.diffuse[2]);
        meuShader.setVec3("light.specular", light.specular[0], light.specular[1], light.specular[2]);

        myScene.draw_scene(meuShader);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);
        
        /* Poll for and process events */
        glfwPollEvents();
        
    }

    glfwTerminate();
    
    return 0;
}