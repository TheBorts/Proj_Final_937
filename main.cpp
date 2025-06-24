#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <filesystem>
#include <string>
#include <fstream>
#include <iostream>

#define STB_IMAGE_IMPLEMENTATION
#include "include/stb_image.h"
#include "include/shaders.h"
#include "include/camera.h"
#include "include/clas.h"
#include "include/physics.h"

using namespace std;

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

camera cam(glm::vec3(4.0f, 1.0f, 0.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;


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

    world myWorld = world();

    for (int i = 1; i < argc; i++) {
        myWorld.add_object(argv[i]);
    }
    
    lightSource light = lightSource(30.0f, 0.0f, 0.0f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 1.0f, 0.0f, 1.0f);
    
    Shader meuShader("../shaders/myshader.vs", "../shaders/myshader.fs");
    //Shader meuShader("../shaders/light.vs", "../shaders/light.fs");
    
    myWorld.objs[1]->tree->m->mod->translate(0.0f, 5.0f, 0.0f);
    myWorld.objs[1]->set_positions();
    myWorld.objs[1]->isVertStatic = std::vector<bool>(myWorld.objs[1]->position.size(), false); // Initialize all vertices as static
    myWorld.objs[1]->isVertStatic[0] = true; // Make the first vertex static
    myWorld.objs[1]->isStatic = false;
    myWorld.start_world();
    
    // Set up the viewport
    glViewport(0, 0, 800, 600);
    int width, height;

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);  // Fragments closer to the camera overwrite farther ones

    // Enable wireframe mode for debugging
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Uncomment to enable wireframe

    /* Loop until the user closes the window */
    long long frameCount = 0;
    bool saveFrames = true;
    
    int sceneIndex = 1;
    std::string path = "../scenes/scene" + std::to_string(sceneIndex);
    std::filesystem::create_directory(path);
    std::string pseudoObjPath = path + "/frame";
    std::string mtlPath = "../scenes/scene" + std::to_string(sceneIndex) + "/scene.mtl";
    writer saver(pseudoObjPath + std::to_string(frameCount) + ".obj", mtlPath);

    if (saveFrames){
        saver.writeMTL(myWorld);
    }

    std::cout << "Starting simulation..." << std::endl;

    while (!glfwWindowShouldClose(window))
    {
        if (saveFrames && frameCount < 10){
            saver.setPaths(pseudoObjPath + std::to_string(frameCount) + ".obj", mtlPath);
            saver.writeOBJ(myWorld);
        }
        frameCount++;
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;


        processInput(window);
        glfwGetFramebufferSize(window, &width, &height);

        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        meuShader.use();
        
        glm::mat4 view = cam.getViewMatrix();
        glm::mat4 projection = glm::perspective(glm::radians(cam.Zoom), (float)width / (float)height, 0.1f, 3000.0f);
        glm::mat4 model = glm::mat4(1.0f);

        meuShader.setMat4("view", view);
        meuShader.setMat4("projection", projection);
        meuShader.setMat4("model", model);
        meuShader.setVec3("viewPos", cam.Position);
        
        meuShader.setVec3("light.position", light.position[0], light.position[1], light.position[2]);
        meuShader.setVec3("light.ambient", light.ambient[0], light.ambient[1], light.ambient[2]);
        meuShader.setVec3("light.diffuse", light.diffuse[0], light.diffuse[1], light.diffuse[2]);
        meuShader.setVec3("light.specular", light.specular[0], light.specular[1], light.specular[2]);
        
        myWorld.draw_world(meuShader);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);
        
        /* Poll for and process events */
        glfwPollEvents();
        
    }

    glfwTerminate();

    //free resources

    return 0;
}