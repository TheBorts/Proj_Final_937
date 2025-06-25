#ifndef CAMERA_H
#define CAMERA_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

enum camera_Movement {
    FORWARD,
    BACKWARD,
    LEFT,
    RIGHT
};

const float YAW = -180.0f;
const float PITCH = 0.0f;
const float SPEED = 90.0f;
const float SENSITIVITY = 0.6f;
const float ZOOM = 45.0f;

class camera{
    public:
        glm::vec3 Position;
        glm::vec3 Front;
        glm::vec3 Up;
        glm::vec3 Right;
        glm::vec3 WorldUp;

        float Yaw;
        float Pitch;
        
        float MoveSpeed;
        float MouseSense;
        float Zoom;

        // Default constructor
        camera(glm::vec3 position = glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f), float yaw = YAW, float pitch = PITCH) : MoveSpeed(SPEED), MouseSense(SENSITIVITY), Zoom(ZOOM){
            Position = position;
            WorldUp = up;
            Yaw = yaw;
            Pitch = pitch;
            Front = -1.0f * glm::normalize(Position); // Camera looking at the origin
            updateCameraVectors();
        }

        // Constructor with parameters for position, up vector, yaw, and pitch and camera looking at the origin
        camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) : MoveSpeed(SPEED), MouseSense(SENSITIVITY), Zoom(ZOOM){
            Position = glm::vec3(posX, posY, posZ);
            WorldUp = glm::vec3(upX, upY, upZ);
            Yaw = yaw;
            Pitch = pitch;
            Front = -1.0f * glm::normalize(Position); // Camera looking at the origin
            updateCameraVectors();
        }

        // Returns the view matrix
        glm::mat4 getViewMatrix(){
            return glm::lookAt(Position, Position + Front, Up);
        }

        // Processes input from mouse movement
        void ProcessMouseMove(float xoffset, float yoffset, GLboolean constrainPitch = true) {
            xoffset *= MouseSense;
            yoffset *= MouseSense;

            Yaw += xoffset;
            Pitch += yoffset;

            if (constrainPitch) {
                if (Pitch > 89.0f)
                    Pitch = 89.0f;
                if (Pitch < -89.0f)
                    Pitch = -89.0f;
            }

            updateCameraVectors();
        }

        // Processes input from mouse scroll
        void ProcessKeyboard(camera_Movement direction, float deltaTime){
            float velocity = MoveSpeed * deltaTime;
            if (direction == FORWARD)
                Position += Front * velocity;
            if (direction == BACKWARD)
                Position -= Front * velocity;
            if (direction == LEFT)
                Position -= Right * velocity;
            if (direction == RIGHT)
                Position += Right * velocity;

            Position.y = 1.0f; // Prevent camera from going below ground level
        }

    private:
        
        // Updates the camera vectors based on the current yaw and pitch
        void updateCameraVectors(){
            glm::vec3 front;
            front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
            front.y = sin(glm::radians(Pitch));
            front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
            Front = glm::normalize(front);

            Right = glm::normalize(glm::cross(Front, WorldUp));
            Up = glm::normalize(glm::cross(Right, Front));
        }
};


#endif
