#version 330 core
precision highp float;
out vec4 FragColor;

struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;    
    float shininess;
    int hasTexture;
}; 

struct Light {
    vec3 position;
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord; 

uniform sampler2D OurTexture;
uniform vec3 viewPos;
uniform Material material;
uniform Light light;

void main()
{
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(light.position - FragPos);
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);  
    
    float diff = max(-dot(norm, lightDir), 0.0);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);

    vec3 ambient, diffuse;
    if (material.hasTexture == 1) {
        ambient = light.ambient * texture(OurTexture, TexCoord).rgb;
        diffuse = light.diffuse * (diff * texture(OurTexture, TexCoord).rgb);
    } else {
        ambient = light.ambient * material.ambient;
        diffuse = light.diffuse * (diff * material.diffuse);
    }
    
    vec3 specular = light.specular * (spec * material.specular); 
    vec3 result = ambient + diffuse + specular;

    FragColor = vec4(result, 1.0);    
}