#ifndef __WOLFF___RENDERING_HPP__
#define __WOLFF___RENDERING_HPP__

/*
The rendering code is largely instired by the modern OpenGL tutorial https://learnopengl.com/
*/
#include <string>
#include <vector>

#include "glad/glad.h"
#include <GLFW/glfw3.h>

struct Vector2f
{
    float x;
    float y;
};

struct Vector4f
{
    float x;
    float y;
    float z;
    float w;
};

class Renderer
{
private:
    unsigned int _window_size;
    GLFWwindow *_window;

    GLuint _vbo;
    GLuint _vbo_instances;
    GLuint _vao;
    GLuint _ebo;
    GLuint _shader_program;

    std::vector<Vector4f> _instance_data;

    const std::string _vs_source = R"""(
        #version 330 core
        layout(location = 0) in vec2 vertex;
        layout(location = 1) in vec4 transform;
        layout(location = 2) in vec4 color;

        uniform vec4 projection;

        out vec4 vertex_color;

        void main()
        {
            vec2 position = vec2(vertex.x * transform.z + transform.x, vertex.y * transform.w + transform.y);
            gl_Position = vec4(
                position.x * 2.0 / (projection.z - projection.x) - (projection.z + projection.x) / (projection.z - projection.x), 
                position.y * 2.0 / (projection.w - projection.y) - (projection.w + projection.y) / (projection.w - projection.y), 
                0.0, 
                1.0
            );

            vertex_color = color;
        }
    )""";

    const std::string _fs_source = R"""(
        #version 330 core
        out vec4 FragColor;

        in vec4 vertex_color;

        void main()
        {
            FragColor = vertex_color;
        } 
    )""";

    int create_window(const std::string title)
    {
        int init_result = glfwInit();
        if (!init_result)
        {
            std::cerr << "Failed to initialize GLFW\n";

            glfwTerminate();
            return 0;
        }

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

        _window = glfwCreateWindow(_window_size, _window_size, title.c_str(), nullptr, nullptr);
        if (!_window)
        {
            std::cerr << "Unable to create GLFW window\n";

            glfwTerminate();
            return 0;
        }

        glfwMakeContextCurrent(_window);
        // glfwSwapInterval(0);

        int gladInitRes = gladLoadGL();
        if (!gladInitRes)
        {
            std::cerr << "Unable to initialize glad\n";
            return 0;
        }

        return 1;
    }

    void create_buffers()
    {
        static float vertices[] = {
            1.0f, 1.0f, // top right
            1.0f, 0.0f, // bottom right
            0.0f, 0.0f, // bottom left
            0.0f, 1.0f, // top left
        };
        static unsigned int indices[] = {
            0, 1, 3, // first triangle
            1, 2, 3  // second triangle
        };

        // VAO
        glGenVertexArrays(1, &_vao);
        glBindVertexArray(_vao);

        // VBO
        glGenBuffers(1, &_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, _vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

        // VBO instances
        glGenBuffers(1, &_vbo_instances);

        // EBO
        glGenBuffers(1, &_ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, _vbo);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void *)0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, _vbo_instances);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 2 * sizeof(Vector4f), (void *)0);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 2 * sizeof(Vector4f), (void *)(1 * sizeof(Vector4f)));
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glVertexAttribDivisor(1, 1);
        glVertexAttribDivisor(2, 1);

        glBindVertexArray(0);
    }

    void create_shaders()
    {
        GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
        const char *vs_code = _vs_source.c_str();
        glShaderSource(vertex_shader, 1, &vs_code, NULL);
        glCompileShader(vertex_shader);
        check_shader_error(vertex_shader);

        GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
        const char *fs_code = _fs_source.c_str();
        glShaderSource(fragment_shader, 1, &fs_code, NULL);
        glCompileShader(fragment_shader);
        check_shader_error(fragment_shader);

        _shader_program = glCreateProgram();
        glAttachShader(_shader_program, vertex_shader);
        glAttachShader(_shader_program, fragment_shader);
        glLinkProgram(_shader_program);

        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);
    }

    int check_shader_error(GLuint shader) const
    {
        int success;
        char infoLog[512];
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

        if (!success)
        {
            glGetShaderInfoLog(shader, 512, NULL, infoLog);
            std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n"
                      << infoLog << std::endl;
        }
        return success;
    }

public:
    Renderer() = default;

    void create(unsigned int window_size, std::string title)
    {
        _window_size = window_size;

        if (!create_window(title))
        {
            exit(EXIT_FAILURE);
        }
        create_buffers();
        create_shaders();
    }

    void set_projection(const Vector4f &projection) const
    {
        glUseProgram(_shader_program);
        glUniform4f(glGetUniformLocation(_shader_program, "projection"), projection.x, projection.y, projection.z, projection.w);
    }

    void draw_rect(const Vector2f &position, const Vector2f &scale, const Vector4f &color)
    {
        _instance_data.push_back({position.x, position.y, scale.x, scale.y});
        _instance_data.push_back(color);
    }

    void end_draw()
    {
        glBindBuffer(GL_ARRAY_BUFFER, _vbo_instances);
        glBufferData(GL_ARRAY_BUFFER, sizeof(Vector4f) * _instance_data.size(), &_instance_data[0], GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glUseProgram(_shader_program);

        glBindVertexArray(_vao);
        // glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, _instance_data.size() / 2);
        glBindVertexArray(0);

        glfwSwapBuffers(_window);
        _instance_data.clear();
    }

    bool is_running() const
    {
        glfwPollEvents();
        return !glfwWindowShouldClose(_window);
    }
};

#endif