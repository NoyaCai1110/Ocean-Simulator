#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

#include "shader.h"
 
#include <iostream>
#include <cmath>
#include "CImg.h"
using namespace std;
using namespace cimg_library;

#define WindowSizeX 1000
#define WindowSizeY 1000
//#define _USE_MATH_DEFINES 
#define E     2.71828182845904523536

#pragma pack(2)

typedef unsigned char  BYTE;
typedef unsigned short WORD;
typedef unsigned long  DWORD;
typedef long    LONG;


//BMP File Header
struct MYBITMAPFILEHEADER
{
    WORD  bfType;		
    DWORD bfSize;		
    WORD  bfReserved1;	
    WORD  bfReserved2;	
    DWORD bfOffBits;	
};

//BMP Info Header
struct MYBITMAPINFOHEADER
{
    DWORD biSize;			
    LONG  biWidth;			
    LONG  biHeight;			
    WORD  biPlanes;			
    WORD  biBitCount;		
    DWORD biCompression;    
    DWORD biSizeImage;      
    LONG  biXPelsPerMeter;  
    LONG  biYPelsPerMeter;  
    DWORD biClrUsed;        
    DWORD biClrImportant;	
};

//Pixel color info
struct RGBColor
{
    char B;		
    char G;		
    char R;		
};

void WriteBMP(const char* FileName, RGBColor* ColorBuffer, int ImageWidth, int ImageHeight)
{
    
    const int ColorBufferSize = ImageHeight * ImageWidth * sizeof(RGBColor);

    MYBITMAPFILEHEADER fileHeader;
    fileHeader.bfType = 0x4D42;	
    fileHeader.bfReserved1 = 0;
    fileHeader.bfReserved2 = 0;
    fileHeader.bfSize = sizeof(MYBITMAPFILEHEADER) + sizeof(MYBITMAPINFOHEADER) + ColorBufferSize;
    fileHeader.bfOffBits = sizeof(MYBITMAPFILEHEADER) + sizeof(MYBITMAPINFOHEADER);

    MYBITMAPINFOHEADER bitmapHeader = { 0 };
    bitmapHeader.biSize = sizeof(MYBITMAPINFOHEADER);
    bitmapHeader.biHeight = ImageHeight;
    bitmapHeader.biWidth = ImageWidth;
    bitmapHeader.biPlanes = 1;
    bitmapHeader.biBitCount = 24;
    bitmapHeader.biSizeImage = ColorBufferSize;
    bitmapHeader.biCompression = 0; 


    FILE* fp;

    fopen_s(&fp, FileName, "wb");

    fwrite(&fileHeader, sizeof(MYBITMAPFILEHEADER), 1, fp);
    fwrite(&bitmapHeader, sizeof(MYBITMAPINFOHEADER), 1, fp);

    fwrite(ColorBuffer, ColorBufferSize, 1, fp);

    fclose(fp);
}

void ScreenShot()
{
    RGBColor* ColorBuffer = new RGBColor[WindowSizeX * WindowSizeY];

    glReadPixels(0, 0, WindowSizeX, WindowSizeY, GL_BGR, GL_UNSIGNED_BYTE, ColorBuffer);

    WriteBMP("./output/Water.bmp", ColorBuffer, WindowSizeX, WindowSizeY);
    
    delete[] ColorBuffer;
}

void tone_reproduction(){
    CImg<unsigned char> img("./output/GI.bmp");
    CImg<unsigned char> TRImgWard[3] = {img, img, img};
    CImg<unsigned char> TRImgReinhard[3] = {img, img, img};
    double Ldmax = 500;
    double Lmax[3] = {10, 100, 1000};
    double Lwa[3] = {0, 0, 0};
    double Lxy[3] = {0, 0, 0};
    double sf[3] = {0, 0, 0};
    double test = 0;
    cimg_forXY(img, x, y)
    {
        for(int i = 0; i < 3; i ++)
        {
            Lxy[i] = (0.27 * (double)img(x, y, 0) + 0.67 * (double)img(x, y, 1) + 0.06 * (double)img(x, y, 2))* Lmax[i]  / (255);
            Lwa[i] += log(0.000000000001 + Lxy[i])/((double)img._height * (double)img._width);
        }      
        if(img(x, y, 0) > test) test = img(x, y, 0);
    }
    //cout << test << endl;
    for(int i = 0; i < 3; i ++){
        cout << "Lwa = " << Lwa[i] << endl;
        Lwa[i] = exp(Lwa[i]);
        //cout << "Lwa = " << Lwa[i] << endl;
        sf[i] = pow((1.219 + pow((Ldmax/2), 0.4))/(1.219 + pow(Lwa[i], 0.4)), 2.5);
        //cout << "sf = " << sf[i] << endl;
        cimg_forXY(TRImgWard[i], x, y){
            for(int j = 0; j < 3; j ++){
                if(Lmax[i] * sf[i] * TRImgWard[i](x, y, j) / Ldmax > 255)
                    TRImgWard[i](x, y, j) = 255;
                else
                    TRImgWard[i](x, y, j) = Lmax[i] * sf[i] * TRImgWard[i](x, y, j) / Ldmax;         
            }        
        }
    }
    // choose the keyvalue that you want to use for reinhard tone reproduction
    double a = 0.18;
    //double a = 0.09;
    //double a = 0.36;
    //double a = 0.72;
    for(int i = 0; i < 3; i ++){
        cimg_forXY(TRImgReinhard[i], x, y){
            double RGB_scaled[3] = {0, 0, 0}; 
            double RGB_target[3] = {0, 0, 0};
            for(int j = 0; j < 3; j ++){
                RGB_scaled[j] = Lmax[i] * a * TRImgReinhard[i](x, y, j) / (Lwa[i] * 255);
                RGB_target[j] = RGB_scaled[j] / (RGB_scaled[j] + 1);
                if(RGB_target[j] > 1)
                    TRImgReinhard[i](x, y, j) = 255;
                else
                    TRImgReinhard[i](x, y, j) = RGB_target[j] * 255;         
            }        
        }
    }
    TRImgWard[0].save("./output/TR_Ward_Low.bmp");
    TRImgWard[1].save("./output/TR_Ward_Mid.bmp");
    TRImgWard[2].save("./output/TR_Ward_High.bmp");
    TRImgReinhard[0].save("./output/TR_Reinhard_Low.bmp");
    TRImgReinhard[1].save("./output/TR_Reinhard_Mid.bmp");
    TRImgReinhard[2].save("./output/TR_Reinhard_High.bmp"); 
}

void tone_reproduction_Advanced(){  //Adaptive Logarithmic Mapping
    CImg<unsigned char> img("./output/GI.bmp");
    CImg<unsigned char> TR_ALM[3] = {img, img, img};
    double Ldmax = 500;
    double Lmax[3] = {10, 100, 1000};
    double Lwmax[3] = {0, 0, 0};
    double Lwa[3] = {0, 0, 0};
    double Lxy[3] = {0, 0, 0};

    cimg_forXY(img, x, y)
    {
        for(int i = 0; i < 3; i ++)
        {
            Lxy[i] = (0.27 * (double)img(x, y, 0) + 0.67 * (double)img(x, y, 1) + 0.06 * (double)img(x, y, 2))* Lmax[i]  / (255);
            Lwa[i] += log(0.000000000001 + Lxy[i])/((double)img._height * (double)img._width);
        }      
    }
    double b = 0.85;
    for(int i = 0; i < 3; i ++){
        //cout << "Lwa = " << Lwa[i] << endl;
        Lwa[i] = exp(Lwa[i]);
        //cout << "Lwa = " << Lwa[i] << endl;
        Lwmax[i] = Lmax[i] / Lwa[i];
        cimg_forXY(TR_ALM[i], x, y){
            double Lw = (0.27 * (double)TR_ALM[i](x, y, 0) + 0.67 * (double)TR_ALM[i](x, y, 1) + 0.06 * (double)TR_ALM[i](x, y, 2))* Lmax[i]  / (255 * Lwa[i]);
            double Ld = log(Lw + 1) / (log10(Lwmax[i] + 1) * log(2 + 8 * (pow(Lw/Lwmax[i], log(b)/log(0.5)))));
            for(int j = 0; j < 3; j ++){
                if(Ld * TR_ALM[i](x, y, j) > 255){
                    TR_ALM[i](x, y, j) = 255;
                }
                else{
                    TR_ALM[i](x, y, j) = Ld * TR_ALM[i](x, y, j);
                }
            }
        }
    }
    
    TR_ALM[0].save("./output/TR_ALM_Low.bmp");
    TR_ALM[1].save("./output/TR_ALM_Mid.bmp");
    TR_ALM[2].save("./output/TR_ALM_High.bmp");
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

Shader * shaderPtr;
double mouseX = 0, mouseY = 0;
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {

        glfwGetCursorPos(window, &mouseX, &mouseY);
        mouseX /= WindowSizeX;
        mouseY /= WindowSizeY;
        mouseY = 1 - mouseY;//in shader, coordinate starts at bottom left, while mouse coordinate starts at top left
        cout << "Left pressed, " << mouseX << " " << mouseY << endl;

        float time = glfwGetTime();
        shaderPtr->setFloat("timeStamp", time);
    }
}

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(1000, 1000, "RayTracer", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    // glad: load all OpenGL function pointers
    // ---------------------------------------
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
 
    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    
    
    /*
    float vertices[] = {
        // positions         // colors
        -2.0f, -1.0f, -1.0f,  1.0f, 0.0f, 0.0f,  // bottom right
        -2.0f, -1.0f, -4.5f,  0.0f, 1.0f, 0.0f,  // bottom left
        0.9f, -1.0f, -1.0f,   1.0f, 1.0f, 1.0f ,   // top 
        0.9f, -1.0f, -4.5f,   0.0f, 0.0f, 1.0f
    };
    */
    ///*
    float vertices[] = {
			 -1.0f,  -1.0f, 0.0f,  1.0f, 0.0f, 0.0f,// top right
			 1.0f, -1.0f, 0.0f,  0.0f, 1.0f, 0.0f,// bottom right
			-1.0f,  1.0f, 0.0f,  1.0f, 1.0f, 1.0f,// bottom left
			1.0f,  1.0f, 0.0f,   0.0f, 0.0f, 1.0f// top left 
		};
    //*/
    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);
 
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
 
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
 
    glBindVertexArray(0);
 
    //Shader ColorShader("./src/shaders/GI.vs", "./src/shaders/GI.fs");
    Shader WaterShader("./src/shaders/Water.vert", "./src/shaders/Water.frag");
    shaderPtr = &WaterShader;
    glm::mat4 model = glm::mat4(1.0f);;
    glm::mat4 view = glm::mat4(1.0f);;
    glm::mat4 projection = glm::mat4(1.0f);;
    glm::mat4 mvp;
    
    /*mvp = projection * view * model;
    cout << "MVP" << endl;
    for(int i = 0; i < 4; i ++){
        for(int j = 0; j < 4; j ++){
            cout << mvp[i][j] << " ";
        }
        cout << " " <<endl;
    }*/
    
    // render loop
    // -----------
    bool saved = false;
    float time_t = 0;
    float now = 0;
    float delta = 0;
    float previous = 0;
    while (!glfwWindowShouldClose(window))
    {
 
        // render
        // ------
        glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
        //glClearColor(0.0, 0.6, 0.8, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// Enable depth testing since it's disabled when drawing the framebuffer rectangle
		glDisable(GL_DEPTH_TEST); // prevents framebuffer rectangle from being discarded

        /*
        int modelLoc = glGetUniformLocation(ColorShader.ID, "model");
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
        int viewLoc = glGetUniformLocation(ColorShader.ID, "view");
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
        int projectionLoc = glGetUniformLocation(ColorShader.ID, "projection");
        glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
        int mvpLoc = glGetUniformLocation(ColorShader.ID, "mvp");
        glUniformMatrix4fv(mvpLoc, 1, GL_FALSE, glm::value_ptr(mvp));
        */
        ///*
        now = glfwGetTime();
        delta = now - previous;
        previous = now;
        time_t += delta;
        //time_t = time_t - floor(time_t); 
        int timeLoc = glGetUniformLocation(WaterShader.ID, "time_t");
        glUniform1f(timeLoc, time_t);
        WaterShader.setVec3("mousePos", mouseX, mouseY, -1);
        //*/
        //ColorShader.use();
        WaterShader.use();
        
        // render the triangle
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glfwSetMouseButtonCallback(window, mouse_button_callback);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    ScreenShot();
    //Uncomment the lines below to enable tone reproduction

    //tone_reproduction_Advanced();
    //tone_reproduction();

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
 
    glfwTerminate();
    return 0;
}

