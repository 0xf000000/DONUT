#include <stdio.h>
#include <Math.h>

// this Programm will be the ground skelleton for a Donut printing progranmm

void renderFrame(float A, float B); 

int main(void){
 


    // 1920 x 1080 
   renderFrame(0,0);

    return 0;
}

void renderFrame(float A, float B){
    const int screenWidth = 1920; 
    const int screenHeight = 1080; 
    const float theta_spacing= 0.07; 
    const float phi_spacing= 0.02;
    const float R1= 1; 
    const float R2 =2; 
    const float K2 = 5;
    const float K1 = 1920 *K2*3/((R1 + R2));
    const char characters[12] = ".,-~:;=!*#$@";



    float cosA = cos(A), sinA = sin(A); 
    float cosB = cos(B), sinB = sin(B);

    // fill array with ' ' 
    char output[screenWidth][screenHeight];

    for(int O = 0; O < screenWidth -1; O++){
        for(int Q = 0;  Q < screenHeight - 1; Q++){
            output[O][Q] = ' ';
        }
    }

    // fill zbuffer with 0 
    float zbuffer[screenWidth][screenHeight];
    for(int W = 0; W < screenWidth -1; W++){
        for(int D = 0;  D < screenHeight - 1; D++){
            zbuffer[W][D] = 0;
        }
    }






// 
for(float theta=0; theta < 2* M_PI; theta += theta_spacing){
    // cos and sin of thetha
    float cosTheta = cos(theta), sinTheta = sin(theta);
    
    for(float phi =0; phi < 2 * M_PI; phi += phi_spacing){
        // cos and sin of phi
        float cosPhi = cos(phi), sinPhi = sin(phi); 

        float circleX = R2 + R1*cosTheta;
        float circleY = R1 *sinTheta; 


        // 
        float x = circleX*(cosB*cosPhi + sinA * sinB * sinPhi) - circleY*cosA*sinB;
        float y = circleX * (sinB*cosPhi-sinA*cosB*sinPhi) + circleY*cosA*cosB; 
        float z = K2 + cosA * circleX* sinPhi + circleY * sinA;
        float ooz = 1/z; // "one over z"

        // the x and y projection on the screen
        int xp = (int) (screenWidth/2 + K1 *ooz *x); 
        int yp = (int) (screenHeight/2 - K1 * ooz *z); 

        // calculate luminace so this will decide which pixel we will see
        float L = cosPhi * cosTheta*sinB - cosA*cosTheta*sinPhi - sinA * sinTheta + cosB * (cosA* sinTheta - cosTheta * sinA * sinPhi);


        // the range of the luminance will be -sqrt(2) to sqrt(2) so if its  < 0 we will not bother to make it visible

        if(L > 0  ){

            if(ooz > zbuffer[xp][yp]){
                zbuffer[xp][yp] = ooz;
                int luminance_index = L*8;
                output[xp][yp] = characters[luminance_index]; 

            }
        }
    }
}

printf("\x1b[H");
for(int j = 0; j< screenHeight; j++ ){
    for(int i=0; j < screenWidth; i++){

        printf("?=");
        putchar(output[i][j]);
    }
    putchar('\n');
}


}


