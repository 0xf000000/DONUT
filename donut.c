#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// this Programm will be the ground skelleton for a Donut printing progranmm

int printDonut(float A, float B); 

int main(int argc, char** argv){

    float A = 0, B =0;
    while(1){
        printDonut(A,B);

        A -= 0.0009;
        B -= 0.0005;

    }

    // 1920 x 1080 
  

    return 0;
}

 const int screenWidth = 50; 
    const int screenHeight = 50; 
   
    const float theta_spacing= 0.07; 
    const float phi_spacing= 0.02;
    const float R1= 1; 
    const float R2 =2; 
    const float K2 = 5;
    const float K1 = screenWidth*K2*3/(8*(R1+R2));
 

int printDonut(float A, float B){

   


    float cosA = cos(A), sinA = sin(A); 
    float cosB = cos(B), sinB = sin(B);

     char** output = (char**) malloc(screenWidth * sizeof(char*));

     for (int i = 0; i < screenWidth; i++){
        output[i] = (char*)malloc(screenHeight * sizeof(char));
     }
       
        float** zbuffer = (float**) malloc(screenWidth * sizeof(float*));

     for (int i = 0; i < screenWidth; i++){
        zbuffer[i] = (float*)malloc(screenHeight * sizeof(float));
     }

    // fill array with ' ' 
    

    for(int i = 0; i < screenWidth; i++){
        for(int j = 0;  j < screenHeight; j++){
            output[i][j] = ' ';
        }
    }

    // fill zbuffer with 0 
    
    for(int i = 0; i < screenWidth; i++){
        for(int j = 0;  j < screenHeight ; j++){
            zbuffer[i][j] = 0;
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
        int yp = (int) (screenHeight/2 - K1 * ooz *y); 

        // calculate luminace so this will decide which pixel we will see
        float L = cosPhi * cosTheta*sinB - cosA*cosTheta*sinPhi - sinA * sinTheta + cosB * (cosA* sinTheta - cosTheta * sinA * sinPhi);


        // the range of the luminance will be -sqrt(2) to sqrt(2) so if its  < 0 we will not bother to make it visible

        if(L > 0  ){

           
            
        
            if(ooz > zbuffer[xp][yp]){
                
                zbuffer[xp][yp] = ooz;
                int luminance_index = L*8;
                output[xp][yp] = ".,-~:;=!*#$@"[luminance_index]; 

            }
       

        }
    }
}

printf("\x1b[H");
for(int j = 0; j< screenHeight; j++ ){
    for(int i=0; i < screenWidth; i++){

       
        putchar(output[i][j]);
    }
    putchar('\n');
}

return 0;
}


