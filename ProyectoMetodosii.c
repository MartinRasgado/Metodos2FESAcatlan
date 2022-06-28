//Metodo Newton-Raphson
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>
#include <stdbool.h> 

void EvaluaFuncion(int opcion, float *vectorX, int size, float *vectorF){
    float x, y, z = 0;
    switch(opcion){
        case 1:
            x = vectorX[0];
            y = vectorX[1];
            vectorF[0] = pow(x,2)+(x*y)-10;
            vectorF[1] = y+(3*(x*(pow(y,2))))-50;
            break;
        case 2:
            x = vectorX[0];
            y = vectorX[1];
            vectorF[0] = pow(x,2)+pow(y,2)-9;
            vectorF[1] = -1*exp(x)-2*(y)-3;
            break;
        case 3:
            x = vectorX[0];
            y = vectorX[1];
            z = vectorX[2];
            vectorF[0] = 2*pow(x,2)-4*x+pow(y,2)+3*pow(z,2)+6*z+2;
            vectorF[1] = pow(x,2)+pow(y,2)-2*y+2*pow(z,2)-5;
            vectorF[2] = 3*pow(x,2)-12*x+pow(y,2)-3*pow(z,2)+8;
            break;
        case 4:
            x = vectorX[0];
            y = vectorX[1];
            z = vectorX[2];
            vectorF[0] = pow(x,2)-4*x+pow(y,2);
            vectorF[1] = pow(x,2)-x-12*y+1;
            vectorF[2] = 3*pow(x,2)-12*x+pow(y,2)-3*pow(z,2)+8;
            break;
    }
}

void ImprimirM2x2(float Matriz[][2]){
	int i,j;
	printf("\n");
	for(i=0;i< 2;i++){
		for(j=0;j<2;j++){
			printf("%.12f ",Matriz[i][j]);
		}
		printf("\n");
	}
	printf("\n \n");
}
void ImprimirM3x3(float Matriz[][3]){
	int i,j;
	printf("\n");
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			printf("%.12f ",Matriz[i][j]);
		}
		printf("\n");
	}
	printf("\n \n");
}


void Jacobiana2x2(int opcion, float *vectorX, float MatrizJacobiana[][2]){
    float x,y = 0;
    switch (opcion){
    case 1:
        x = vectorX[0];
        y = vectorX[1];
        MatrizJacobiana[0][0] = (2*x) + y;
        MatrizJacobiana[0][1] = x;
        MatrizJacobiana[1][0] = (3*(pow(y,2)));
        MatrizJacobiana[1][1] = (6*(x*y))+1;
        break;
    case 2:
        x = vectorX[0];
        y = vectorX[1];
        MatrizJacobiana[0][0] = 2*x;
        MatrizJacobiana[0][1] = 2*y;
        MatrizJacobiana[1][0] = -1*(exp(x));
        MatrizJacobiana[1][1] = -2;
        break;
    default:
        break;
    }
}
void Jacobiana3x3(int opcion, float *vectorX, float MatrizJacobiana[][3]){
    float x,y,z = 0;
    switch (opcion){
    case 3:
        x = vectorX[0];
        y = vectorX[1];
        z = vectorX[2];
        MatrizJacobiana[0][0] = (4*x)-4;
        MatrizJacobiana[0][1] = 2*y;
        MatrizJacobiana[0][2] = (6*z)+6;
        MatrizJacobiana[1][0] = 2*x;
        MatrizJacobiana[1][1] = (2*y)-2;
        MatrizJacobiana[1][2] = 4*z;
        MatrizJacobiana[2][0] = (6*x)-12;
        MatrizJacobiana[2][1] = 2*y;
        MatrizJacobiana[2][2] = -6*z;
        break;
    case 4:
        x = vectorX[0];
        y = vectorX[1];
        z = vectorX[2];
		MatrizJacobiana[0][0]=(2*x)-4;
		MatrizJacobiana[0][1]=2*y;
		MatrizJacobiana[0][2]=0;
		MatrizJacobiana[1][0]=(2*x)-1;
		MatrizJacobiana[1][1]=-12;
		MatrizJacobiana[1][2]=0;
		MatrizJacobiana[2][0]=(6*x)-12;
		MatrizJacobiana[2][1]=2*y;
		MatrizJacobiana[2][2]=-6*z;
        break;
    default:
        break;
    }
}

void Inversa2x2(float matriz[][2], float inversa[][2]){
	float determinante;
	float adjunta[2][2];
	int i, j = 0;
	adjunta[0][0] = matriz[1][1];
	adjunta[0][1] = matriz[1][0] * -1;
	adjunta[1][0] = matriz[0][1] * -1;
	adjunta[1][1] = matriz[0][0];
	
	determinante = matriz[0][0] * matriz[1][1] - (matriz[0][1] * matriz[1][0]);
	
	for(i=0; i<2; i++){
		for(j=0; j<2; j++){
			inversa[i][j]=(1/determinante) * adjunta[j][i];
		}
	}
}
void Inversa3x3(float matriz[][3], float inversa[][3]){
	float determinante;
	float adjunta[3][3];
	int i, j = 0;
	adjunta[0][0] = matriz[1][1] * matriz[2][2] - (matriz[1][2] * matriz[2][1]);
	adjunta[0][1] =  (matriz[1][0] * matriz[2][2] - (matriz[1][2] * matriz[2][0] ) ) * -1;
	adjunta[0][2] =  matriz[1][0] * matriz[2][1] - (matriz[1][1] * matriz[2][0]);
	adjunta[1][0] = (matriz[0][1] * matriz[2][2] - (matriz[0][2] * matriz[2][1] ) ) * -1;
	adjunta[1][1] = matriz[0][0] * matriz[2][2] - (matriz[0][2] * matriz[2][0]);
	adjunta[2][2] = matriz[0][0] * matriz[1][1] - (matriz[0][1] * matriz[1][0]);
	adjunta[1][2] = (matriz[0][0] * matriz[2][1] - (matriz[0][1] * matriz[2][0] ) ) * -1;
	adjunta[2][0] = matriz[0][1] * matriz[1][2] - (matriz[0][2] * matriz[1][1]);
	adjunta[2][1] = (matriz[0][0] * matriz[1][2] - (matriz[0][2] * matriz[1][0] ) ) * -1;
	
	
	determinante = ( matriz[0][0] * matriz[1][1] * matriz[2][2] ) + ( matriz[0][1] * matriz[1][2] * matriz[2][0] ) + ( matriz[1][0] * matriz[2][1] * matriz[0][2] )
	- ( (matriz[0][2] * matriz[1][1] * matriz[2][0] ) + ( matriz[0][1] * matriz[1][0] * matriz[2][2] ) + ( matriz[1][2] * matriz[2][1] * matriz[0][0] ));
	
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			inversa[i][j]=(1/determinante) * adjunta[j][i];
		}
	}
}

float MultMat2x2(float JacobianaInv[][2], float *vectorF, float *MatrizResult){
    float temporal = 0;
    int i, j, k = 0;
    for (i = 0; i < 2; i++ ){ //i para las filas de la matriz resultante    
        for (k = 0 ; k < 1 ; k++ ){ // k para las columnas de la matriz resultante
            temporal = 0 ;
            for (j = 0 ; j < 2 ; j++ ){ //j para realizar la multiplicacion de los elementos   de la matriz
                temporal += JacobianaInv[i][j] * vectorF[j];
                MatrizResult[i] = temporal;
            }
        }
    }
}
float MultMat3x3(float JacobianaInv[][3], float *vectorF, float *MatrizResult){
    float temporal = 0;
    int i, j, k = 0;
    for (i = 0 ; i < 3 ; i++ ) //i para las filas de la matriz resultante
    {
        for (k = 0 ; k < 1 ; k++ ) // k para las columnas de la matriz resultante
        {
            temporal = 0 ;
            for(j = 0 ; j < 3 ; j++ ) //j para realizar la multiplicacion de 
            {                          //los elementos   de la matriz
                temporal += JacobianaInv[i][j] * vectorF[j];
                MatrizResult[i] = temporal;
            }
        }
    }
}

void ActualizacionVectorX(float *oldVectorX, float *vectorResta, int size, float *newVectorX){
	int i = 0;
	for(i = 0; i< size; i++){
		newVectorX[i] = oldVectorX[i] - vectorResta[i];
	}
}

void VectorDiferencias(float *oldVectorX, float *newVectorX, int size, float *vectorDiferencias){
	int i=0;
	for(i=0; i<size; i++){
		vectorDiferencias[i] = newVectorX[i] - oldVectorX[i];
	}
}

float NormaEspectral(float vectorDiferencias[], int size){
	float norma;
	int i=0;
	for(i=0; i<size-1; i++){
		if(fabs(vectorDiferencias[i]) < fabs(vectorDiferencias[i+1])){
			norma = fabs(vectorDiferencias[i]+1);
		}
		else{
			norma = fabs(vectorDiferencias[i]);
		}
	}
	return norma;
}

void ImprimirVector(float *vector, int size){
	int i = 0;
	printf("\n");
	for(i=0; i<size; i++){
			printf("%.12f ",vector[i]);
		printf("\n");
	}
	printf("\n \n");
}

void Bucle(int opcion, int size){
	int i, iteraciones, manejaIt=1;
	bool continuar = true;
	float auxiliar;
	float tolerancia, normaEspectral;
	float vectorX[size];
	float vectorXNew[size];
	float vectorDif[size];
	float vectorF[size];
	float MatrizInv[size][size];
	float MatrizJ[size][size];
	float vectorR [size];
	for(i=0; i<size; i++){
		fflush(stdin);
		printf("X_%d=", i+1); scanf("%f",&auxiliar);
		vectorX[i] = auxiliar;
	}
	printf("\nIngresa número maximo de iteraciones :"); scanf("%d",&iteraciones);
	printf("\nIngresa la tolerancia: "); scanf("%f",&tolerancia);
	system("cls");
	do{
		EvaluaFuncion(opcion, vectorX, size, vectorF);		
		if(opcion == 1 || opcion == 2){		
			Jacobiana2x2(opcion, vectorX, MatrizJ);	
			Inversa2x2(MatrizJ, MatrizInv);
			MultMat2x2(MatrizInv, vectorF, vectorR);
		}
		else{
			Jacobiana3x3(opcion, vectorX, MatrizJ);
			Inversa3x3(MatrizJ, MatrizInv);	
			MultMat3x3(MatrizInv, vectorF, vectorR);
		}
		ActualizacionVectorX(vectorX, vectorR, size, vectorXNew);
		VectorDiferencias(vectorX, vectorXNew, size, vectorDif);
		normaEspectral = NormaEspectral(vectorDif, size);
		if(opcion == 1 || opcion == 2){
			printf("Inversa");
			ImprimirM2x2(MatrizInv);
		}else{
			printf("Inversa");
			ImprimirM3x3(MatrizInv);
		}
		
		for(i=0; i<size; i++){
			vectorX[i] = vectorXNew[i];
		}	
			
		printf("Vector F");
		ImprimirVector(vectorF, size);
		printf("Vector X");
		ImprimirVector(vectorX, size);
		printf("Error \n");
		printf("%f", normaEspectral);
		printf("\n \n");
			
		manejaIt++;
		if(normaEspectral < tolerancia){
			continuar = false;
			printf("Tolerancia alcanzada \n");
		}else if(iteraciones < manejaIt){
			continuar = false;
			printf("Iteraciones Maximas alcanzadas \n");
		}else{
			continuar = true;
		//	printf("Continuar");
		}
		
	}while(continuar == true);
	
	printf("Vector Alcanzado");
	ImprimirVector(vectorX, size);
	printf("\n \n");
	int unaMas;
	printf("Desea evaluar otro sistema\n 1-si 2-no\n");
	scanf("%d", &unaMas);
	if(unaMas == 1){
		Menu1();
	}
	else{
		return 0;
	}
	
}

void Menu1(){
    int opcion1;
    system("cls");
    printf("\t\t Método de Newton-Raphson \n\n");
	printf("\t Menú\n\n");
	printf("\t Sistemas\n\n");
	printf("\t1.- f1(x,y)=x^2+xy-10=0 \t\t f2(x,y)=y+3xy^2-50=0 \n\n");
	printf("\t2.- f1(x,y)=x^2+y^2-9=0 \t\t f2(x,y)=-e^x-2y-3=0; \n\n");
	printf("\t3.- f1(x,y,z)=2x^2-4x+y^2+3z^2+6z+2=0 \t f2(x,y,z)=x^2+y^2-2y+2z^2-5=0 \tf3(x,y,z)=3x^2-12x+y^2-3z^2+8=0\n\n");
	printf("\t4.- f1(x,y,z)=x^2-4x+y^2=0 \t\t f2(x,y,z)=x^2-x-12y+1=0 \tf3(x,y,z)=3x^2-12x+y^2-3z^2+8=0\n\n");
	printf("\t5.- Salir\n\n");
	printf("Ingresa el número de la función que deseas resolver: ");
	scanf("%d",&opcion1);
    switch(opcion1){
        case 1:
        	system("cls");
			printf("\t\t Método de Newton-Raphson \n\n");
			printf("\t1.- f1(x,y)=x^2+xy-10=0 \t\t f2(x,y)=y+3xy^2-50=0 \n\n");
			printf("Ingresa los valores iniciales de X-0 \n");
			Bucle(1, 2);
            break;
        case 2:
        	system("cls");
			printf("\t\t Método de Newton-Raphson \n\n");
			printf("\t2.- f1(x,y)=x^2+y^2-9=0 \t\t f2(x,y)=-e^x-2y-3=0; \n\n");
			printf("Ingresa los valores iniciales de x-0: \n");
			Bucle(2,2);
            break;
        case 3:
        	system("cls");
			printf("\t\t Método de Newton-Raphson \n\n");
			printf("\t3.- f1(x,y,z)=2x^2-4x+y^2+3z^2+6z+2=0 \t f2(x,y,z)=x^2+y^2-2y+2z^2-5=0 \tf3(x,y,z)=3x^2-12x+y^2-3z^2+8=0\n\n");
			printf("Ingresa los valores iniciales de x-0: \n");
			Bucle(3,3);
            break;
        case 4:
        	system("cls");
			printf("\t\t Método de Newton-Raphson \n\n");
			printf("\t4.- f1(x,y,z)=x^2-4x+y^2=0 \t\t f2(x,y,z)=x^2-x-12y+1=0 \tf3(x,y,z)=3x^2-12x+y^2-3z^2+8=0\n\n");
			printf("Ingresa los valores iniciales de x_0: \n");
			Bucle(4,3);
            break;
        case 5:
        	printf("Adios");
        	return 0;
            break;
        default:
        	printf("Error, esa no es una opcion");
        	Menu1();
            break;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////DiferenciasDivididas

void DiferenciaD(float *fx,float *x, int i){
    float punto_I,**fn, error,Pn=0,Pauxiliar=0, Paux=0, auxOrd, auxFord, cicloError;
    int ciclo=0,Grado, w,y,z, ajuste, contador,contador2, ifor, xfor, opcion2;

	do{
    	system("cls");
		if(ciclo==0){
        	printf("\nTabla ordenada\n");
        	//Ordenar los datos aqui
        	for(ifor=0; ifor<i - 1; ifor++){
        		for(xfor=ifor+1; xfor<i; xfor++){
				    if(x[ifor]>x[xfor]){
            			auxOrd=x[ifor];
            			auxFord=fx[ifor];
            			x[ifor]=x[xfor];
            			fx[ifor]=fx[xfor];
            			x[xfor]=auxOrd;
            			fx[xfor]=auxFord;
        			}
				}
        	}
            for(w=0; w < i; w++){
                printf("\n   %d  %f  %f \n",w,x[w],fx[w]);
            }
            ciclo=1;
            printf("\n\n");
            system("pause");
        }
	    else{
			printf("\n  ¿Desea interpolar otro punto con la misma tabla?\n");
            printf("\n\n\t  1.Si 2.No\n");
            printf("\n");
            scanf("%d",&ciclo);
            printf("\n");
            system("pause");
		}
		if( ciclo<1 || 2<ciclo){
			system("cls");
            printf("\n\tOpción inválida");
            printf("\n\tIntente de nuevo");
            printf("\n\n");
            system("pause");
        }
		else{
			do{
			if(ciclo==1){
				system("cls");
                printf("\n\n  Punto a interpolar: ");
                scanf("%f",&punto_I);
                if(punto_I>= x[0] && punto_I<=x[i-1] ){
                	cicloError = 1;
					printf("\n\n Grado del Polinomio: ");
                    scanf("%d",&Grado);
					if(Grado>0){
                    	for(w=1; w<i ;w++){
                            if(x[w] > punto_I){
                            	ajuste = w-1;
                                break;
                            }
                        }
                        //printf("ajuste=%d",ajuste);
                        if(Grado < i- ajuste){
                        	fn = (float**)malloc(Grado*sizeof(float*));//Creo columnas
							for(w=0; w<i; w++){
								fn[w]=(float*)malloc(i*sizeof(float));//Creo renglones
							}
                            for(w=0; w < Grado; w++){
                            	for(z=0; z < i; z++){
                                	fn[w][z]=0;
                                }
                            }
                            y=1;
                            for(w=0; w < Grado; w++){
                                for(z=0; z < i-y; z++){
                                	if(w==0){
                                    	fn[w][z]= (fx[z+1] - fx[z])/(x[z+1] - x[z]);
                                    }
									else{
                                    	fn[w][z]= (fn[w-1][z+1] - fn[w-1][z])/(x[z+1+w] - x[z]);
                                    }
                                }
                                y++;
                            }
							printf("\n    i      x        fx    \n");
                            for(w=ajuste; w < i; w++){
                            	printf("\n   %d  %f  %f \n",w,x[w],fx[w]);
                            }
                            printf("\n");
                            for(w=0; w < Grado; w++){
                            	printf("    f[%i]    ",w+1);
                            }
                            printf("\n");
                            for(w=ajuste; w < i; w++){
                            	for(z=0;z < Grado; z++){
                                	printf("   %.6f  ",fn[z][w]);
                                }
                                printf("\n");
                            }
                            Pn=fx[ajuste];
                            contador2=ajuste;
                            for(z=0; z < Grado; z++){
                            	Pauxiliar=1;
                                contador=ajuste;
                                for(w=0; w <=z; w++){
                                	Pauxiliar= Pauxiliar*(punto_I - x[contador]);
                                    contador++;
                                }
                                Paux= fn[z][contador2]*Pauxiliar;
                                Pn= Paux + Pn ;
                                printf("\n P%d(x= %f) = %.6f ",z+1,punto_I,Pn);
                                if(z > 0){
                                	printf("\n\n ERTS %d = %f \n",z,Paux);
                                }
                                contador2;
                            }
                            printf("\n\n");
                            system("pause");
                            free(fn);
                            }
							else{
                                if(ajuste != 0){
                                	while(ajuste > 0 && i-ajuste <= Grado){
                                    	ajuste--;
                                    }
                                    // printf("ajuste=%d",ajuste);
                                    if(Grado < i- ajuste){
                                    	fn = (float**)malloc(Grado*sizeof(float*));//Creo columnas
                                        for(w=0; w<i; w++){
                                            fn[w]=(float*)malloc(i*sizeof(float));//Creo renglones
                                        }
                                        for(w=0; w < Grado; w++){
                                        	for(z=0;z < i; z++){
                                            	fn[w][z]=0;
                                            }
                                        }
                                        y=1;
                                        for(w=0; w < Grado; w++){
                                        	for(z=0; z < i-y; z++){
                                            	if(w==0){
                                                	fn[w][z]= (fx[z+1] - fx[z])/(x[z+1] - x[z]);
                                                }
												else{
                                                	fn[w][z]= (fn[w-1][z+1] - fn[w-1][z])/(x[z+1+w] - x[z]);
                                                }
                                            }
                                            y++;
                                        }
                                        printf("\n    i      x        fx    \n");
                                        for(w=ajuste; w < i; w++){
                                        	printf("\n   %d  %f  %f \n",w,x[w],fx[w]);
                                        }
                                        printf("\n");
                                        for(w=0; w < Grado; w++){
                                        	printf("    f[%i]    ",w+1);
                                        }
                                        printf("\n");
                                        for(w=ajuste; w < i; w++){
                                        	for(z=0;z < Grado; z++){
                                            	printf("   %.6f  ",fn[z][w]);
                                            }
                                            printf("\n");
                                        }
                                        Pn=fx[ajuste];
                                        contador2=ajuste;
                                        for(z=0; z < Grado; z++){
                                        	Pauxiliar=1;
                                            contador=ajuste;
                                            for(w=0; w <=z; w++){
                                            	Pauxiliar= Pauxiliar*(punto_I - x[contador]);
                                                contador++;
                                            }
                                            Paux= fn[z][contador2]*Pauxiliar;
                                            Pn= Paux + Pn ;
                                            printf("\n P%d(x= %f) = %.6f ",z+1,punto_I,Pn);
                                            if(z > 0){
                                            	printf("\n\n Error %d = %f \n",z,Paux);
                                            }
                                            contador2;
                                        }
                                        printf("\n\n");
                                        system("pause");
                                        free(fn);
                                    }
									else{
                                    	printf("\n Los valores que proporcionaste no alcanzan para realizar polinomio de este número de grado\n\n");
                                        system("pause");
                                    }
                            	}
								else{
                                	printf("\n Los valores que proporcionaste no alcanzan para realizar polinomio de este número de grado\n\n");
                                    system("pause");
                                }
                            }
                        }
						else{
                        	printf("\n Este número es inválido\n\n");
                            system("pause");
                        }
                }
				else{
                	system("cls");
                    printf("\n\tEste punto no está en el intervalo");
                    printf("\n\tIntente de nuevo");
                    printf("\n\n");
                    system("pause");
                    cicloError=0;
                }
            }        	
            else{
            	system("cls");
            	printf("\n\t ¿Que desea? \n\t 1.Cerrar el programa \n\t 2.Ingresar otra Tabla");
            	scanf("%d", &opcion2);
            	if(opcion2 == 2){
            		Menu_secundario();
				}
				else{
					printf("\n\n\tAdios");
				}
			}
			}while(cicloError!=1);
        }     
    }while(ciclo!=2);
 }

void portada(){
	printf("\n\tInterpolacion Polinomial");
		printf("\n\t\t Bautista Padron Oscar Raymundo\n");
	printf("\n\t\t Lechuga Ochoa Daniel\n");
	printf("\n\t\t Moreno Cruz Tatiana\n");
	printf("\n");
	system("pause");
	system("cls");
}

void Menu_secundario(){
    int Eleccion_menup=0;
    float *x,*fx;
    int w,y,tam,Correcto;
    system("cls");
    printf(" \n Ingresa el tamaño de la tabla (cantidad de x):  ");
    scanf("%d",&tam);

    if(tam>0){
        x=(float*)malloc((tam)*sizeof(float));
        fx=(float*)malloc((tam)*sizeof(float));
        for(w=0;w<tam;w++){
            x[w]=0;
            fx[w]=0;
        }
        system("cls");
        printf("\n Ingresa los valores \n");
        for(w=0; w<tam; w++){
            printf("\n  X%d = ",w);
            scanf("%f",&x[w]);
            printf("  f(X%d) = ",w);
            scanf("%f",&fx[w]);
        }

        do{
        	system("cls");
            printf("\n    i      x        fx    \n");
            for(w=0; w<tam; w++){
                printf("\n   %d  %f  %f \n",w,x[w],fx[w]);
            }
            printf("\n\n Estos datos son correctos \n\t\t 1.Si 2.No: " );
            scanf("%d", &Correcto);

            if(Correcto !=2 && Correcto !=1){
                printf("\n Opción inválida\n\n");
                system("pause");
            }else{
                if(Correcto==2){
                    system("cls");
                    printf("\n A continuación haz la corrección de la tabla\n");
					for(w=0; w<tam; w++){
                        printf("\n  X%d = ",w);
                        scanf("%f",&x[w]);
                        printf("  f(X%d) = ",w);
                        scanf("%f",&fx[w]);
                    }
            	}
            }
		}while(Correcto!=1);

			system("cls");
            printf("\n\n   1. Diferencias Divididas\n");
            printf("\n\n   2. Ajuste de Curvas Spline Cubico\n");
            printf("\n\n   3. Salir\n");
            printf("\n");
            scanf("%d",&Eleccion_menup);
            printf("\n");
            system("pause");

            if( Eleccion_menup<1 && 3<Eleccion_menup){
                system("cls");
                printf("\n\tOpción inválida");
                printf("\n\tIntente de nuevo");
                printf("\n\n");
                system("pause");
            }
			else{
                switch(Eleccion_menup){
                    case 1:
                        DiferenciaD(fx,x,tam);
                        break;
                    case 2:
                        Spline(x, fx, tam);
                        break;
                    }
                }
    }
	else{
        printf("\n Este número es inválido\n\n");
        system("pause");
    }

    free(x);
    free(fx);
}

void Menu_principal(){
    int E_M=0;
    do{
		system("cls");
        if(E_M==0){
           // printf("\n  BIENVENIDO AL PROGRAMA\n");
            E_M =1;
           // printf("\n\n");
           // system("pause");
		}
		else{
            printf("\n  ¿Desea realizar otra interpolación con otra tabla?\n");
            printf("\n\n   1.Si\n");
            printf("\n\n   2.No\n");
            scanf("%d",&E_M);
            printf("\n");
            system("pause");
        }

        if( E_M<1 || 2<E_M){
            system("cls");
            printf("\n\tOpción inválida");
            printf("\n\tIntente de nuevo");
            printf("\n\n");
            system("pause");
        }
		else{
            if(E_M==1){
                Menu_secundario();
            }
        }
    }while(E_M!=2);
}




void Spline(float *x, float *y, int tama){
	float *h, *a, *b, *c, *d, *F, *S;
	float **matriz, *vectorB, **matrizC; 
	int i, j;
	
	h = malloc((tama)*sizeof(float));
	a = malloc((tama)*sizeof(float));
	b = malloc((tama)*sizeof(float));
	c = malloc((tama)*sizeof(float));
	d = malloc((tama)*sizeof(float));
	F = malloc((tama)*sizeof(float));
	S = malloc((tama)*sizeof(float));
	
	for(i=0; i<tama-1; i++ ){
		h[i] = x[i+1] - x[i];
		
		F[i] = (y[i+1] - y[i]) / h[i];
	}
	
	matriz = malloc((tama-2)*sizeof(float*));
	matrizC = malloc((tama-2)*sizeof(float*));
	
	for(i=0; i < tama-2; i++){
		matriz[i] = malloc((tama-2)*sizeof(float));
	}
	for(i=0; i < tama-2; i++){
		matrizC[i] = malloc((tama)*sizeof(float));
	}
	
	for(i=0; i<tama-2; i++){
		for(j=0; j<tama; j++){
			if(i == j){
				matrizC[i][j] = h[i];
			}
			else{			
				if(j == i+1){
					matrizC[i][j] = 2*(h[i]+h[i+1]);
				}else{			
					if(j == i+2){
					matrizC[i][j] = h[i+1];
					}
					else{
					matrizC[i][j] = 0;
					}
				}	
			}
		}
	}
	
	for(i=0; i<tama-2; i++){
		for(j=0; j<tama-1; j++){
			if(j == tama-2){
				matriz[i][j] = 6*(F[i+1]-F[i]);
			}
			else{
				matriz[i][j] = matrizC[i][j+1]; 
			}			
			printf("\t %f", matriz[i][j]);
		}
	//	printf("\n");
	}
	system("pause");
	
}


void SolucionGauss(float **m, float *S, int tam){


}

void poli(float **x1, float *fx, int tamanio, float *x){
	printf("llamo a poli");
	system("pause");
    float **Evaluacion, **Diferencia,*Suma1;
    int i,j,w;
    printf("\n Polinomio de tercer grado:\n\n P3(x)= %.6f + (%.6f)x + (%.6f)x^2 + (%.6f)x^3 \n",x1[0][2], x1[1][2], x1[2][2], x1[3][2]);
    printf("\n");
    printf("Polinomio de segundo grado:\n\n P2(x)= %.6f + (%.6f)x + (%.6f)x^2\n",x1[0][1], x1[1][1], x1[2][1]);
    printf("\n");
    printf("Polinomio lineal:\n\n P(x)= %.6f + (%.6f)x \n",x1[0][0], x1[1][0]);

    Evaluacion=malloc((tamanio)*sizeof(float*));
    Diferencia=malloc((tamanio)*sizeof(float*));
    Suma1=malloc((3)*sizeof(float));
    for(i=0;i<tamanio;i++){
        Evaluacion[i]=malloc((3)*sizeof(float));
        Diferencia[i]=malloc((3)*sizeof(float));
    }

    for(i=0;i<tamanio;i++){
        for(j=0;j<3;j++){
            if(j==0){
            	Evaluacion[i][j]=x1[0][0]+(x1[1][0]*x[i]);
            }
            if(j==1){
            	Evaluacion[i][j]=x1[0][1]+(x1[1][1]*x[i])+( x1[2][1]*pow(x[i],2));
            }

            if(j==2){
            	Evaluacion[i][j]=x1[0][2]+(x1[1][2]*x[i])+( x1[2][2]*pow(x[i],2))+(x1[3][2]*pow(x[i],3));
            }
        }
    }
    printf("\n");
	for(w=0; w < 3; w++){
          printf("     P(%d)      |",w+1);
      }
    printf("\n\n");
	for(i=0;i<tamanio;i++){
    	for(j=0;j<3;j++){
        	printf("   %.6f    |",Evaluacion[i][j]);//Es la evaluación de los puntos X en el polinomio
        }
        printf("\n");
    }
     for(i=0;i<tamanio;i++){
        for(j=0;j<3;j++){
        	if(j==0){
            	Diferencia[i][j]=pow(fabs(fx[i]-Evaluacion[i][j]),2);
            }
            if(j==1){
            	Diferencia[i][j]=pow(fabs(fx[i]-Evaluacion[i][j]),2);
            }
            if(j==2){
            	Diferencia[i][j]=pow(fabs(fx[i]-Evaluacion[i][j]),2);
            }
        }
    }
    for(i=0;i<3;i++){
        Suma1[i]=0;
        for(j=0;j<tamanio;j++){
            Suma1[i]=Suma1[i]+Diferencia[j][i];
        }
    }
	printf("\n");
    for(w=0; w < 3; w++){
          printf("     ( y - P(%d))^2  |",w+1);
    }
	printf("\n\n");

    for(i=0;i<tamanio;i++){
        for(j=0;j<3;j++){
        	printf("        %.6f    |",Diferencia[i][j]);//Es la diferencia de los puntos evaluados en el polinomio con Y
        }
        printf("\n");
    }

	printf("\n\n Suma");
    for(j=0;j<3;j++){
        printf("   %.6f     |",Suma1[j]);
    }
	if(Suma1[0] < Suma1[1]){
        if(Suma1[0] < Suma1[2]){
            printf("\n\n El polinomio que nos conviene es el Lineal");
        }
		else{
            printf("\n\n El polinomio que nos conviene es el Cuadrático");
        }
    }
	else{
        if(Suma1[1] < Suma1[2]){
            printf("\n\n El polinomio que nos conviene es el Cuadrático");
        }
		else{
            if(Suma1[2]<Suma1[0]){
                printf("\n\n El polinomio que nos conviene es el Cúbico");
            }
			else{
                printf("\n\n El polinomio que nos conviene es el Lineal");
            }
        }
    }
	printf("\n\n");
	free(Evaluacion);
	free(Diferencia);
	free(x1);
	free(Suma1);
}

void gauss(float **MatrizPoli, float *VectorSol, float *fx, int tamanio, float *x){
	printf("Entramos a gauss");
    float **matriz2,*resultado2,**x1,**matriz,*resultado;
    int l=4,m;
    int i,j,k;
	matriz2=malloc((4)*sizeof(float*)); //Creo memoria para cada matriz que va a ser utilizada
	matriz=malloc((4)*sizeof(float*));
	resultado2=malloc((4)*sizeof(float));
	resultado=malloc((4)*sizeof(float));
	x1 = malloc((4)*sizeof(float*));

    for(i=0;i<4;i++){
        matriz2[i]=malloc((4)*sizeof(float));
    }
    for(i=0;i<4;i++){
        matriz[i]=malloc((4)*sizeof(float));
    }

	for(i=0;i<4;i++){
        x1[i]=malloc((4)*sizeof(float));
    }
     for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			x1[i][j]=0;
		}
    }
	//valor matriz2
    for(i=0;i<4;i++){
		for(j=0;j<4;j++){
			matriz2[i][j]=MatrizPoli[i][j];
			matriz[i][j]=MatrizPoli[i][j];
		}
    }

 //VALORES VECTOR SOLUCIÓN
 	for(i=0;i<4;i++){
        resultado2[i]=VectorSol[i];
        resultado[i]=VectorSol[i];
	}

//MÉTODO GAUSS
	for(i=0;i<3;i++)
		/*{
		    printf("[ ");
			for(j=0;j<3;j++)
			{
				printf(" %.2f ",x[i][j]);
			}
			printf(" ]");
        }*/
        //Hasta esta parte todo está bien
	for(m=0;m<3;m++){
        for(k=0;k<l-1;k++){
            for(i=k+1;i<l;i++){
                for(j=k;j<l;j++){
                   /* printf(" %f - ",matriz[i][j]);
                    printf("( %f ",matriz[k][j]);
                    printf(" *%f )",matriz[i][k]);
                    printf(" /%f ",matriz[k][k]);*/
                    matriz2[i][j]=MatrizPoli[i][j]-(MatrizPoli[k][j]*MatrizPoli[i][k]/MatrizPoli[k][k]);
                	if(j==l-1){
						resultado2[i]=VectorSol[i]-(VectorSol[k]*MatrizPoli[i][k]/MatrizPoli[k][k]);
                    }
                    /*printf("++%f ",matriz2[i][j]);
                    system("pause");*/
                }
            }
            for(i=0;i<4;i++){
                for(j=0;j<4;j++){
				MatrizPoli[i][j]=matriz2[i][j];
                }
			}
        	for(i=0;i<4;i++){
            	VectorSol[i]=resultado2[i];
			}
    	}
        /*printf("\n  Valores de la matriz\n\n");
		for(i=0;i<fila;i++)
		{
		    printf(" |");
			for(j=0;j<fila;j++)
			{
				printf(" %.2f ",matriz2[i][j]);
			}
			printf("|= ");
			printf("%.2f",resultado2[i]);
			printf("\n");
			}*/

		for(i=l-1;i>=0;i--){
           // printf("++%f ",resultado2[i]);
            x1[i][l-2]= resultado2[i];

 // system("pause");
        	for(j=l-1; j>i; j--){
            	x1[i][l-2]=x1[i][l-2]-matriz2[i][j]*x1[j][l-2];
        	}
//printf("%f",x[i]);
        	x1[i][l-2]=x1[i][l-2]/matriz2[i][i];
		}
        l--;
        for(i=0;i<4;i++){
            for(j=0;j<4;j++){
				MatrizPoli[i][j]=matriz[i][j];
            }
		}
        for(i=0;i<4;i++){
            VectorSol[i]=resultado[i];
		}
	}
	
	poli(x1,fx,tamanio,x);
    free(resultado2);
    free(matriz2);
    free(matriz);
    free(resultado);
}

void minimos(float *fx, float *x, int i){
    float **matriz1,**matriz2,**matriz3,**MatrizPot,**MatrizMult,**MatrizPoli,*VectorSol;
    float *Suma;
    int j,k,cont,tamanio,w;
    MatrizPot=malloc((i)*sizeof(float*));
    MatrizMult=malloc((i)*sizeof(float*));
    MatrizPoli=malloc((4)*sizeof(float*));
    Suma=malloc((10)*sizeof(float));
    VectorSol=malloc((3)*sizeof(float));
    tamanio=i;
    for(j=0;j<i;j++){
        MatrizPot[j]=(float*)malloc((5)*sizeof(float));
    }
    for(j=0;j<i;j++){
        MatrizMult[j]=(float*)malloc((2)*sizeof(float));
    }
    for(j=0;j<4;j++){
        MatrizPoli[j]=(float*)malloc((4)*sizeof(float));
    }

    for(j=0;j<i;j++){    //Matriz potencias de X
        for(k=0;k<5;k++){
            MatrizPot[j][k]=pow(x[j],k+2);
        }
    }
    for(j=0;j<i;j++){    //Matriz potencias de X*Y
        for(k=0;k<3;k++){
            MatrizMult[j][k]=pow(x[j],k+1)*fx[j];
        }
    }

    system("cls");

    printf("\n   | i |     x    |    fx   | \n");
    for(w=0;w < i;w++){
         printf("\n   |%d  |%f  |%f |\n",w+1,x[w],fx[w]);
     }
     printf("\n");
//Aquí poner las sumas de X y fx.Imprimir las tablas x y fx y sus respectivas sumas

    Suma[8]=0;
    for(j=0;j<i;j++){
        Suma[8]=Suma[8]+x[j];
        //Suma[8] tiene el valor de la suma del vector de las X
    }

    Suma[9]=0;
    for(j=0;j<i;j++){
        Suma[9]=Suma[9]+fx[j];
        //Suma[9] tiene el valor de la suma del vector de las fx
    }

    for(j=0;j<5;j++){
        Suma[j]=0;
        for(k=0;k<i;k++){
            Suma[j]=Suma[j]+MatrizPot[k][j];
        }
    }

    cont=0;
    for(j=5;j<8;j++){
        Suma[j]=0;
        for(k=0;k<i;k++){
            Suma[j]=Suma[j]+MatrizMult[k][cont];
        }
    	cont++;
    }

    VectorSol[0]=Suma[9];
    cont=1;
    for(j=5;j<8;j++){
    	VectorSol[cont]=Suma[j];
        cont++;
    }

    printf("\n Suma ");
    printf("  |%.6f  |%.6f |\n\n",Suma[8],Suma[9]);
	printf("\n");

    for(w=2; w < 7; w++){
          printf("   Xi^%d     |",w);
    }
    printf("\n");

    for(j=0;j<i;j++){
    	for(k=0;k<5;k++){
            printf(" %.6f  |",MatrizPot[j][k]);//Esta matriz muestra al vector X elevado a cierta potencia
        }
    	printf("\n");
    }
    printf("\n");
    printf("\n Suma");
    for(j=0;j<5;j++){
		printf("  %.6f   |",Suma[j]);//Imprime el valor de las Sumas pero aqui Suma[8] tiene el valor de las sumas de X
		                               //Y el valor de las sumas de Y está en Suma[9]
        }

    printf("\n\n");

    for(w=0; w < 3; w++){
    	printf(" Xi^%d Y   |",w+1);
    }

    printf("\n");
    printf("\n");
    for(j=0;j<i;j++){
    	for(k=0;k<3;k++){
        	printf(" %.6f |",MatrizMult[j][k]);//Esta matriz muestra al vector X multiplicando al vector Y
        }
		printf("\n");
    }

    printf("\n Suma");
    for(j=5;j<8;j++){
        printf("  %.6f   |",Suma[j]);//Imprime el valor de las Sumas pero aqui Suma[8] tiene el valor de las sumas de X
                                       //Y el valor de las sumas de Y está en Suma[9]
    }

    MatrizPoli[0][0]=i;
    MatrizPoli[0][1]=Suma[8];
    MatrizPoli[0][2]=Suma[0];
    MatrizPoli[0][3]=Suma[1];
    MatrizPoli[1][0]=Suma[8];
    MatrizPoli[1][1]=Suma[0];
    MatrizPoli[1][2]=Suma[1];
    MatrizPoli[1][3]=Suma[2];
    MatrizPoli[2][0]=Suma[0];
    MatrizPoli[2][1]=Suma[1];
    MatrizPoli[2][2]=Suma[2];
    MatrizPoli[2][3]=Suma[3];
    MatrizPoli[3][0]=Suma[1];
    MatrizPoli[3][1]=Suma[2];
    MatrizPoli[3][2]=Suma[3];
    MatrizPoli[3][3]=Suma[4];
    printf("\n\n");
    system("pause");

    gauss(MatrizPoli,VectorSol, fx,tamanio,x);

    free(matriz1);
    free(matriz2);
    free(matriz3);
     free(MatrizPot);
    free(MatrizMult);
    free(MatrizPoli);
    free(VectorSol);
    free(Suma);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Integración

void Menu_principal2()
{
    int Eleccion_menup=0;
  do
   {

    
    Eleccion_menup=1;
    printf("\n");
    getch();

    if( Eleccion_menup<1 || 4<Eleccion_menup)
    {
        system("cls");
        printf("\n\tOpción inválida");
        printf("\n\tIntente de nuevo");
        printf("\n");
        
    }else
          {
             switch(Eleccion_menup)
             {

             
             case 1:
                Menu_principal_DI();
                break;
             }
          }
  }while(Eleccion_menup!=4);
}
void Integracion(float *fx,float *x,int i)
{
     int Com=0;
     float h,SimpsonU=0,SimpsonCom=0,SimpsonComD = 0;
     int w;

     h = x[1] - x[0];

    system("cls");
    if((i-1)%2 == 0)
    {
           //printf("%f,%f\n",fx[0],fx[i-1]);
           SimpsonU=fx[0]+fx[i-1];

           for(w=1; w<i-1 ;w++)
           {
               if((w%2)== 0)
               {
                   //printf("\n+%f",SimpsonU);system("pause");
                   SimpsonU= SimpsonU + (2*fx[w]);
               }else
                 {
                    // printf("\n-%f",SimpsonU);system("pause");
                   SimpsonU=SimpsonU + (4*fx[w]);
                 }
           }

           //printf("1/3= %f\n", SimpsonU);system("pause");
           SimpsonU = (SimpsonU*h)/3;

           printf("\n    i |     x    |    fx    \n");
                             for(w=0;w<i;w++)
                             {
                                 printf("\n   %d  |%f  |%f \n",w,x[w],fx[w]);
                             }

           printf("\n Simpson 1/3 = %.6f \n",SimpsonU);
           printf("\n\n");
           system("pause");
    }else
       {
            Com= (i-1) - 3;
           // printf("\n com= %d",Com);
            SimpsonCom= fx[0] + fx[Com];
            for(w =1 ;w < Com; w++)
            {
                if((w%2)== 0)
                   {
                       //printf("\n+%f",SimpsonU);system("pause");
                       SimpsonCom= SimpsonCom + (2*fx[w]);
                   }else
                     {
                        // printf("\n-%f",SimpsonU);system("pause");
                       SimpsonCom=SimpsonCom + (4*fx[w]);
                     }
            }

            SimpsonComD= fx[Com] + fx[i-1];
            for(w=Com+1; w < i-1;w++)
            {
                if((w%3)==0)
                   {
                       //printf("\n++%f\n", SimpsonD);system("pause");
                       SimpsonComD = SimpsonComD + (2*fx[w]);
                   }else
                      {
                          //printf("--%f\n", SimpsonD);system("pause");
                          SimpsonComD = SimpsonComD + (3*fx[w]);
                      }
            }


              SimpsonCom = (SimpsonCom*h)/3;

              SimpsonComD = (SimpsonComD*h*3)/8;
            // printf("\n 1/3= %f\n", SimpsonCom);system("pause");
             //printf("3/8= %f\n", SimpsonComD);system("pause");

              SimpsonCom = SimpsonCom + SimpsonComD;

               printf("\n    i |     x    |    fx    \n");
                                 for(w=0;w<i;w++)
                                 {
                                     printf("\n   %d  |%f  |%f \n",w,x[w],fx[w]);
                                 }

               printf("\n Simpson 1/3 con 3/8 = %.6f \n",SimpsonCom);
                   printf("\n\n");
                system("pause");
       }
}

void Derivacion(float *fx,float *x,int i)
{
    float intervaloA, intervaloB;
    int k,s,w,cont=0;
    float h;
    float *PDerivada,*SDerivada;

    PDerivada = (float*)malloc((i)*sizeof(float));
    SDerivada = (float*)malloc((i)*sizeof(float));

    for(k=0;k<i;k++)
    {
        PDerivada[k]=0;
        SDerivada[k]=0;
    }

    h = x[1] - x[0];

    system("cls");
    printf("\n Ingrese el intervalo de valores para calcular sus derivadas :\n");
    printf("\n comienza en x =  ");
    scanf("%f",&intervaloA);
    printf("\n Termina en x = ");
    scanf("%f", &intervaloB);

    for(w=0;w<i;w++)
    {
        if(intervaloB == x[w])
        {
            break;
        }
         cont++;
    }

    for(k=0; k<i; k++)
    {
        //printf("(%f,%f)", intervaloA,x[k]);
        if(intervaloA == x[k])

        {
            //printf("K =%d,** %d",k,cont);
            //system("pause");
            for(s=k; s <= cont; s++)
            {
                if(s!=0 && s!=(i-1))
                {
                    //printf("Derivada h=%f,fx-1=%f,fx+1=%f",h,fx[s+1],fx[s-1]);
                    PDerivada[s]=(fx[s+1]-fx[s-1])/(2*h);
                    SDerivada[s]=(fx[s+1]-2*fx[s]+fx[s-1])/(pow(h,2));
                    //system("pause");
                }
            }
            break;
        }
    }

    printf("\n    i |     x    |    fx   |    f'(x)    |    f''(x)    \n");
                     for(w=0;w<i;w++)
                     {
                        if(w!=0 && w!=(i-1))
                        {
                            printf("\n   %d  |%f  |%f   |%f  | %f  \n",w,x[w],fx[w],PDerivada[w],SDerivada[w]);
                        }else
                           {
                               printf("\n   %d  |%f  |%f   | ---   | ---   \n",w,x[w],fx[w]);
                           }

                     }
    printf("\n\n");
    system("pause");
    free(PDerivada);
    free(SDerivada);
}

void Menu_secundarioDI()
{
     int Eleccion_menup=0;
     float *x,*fx, resta1=0,resta2=0;
     int w,y,i,Correcto,h =0;

     system("cls");
     printf(" \n Ingresa la cantidad de valores de la tabla : ");
     scanf("%d",&i);

     if(i>0)
     {

             x=(float*)malloc((i)*sizeof(float));
             fx=(float*)malloc((i)*sizeof(float));

                 h =0;
                 resta1=0;
                 resta2=0;

                    for(w=0;w<i;w++)
                     {
                         x[w]=0;
                         fx[w]=0;
                     }

                     system("cls");
                     printf("\n Ingresa los datos de la tabla\n");

                     for(w=0;w<i;w++)
                     {
                         printf("\n  X%d = ",w);
                         scanf("%f",&x[w]);
                         printf("  f(X%d) = ",w);
                         scanf("%f",&fx[w]);
                     }

                     for(w=0; w<i-2; w++)
                     {
                         resta1=x[w+1]-x[w];
                         resta2=x[w+2]- x[w+1];

                             printf("\n\n X%d-X%d = %f \n X%d-X%d = %f",w+1,w,(resta1),w+2,w+1, (resta2));
                     }

                 resta1=0;
                 resta2=0;
                     for(w=0; w<i-2; w++)
                     {
                         resta1=x[w+1]-x[w];
                         resta2=x[w+2]- x[w+1];

                            if(resta1 == resta2)
                             {
                                 h=1;
                                // printf("++%d",h);
                             }else
                             {
                                 h=0;
                                 break;
                             }

                            // system("pause");
                     }
                    


             do
             {

                         system("cls");

                                 printf("\n    i |     x    |    fx    \n");
                             for(w=0;w<i;w++)
                             {
                                 printf("\n   %d  |%f  |%f \n",w,x[w],fx[w]);
                             }
                             printf("\n\n Son correctos los datos?  1. SI 2. NO: " );
                             scanf("%d", &Correcto);

                             if(Correcto !=2 && Correcto !=1)
                             {
                                printf("\n Opción inválida\n\n");
                                system("pause");
                             }else{
                                     if(Correcto==2)
                                     {
                                             system("cls");
                                             printf("\n  corrección de la tabla\n");
                                                for(w=0;w<i;w++)
                                                 {
                                                     x[w]=0;
                                                     fx[w]=0;
                                                 }


                                             for(w=0;w<i;w++)
                                             {
                                                 printf("\n  X%d = ",w);
                                                 scanf("%f",&x[w]);
                                                 printf("  f(X%d) = ",w);
                                                 scanf("%f",&fx[w]);
                                             }
                                             resta1=0;
                                             resta2=0;

                                             for(w=0; w<i-2; w++)
                                             {
                                                 resta1=x[w+1]-x[w];
                                                 resta2=x[w+2]- x[w+1];

                                                     printf("\n\n X%d-X%d = %f \n X%d-X%d = %f",w+1,w,(resta1),w+2,w+1, (resta2));
                                             }

                                             resta1=0;
                                             resta2=0;

                                             for(w=0; w<i-2; w++)
                                             {
                                                 resta1=x[w+1]-x[w];
                                                 resta2= x[w+2]-x[w+1];

                                                 printf("\n\n X%d-X%d = %f \n X%d-X%d = %f",w+1,w,(resta1),w+2,w+1, (resta2));


                                                        if(resta1 == resta2)
                                                         {
                                                             h=1;
                                                             //printf("++%d",h);
                                                         }else
                                                         {
                                                             h=0;
                                                         }

                                               //system("pause");
                                             }
                                             if(h==0)
                                             {
                                                // printf("\n\n Estos datos no están igualmente espaciados\n" );
                                                 //printf("\n Intente de nuevo, Estos datos no son correctos\n" );
                                                 //printf("\n Si tus datos están igualmente espaciados(Todas las restas iguales)\n");
                                                 printf("\n IGNORA ESTE MENSAJE.Tus datos son correctos\n\n");
                                                 system("pause");
                                             }

                                     }
                                 }

             }while(Correcto!=1);

            do
               {

                system("cls");

                printf("\n \n  DERIVACION E INTEGRACION\n\n");
                printf("\n\n   1. Derivacion por Diferencias centradas\n");
                printf("\n\n   2. Integracion \n");
                printf("\n\n   3. salir\n");
                printf("\n");

                printf("\n Ingrese una opcion: ");
                scanf("%d",&Eleccion_menup);
                printf("\n");
                system("pause");

                if( Eleccion_menup<1 && 3<Eleccion_menup)
                {
                    system("cls");
                    printf("\n\tOpcion invalida");
                    printf("\n\tIntente de nuevo");
                    printf("\n\n");
                    system("pause");
                }else
                      {
                         switch(Eleccion_menup)
                         {

                         case 1:
                            Derivacion(fx,x,i);
                            break;
                         case 2:
                            Integracion(fx, x, i);
                            break;
                         }
                      }

              }while(Eleccion_menup!=3);


     }else
         {
            printf("\n Este número es inválido\n\n");
            system("pause");
         }

             free(x);
             free(fx);
}

void Menu_principal_DI()
{
    int E_M=0;

          do
               {

                system("cls");
                if(E_M==0)
                   {
                        printf("\n  BIENVENIDO AL PROGRAMA\n");
                        E_M =1;
                        printf("\n\n");
                        system("pause");

                   }else{
                            printf("\n  ¿ Desea ingresar otra tabla?\n");
                            printf("\n\n   1. SI\n");
                            printf("\n\n   2. NO \n");
                            printf("\n");

                            printf("\n Proporcione el número de su elección: ");
                            scanf("%d",&E_M);
                            printf("\n");
                            system("pause");
                     }

                if( E_M<1 || 2<E_M)
                {
                    system("cls");
                    printf("\n\tOpción inválida");
                    printf("\n\tIntente de nuevo");
                    printf("\n\n");
                    system("pause");
                }else
                      {
                         if(E_M==1)
                         {
                             Menu_secundarioDI();
                         }
                      }
              }while(E_M!=2);
}

void main()
{
    int n, opcion=0;
    
    printf("\n\n Proyecto Final \n");
	printf("\n\n Elaborado por los Liscenciados:\n");
    printf("\n\n Lic. Rasgado Celaya Julio Martin");
    printf("\n\n Lic. Hernandez Gonzalez Stephanie");
    printf("\n\n Lic. Perez Ibarra Alfredo de Jesus ");
    printf("\n\n\n");
	system("pause");

    while ( opcion != 4 )
    {       
        system("cls");
            
        printf( "\n   1. Introduccion a la solucion de ecuaciones no lineales.");
        printf( "\n      Metodo de Newton. ");
        printf( "\n\n   2. Interpolacion y aproximacion polinomial. ");
        printf( "\n      Metodo de Diferencias Divididas & Spline Cubico.");
        printf( "\n\n   3. Diferenciacion e integracion numerica.");
        printf( "\n      Metodo de Diferenciacion Centrada e Integracion.");
        printf( "\n\n   4. Salir." );
        printf( "\n\n   Introduzca opci%cn (1-4): ", 162 );

        scanf( "%d", &opcion );

        /* Inicio del anidamiento */

        switch ( opcion )
        {
            case 1: Menu1();
                    break;
            
            case 2: Menu_principal();
			        break;  
			
			case 3: Menu_principal2();
					break;
			default:
					system("cls"); 
					printf("\n\n----- Opcion no valida -----\n\n");
					system("pause");
					break;
         }

         /* Fin del anidamiento */

    }

    return 0;
}
