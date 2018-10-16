#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "./lib/clase_ag.cpp"	//incluyo la libreria que está en carpeta contigua
#include <math.h>


#define Nmax			32	//Nmax es el numero de antenas maximas (tiene que ser potencia de 8)
#define N_BYTES         Nmax*2 + Nmax/8 + Nmax/4  // longitud máxima del string (cromosoma). se utiliza a efectos de reserva de espacio
#define N_INDIVIDUOS    1000	//cantidad de individuos que forman a la población
#define NGENERACIONES	1500

#define Frecuencia_Tx	900*pow(10,6) 	//Frecuencia de transmision
#define C_luz			3*pow(10,8) 	//Velocidad de la luz
#define Srx				-80				//Sensibilidad del receptor en dbm
#define Pi 				3.14
#define PRESUPUESTO		10000			//cantidad de dinero disponible para gastar

typedef unsigned char uchar;
typedef unsigned long ulong;

// variables
float *X;
float *Y;
//ganancia mas potencia 
int EIRP[4]={12,15,18,21};

int Costos[4]={100,200,250,500};
float Area_total = 0;
int Nvert = 0;
float Xmin, Xmax, Ymin, Ymax;


CLASECromosoma a, b, c, d; //varios cromosomas auxiliares
CLASEAGenetico ag;	//población de cromosomas a ser cruzadas en el AG

FILE *archivo;		//archivo de salida donde se guardan variables de monitoreo
FILE *apArchivo;
FILE *Coordenadas;


// DECLARAR funciones (prototipo de las funciones que se definen al final del código)
float evaluacion(uchar *);	//función que recibe el string que forma al cromosoma, lo decodifica, y devuelve su fitness
void mostrar(uchar *);	//función que recibe al cromosoma y lo muestra en pantalla (visualización)
int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy);//funcion que verifica si un punto se encuentra dentro de un poligono
float distancia_pol(int nvert, float *X, float *Y, float Px, float Py);//funcion que mide la menor distancia desde un punto al poligono
float calcularArea(float * X, float *Y, int nvert);//funcion que calcula el area del poligono
float Overlap(float x0, float y0, float radio0, float x1, float y1, float radio1);// funcion que calcula el solapamiento de dos circulos
float calcular_radio(int EIRP);// funcion que calcula el radio
float calcular_minimo(float *X);//funcion que calcula el valor minimo de un vector
float calcular_maximo(float *X);//funcion que calcula el valor maximo de un vector

int main(int argc, char *argv[])
{

	if(argc!=3)
	{
		printf("Estan mal los parametros (./proyecto Numero_de_vertices Nombre_mapa)\n");
		return 0;
	}

	//el primer parametro es el numero de vertices del polinomio
	Nvert=atoi(argv[1]);
	
	X=(float *)calloc(Nvert,sizeof(float));
	Y=(float *)calloc(Nvert,sizeof(float));

	char Nombre_archivo[50];
	sprintf(Nombre_archivo,"Salida_%s",argv[2]);

	char Archivo_Coordenadas[50];
	sprintf(Archivo_Coordenadas,"Coordenadas_%s",argv[2]);

	char Nombre_apArchivo[50];

	strcpy(Nombre_apArchivo, argv[2]);

	//apuntador auxiliar a un cromosoma
    CLASECromosoma *apuntador = NULL;

    ag.InicializarSemillaRand(true);	//incializo la semilla aleatoria
    ag.PoblacionAleatoria(N_INDIVIDUOS, N_BYTES);	//configuración de la población inicial
  
    ag.Inicializar_poblacion(N_INDIVIDUOS,Nmax); //inicializa los individuos

    ag.FuncionEvaluacion(evaluacion);	//indico cual es la función de decodificación/fitness
    ag.dTasaMutacion = 0.01f;			//tasa de mutación
    ag.usarElitismo = true;
    float contadorF = 0;
    int contadorC = 0;

    //archivo de excel
    archivo = fopen(Nombre_archivo, "w+");
    
    apArchivo = fopen (Nombre_apArchivo, "r");
    importTabFile (apArchivo, X, Y, Nvert);
   	fclose(apArchivo);

   	//dejamos un margen de 10% para el acote del espacio de busqueda
   	Xmin = calcular_minimo(X)*1.1; 
   	Xmax = calcular_maximo(X)*1.1; 
   	Ymin = calcular_minimo(Y)*1.1; 
   	Ymax = calcular_maximo(Y)*1.1; 
   	
   	Coordenadas= fopen(Archivo_Coordenadas,"w");

   	Area_total=calcularArea(X,Y,Nvert);

   	if(Area_total>(pow(calcular_radio(EIRP[3]),2)* Pi) * Nmax)
   	{
   		Area_total = pow(calcular_radio(EIRP[3]),2)* Pi * Nmax;
   	}

    fprintf(archivo, "Mejor:\tPoblacion:\n"); // \n\r para windows
   	fprintf(Coordenadas, "Px:\tPy:\tdiametro Mayor:\tDiametro menor\tTipo de antena:\tN_antenas:\tcosto:\n");

    //inicio ciclo de entrenamiento
    do{
		
        ag.Generacion();	//itero una generación: selección -> cruce -> mutación
        apuntador = ag.SeleccionarMejor();	//selecciono al mejor
	
        if (contadorF == apuntador->Fitness())
        {
        	contadorC ++;
        }
        else
        {
        	ag.dTasaMutacion = 0.01f;
        	contadorC = 0;
        }

        if (contadorC % 500 == 0 && contadorC != 0)
        {
        	ag.dTasaMutacion = ag.dTasaMutacion * 1.5;
        	printf("Aumento mutacion a: %.3lf\niteracion: %d\nfitness==%lf\n",ag.dTasaMutacion,(int)ag.Edad(),apuntador->Fitness());
        }

        if (ag.dTasaMutacion >= 0.2)
        {
        	ag.dTasaMutacion = 0.01f;
        }
        contadorF = apuntador->Fitness();

        fprintf (archivo, "%.5lf\t%.5lf\n", apuntador->Fitness(), ag.Fitness());
		//muestro el contenido del cromosoma del mejor
    
    }
    while((ag.Edad() < NGENERACIONES) );    // condicion de parada

    mostrar(apuntador->Cromosoma);
    
    //cierro al archivo correctamente
	fclose(archivo);
	fclose(Coordenadas);
    return 0;
}

float calcular_minimo(float *X){
	float aux=9999;
	for (int i = 0; i < Nvert; i++)
	{
		if (X[i]<aux)
		{
			aux = X[i];
		}
	}

	return aux;
}

float calcular_maximo(float *X){
	float aux=0;
	for (int i = 0; i < Nvert; i++)
	{
		if (X[i]>aux)
		{
			aux = X[i];
		}
	}

	return aux;
}


// DEFINICION de funciones

/*
 ** ----------------------------------------------------------------------------
 **     Nombre      : 
 **     Función     : 
 **     Parámetros  : 
 **     Retorna     : 
 ** ----------------------------------------------------------------------------
*/
float calcular_radio (int EIRP) 
{
	return (pow(10,(float)(-Srx+EIRP)/20.0)*(float)(C_luz/(4*Pi*Frecuencia_Tx)))/1000.0;
}

/*
 ** ----------------------------------------------------------------------------
 **     Nombre      : 
 **     Función     : 
 **     Parámetros  : 
 **     Retorna     : 
 ** ----------------------------------------------------------------------------
*/
int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) && 
    	(testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

/*
 ** ----------------------------------------------------------------------------
 **     Nombre      : 
 **     Función     : 
 **     Parámetros  : 
 **     Retorna     : 
 ** ----------------------------------------------------------------------------
*/
float distancia_pol(int nvert, float *X, float *Y, float Px, float Py)
{
	float t=0,Proy_x,Proy_y, distancia=0, min=999999;
	for (int i = 0; i < nvert ; i++)
	{
		if(i==3)
		{	
			t=((Px-X[i])*(X[0]-X[i])+(Py-Y[i])*(Y[0]-Y[i]))/(pow(X[i]-X[0],2)+pow(Y[i]-Y[0],2));

			if(t > 1)
				t=1;
			if(t < 0)
				t=0;
			
			Proy_x= X[i]+t*(X[0]-X[i]);
			Proy_y= Y[i]+t*(Y[0]-Y[i]);
	
			distancia=sqrt(pow(Px - Proy_x,2)+pow(Py - Proy_y,2));
			if(distancia<min)
				min = distancia;
		}
		else
		{
			t=((Px-X[i])*(X[i+1]-X[i])+(Py-Y[i])*(Y[i+1]-Y[i]))/(pow(X[i]-X[i+1],2)+pow(Y[i]-Y[i+1],2));
			
			if(t > 1)
				t=1;
			if(t < 0)
				t=0;
			
			Proy_x= X[i]+t*(X[i+1]-X[i]);
			Proy_y= Y[i]+t*(Y[i+1]-Y[i]);
	
			distancia=sqrt(pow(Px - Proy_x,2)+pow(Py - Proy_y,2));
			if(distancia<min)
				min = distancia;
		}
	}
	return min;
}

/*
 ** ----------------------------------------------------------------------------
 **     Nombre      : 
 **     Función     : 
 **     Parámetros  : 
 **     Retorna     : 
 ** ----------------------------------------------------------------------------
*/
float calcularArea(float * X, float *Y, int nvert)
{
	float Area=0;
	for (int i = 0; i < nvert-1; i++)
	{
		Area=Area+(X[i]*Y[i+1]-Y[i]*X[i+1]);
	}
	Area=Area+X[nvert-1]*Y[0]-Y[nvert-1]*X[0];
	Area=fabsf(Area/2);
	return Area;
}

/*
 ** ----------------------------------------------------------------------------
 **     Nombre      : 
 **     Función     : 
 **     Parámetros  : 
 **     Retorna     : 
 ** ----------------------------------------------------------------------------
*/
float Overlap(float x0, float y0, float radio0, float x1, float y1, float radio1)
{
	float Area_solapada = 0;
	
	float c = sqrt(pow(x1-x0,2)+pow(y1-y0,2));
	float theta0 = 2*acos((c*c+radio0*radio0-radio1*radio1)/(2*c*radio0));
	float theta1 = 2*acos((c*c+radio1*radio1-radio0*radio0)/(2*c*radio1));

	Area_solapada=	0.5*pow(radio0,2)*theta0-sin(theta0/2)*cos(theta0/2)*pow(radio0,2) +
				  	0.5*pow(radio1,2)*theta1-sin(theta1/2)*cos(theta1/2)*pow(radio1,2);

	return Area_solapada;
}

/*
 ** ----------------------------------------------------------------------------
 **     Nombre      : 
 **     Función     : 
 **     Parámetros  : 
 **     Retorna     : 
 ** ----------------------------------------------------------------------------
*/

float Circulo_Fuera_P(float radio, float px, float py)
{
	float penalizacion=0;
	if(!pnpoly(Nvert, X, Y, px+radio, py+radio))
		penalizacion+=0.25;
	if(!pnpoly(Nvert, X, Y, px+radio, py-radio))
		penalizacion+=0.25;
	if(!pnpoly(Nvert, X, Y, px-radio, py+radio))
		penalizacion+=0.25;
	if(!pnpoly(Nvert, X, Y, px-radio, py-radio))
		penalizacion+=0.25;

	return penalizacion;
}

/*
 ** ----------------------------------------------------------------------------
 **     Nombre      : 
 **     Función     : 
 **     Parámetros  : 
 **     Retorna     : 
 ** ----------------------------------------------------------------------------
*/
float evaluacion(uchar* cromosoma)
{
    float x_g[Nmax];
	float y_g[Nmax];
	float radio[Nmax]={};
	int auxiliarP;
	char mask = 3;
	float Area_cubierta=0;
	float Area_solapada=0;
	float Fitness=0;
	int Costo_acumulado=0;
	int Antenas_usadas=0;
	float Penalizacion_FueraP =0;
	float Penalizacion_Total=0;
	float Penalizacion_Circulo=0;
	uchar gray;

	//decodificamos el cromosoma a gray
	for (int j = Nmax/8; j < Nmax*2+Nmax/8-1 ; j++)
	{
		gray=exGrayBinario(cromosoma[j]);
		x_g[j-4]=crDecodificar8(&gray,Xmin,Xmax);
		gray=exGrayBinario(cromosoma[j+1]);
		y_g[j-4]=crDecodificar8(&gray,Ymin,Ymax);
	}

	for (int i = 0; i < Nmax; i++)
	{	
		//revisamos si la antena esta presente o no 
		if(!leerBit(cromosoma,i))
			continue;
		//verificamos si esta dentro del poligono
		if(pnpoly(Nvert, X, Y, x_g[i], y_g[i]) == 1)
		{
			Antenas_usadas++;
			//leer potencia del cromosoma
			auxiliarP = (cromosoma[(int)((i*2)/8+Nmax*2+Nmax/8)] & (mask << (i*2%8)))>>(i*2%8);
			//calculamos el radio y lo guardamos en un arreglo
			radio[i] = calcular_radio(EIRP[auxiliarP]);

			//Calculamos el area cubierta y la acumulamos
			Area_cubierta = Area_cubierta+Pi*pow(radio[i],2);
			//Calculamos el costo y lo acumulamos
			Costo_acumulado = Costo_acumulado + Costos[auxiliarP];
			
			Penalizacion_Circulo = Penalizacion_Circulo + Circulo_Fuera_P(radio[i],x_g[i],y_g[i]);

			for (int j = 0; j < i ; j++)
			{
				//Verifico que el radio no es nulo
				if(radio[j]==calcular_radio(EIRP[0]) || radio[j]==calcular_radio(EIRP[1]) || radio[j]==calcular_radio(EIRP[2]) || radio[j]==calcular_radio(EIRP[3]))
				{
					//Verifico si estan solapadas
					if( ( pow(x_g[i]-x_g[j],2) + pow(y_g[i]-y_g[j],2) ) <= pow(radio[i]+radio[j],2) )
					{
						//verifico si esta una dentro de otra
						if( sqrt( pow(x_g[i]-x_g[j],2) + pow(y_g[i]-y_g[j],2) ) <= fabs(radio[i]-radio[j]) )
						{
							//calculo el area con solapada si estan una dentro de la otra
							if(radio[i]>=radio[j])
								Area_solapada=Area_solapada + Pi*pow(radio[j],2);
							//calculo el area con solapada si estan una dentro de la otra
							else
								Area_solapada=Area_solapada + Pi*pow(radio[i],2);
						}
						//calculo el area solapada si no estan una dentro de otra
						else
						{
							Area_solapada = Area_solapada + Overlap(x_g[i],y_g[i],radio[i],x_g[j],y_g[j],radio[j]);
						}
					}
				}
			}
		}
		//si esta fuera del poligono calculamos la distancia minima hasta el
		else
		{
			Penalizacion_FueraP = Penalizacion_FueraP + distancia_pol(Nvert,X,Y,x_g[i],y_g[i])/100;
		}
	}

	if(Area_solapada > Area_cubierta*0.1)
		Penalizacion_Total = Penalizacion_Total + (Area_solapada - Area_cubierta*0.1)*10;
	
	if(Costo_acumulado > PRESUPUESTO)
		Penalizacion_Total =Penalizacion_Total + (Costo_acumulado - PRESUPUESTO)*10;
	
	if(Antenas_usadas != 0)
		Penalizacion_Circulo= Penalizacion_Circulo/(Antenas_usadas);
	else
		Penalizacion_Circulo= Penalizacion_Circulo/(1+Antenas_usadas);

	Area_cubierta = Area_cubierta - Area_solapada;

	Penalizacion_Total= Penalizacion_Total + (float)Antenas_usadas/(float)(1+Nmax) + Area_total/(1+Area_cubierta) + 
						Penalizacion_FueraP + Penalizacion_Circulo ;

	Fitness = 100/(1+Penalizacion_Total);
	return Fitness;
}

/*
 ** ----------------------------------------------------------------------------
 **     Nombre      : 
 **     Función     : 
 **     Parámetros  : 
 **     Retorna     : 
 ** ----------------------------------------------------------------------------
*/
void mostrar(uchar* cromosoma)
{
    float x_g[Nmax];
	float y_g[Nmax];
	float radio[Nmax]={};
	int auxiliarP;
	char mask = 3;
	float Area_cubierta=0;
	float Area_solapada=0;
	int Costo_acumulado=0;
	int Antenas_usadas=0;
	float Penalizacion_FueraP =0;
	float Penalizacion_Total=0;
	float Penalizacion_Circulo=0;
	uchar gray;

	//decodificamos el cromosoma a gray
	for (int j = Nmax/8; j < Nmax*2+Nmax/8-1 ; j++)
	{
		gray=exGrayBinario(cromosoma[j]);
		x_g[j-4]=crDecodificar8(&gray,Xmin,Xmax);
		gray=exGrayBinario(cromosoma[j+1]);
		y_g[j-4]=crDecodificar8(&gray,Ymin,Ymax);
	}

	for (int i = 0; i < Nmax; i++)
	{	
		//revisamos si la antena esta presente o no 
		if(!leerBit(cromosoma,i))
			continue;
		//verificamos si esta dentro del poligono
		if(pnpoly(Nvert, X, Y, x_g[i], y_g[i]) == 1)
		{
			Antenas_usadas++;
			//leer potencia del cromosoma
			auxiliarP = (cromosoma[(int)((i*2)/8+Nmax*2+Nmax/8)] & (mask << (i*2%8)))>>(i*2%8);
			//calculamos el radio y lo guardamos en un arreglo
			radio[i] = calcular_radio(EIRP[auxiliarP]);
			//Calculamos el area cubierta y la acumulamos
			Area_cubierta = Area_cubierta+Pi*pow(radio[i],2);
			//Calculamos el costo y lo acumulamos
			Costo_acumulado = Costo_acumulado + Costos[auxiliarP];
		
			Penalizacion_Circulo = Penalizacion_Circulo + Circulo_Fuera_P(radio[i],x_g[i],y_g[i]);
			
			for (int j = 0; j < i ; j++)
			{
				//Verifico que el radio no es nulo
				if(radio[j]==calcular_radio(EIRP[0]) || radio[j]==calcular_radio(EIRP[1]) || radio[j]==calcular_radio(EIRP[2]) || radio[j]==calcular_radio(EIRP[3]))
				{  
					//Verifico si estan solapadas
					if( ( pow(x_g[i]-x_g[j],2) + pow(y_g[i]-y_g[j],2) ) <= pow(radio[i]+radio[j],2) )
					{	
						//verifico si esta una dentro de otra
						if( sqrt( pow(x_g[i]-x_g[j],2) + pow(y_g[i]-y_g[j],2) ) <= fabs(radio[i]-radio[j]) )
						{	
							//calculo el area con solapada si estan una dentro de la otra
							if(radio[i]>=radio[j]) 
							{
								Area_solapada=Area_solapada + Pi*pow(radio[j],2);
							}
							//calculo el area con solapada si estan una dentro dentro la otra
							else 
							{
								Area_solapada=Area_solapada + Pi*pow(radio[i],2);
							}
						}
						//calculo el area solapada si no estan una dentro de otra
						else
						{	
							Area_solapada = Area_solapada + Overlap(x_g[i],y_g[i],radio[i],x_g[j],y_g[j],radio[j]);
						}
					}
				}
			}
			fprintf(Coordenadas, "%.2lf\t%.2lf\t%.2lf\t%.2lf\t%d\t%d\t%d\n",x_g[i],y_g[i],2*radio[i],2*radio[i],EIRP[auxiliarP],Antenas_usadas, Costo_acumulado);
		}
		else
		{
			Penalizacion_FueraP = Penalizacion_FueraP + distancia_pol(Nvert,X,Y,x_g[i],y_g[i])/100;
		}
	}

}
