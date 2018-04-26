#include <iostream>
#include <cmath>
#define PI 3.14159265

using std::cout;
using std::endl;

//Se definen funciones
void cond_in(double *y, double x_min, double dx, int n, double A);
void evolucion(double *Y, double *y_minus, double *y_central, int n, double dx, double dt, double c);
void copia(double *Y, double *y_old, int n);
bool print(double*y, int n ,double A );

//Condicion inicial
void cond_in(double *y, double x_min, double dx, int n, double A){
  int i;
  double x;
  for(i=0;i<n;i++){
    x = x_min +  i * dx;
    y[i] = A*sin(x/(2*PI));
  }
}

//Evoluciona la funcion. Se utiliza un arreglo diferente para guardar los valores modificados
void evolucion(double *Y, double *y_minus, double *y_central, int n, double dx, double dt, double c){
  int i;
  double M = pow(c*dt/dx , 2);

  for(i=1;i<n-1;i++){   

    Y[i] = M*(y_central[i+1]+y_central[i-1])+y_central[i]*(2-2*M) - y_minus[i];   
 
  }
}

//actualiza los arreglos para seguir en el siguiente indice temporal
void copia(double *Y, double *y_central, double *y_minus, int n){
  int i;
  for(i=0;i<n;i++){
    y_minus[i]=y_central[i];
    y_central[i]=Y[i];    
  }
}

bool print(double*y, int n, double A){
  double max=0.0;
  for(int i=0; i<n;i++){
    if(y[i]>max){
      max=y[i];    
    }  
  }
  
  if(max == A/2 ){
    return true;   
  }else if(max = A){
    return true;
  }else{
    return false;
  }
}
 //Funcion main
int main(){
  //Definicion de variables
  double A = 1.0;
  double x_min=0.0, x_max=2*PI;
  double dt = 0.001, dx=0.002; 
  int Nx = (x_max -x_min)/dx + 1;
  double t=0.0;
  double c=0.1;  
  double M = pow(c*dt/dx , 2);
  double * y_minus = new double[Nx];
  double * y_central = new double[Nx];
  double * Y = new double[Nx];
  //inicializa la condicion inicial en u_old
  cond_in(y_central, x_min, dx, Nx, A);
  
  for(int j = 0;j<1000;j++){
    
    //en j=0 se tiene la condicion que u(+1) = u(-1) temporalmente
    if(j=0){
      for(int k=0;k<Nx;k++){
        cout<<t<<" "<<x_min+k*dx<<" "<<y_central[k]<<endl;
      }
      Y[j] = (M*(y_central[j+1]+y_central[j-1])+y_central[j]*(2-2*M))/2;
      copia(Y, y_central, y_minus, Nx);   
      t += dt;
      
    }else{
      if(print(y_central,Nx,A) == true){
        for(int k=0;k<Nx;k++){
          cout<<t<<" "<<x_min+k*dx<<" "<<y_central[k]<<endl;
        }
      }
      evolucion(Y, y_minus, y_central, Nx, dx, dt, c);    
      copia(Y, y_central, y_minus, Nx);    
      t += dt;    
    }
  }
  
  return 0;
}
