// OTM.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <iostream>
#include <math.h> 
#include <stdlib.h>


void printMatrix(long double **T)
{
	for (int x = 0; x < 2; x++)
	{
		for (int y = 0; y < 2; y++)  // loop for the three elements on the line
		{
			std::cout << T[x][y] << " ";  // display the current element out of the array
		}
		std::cout << std::endl;  // when the inner loop is done, go to a new line
	}
}

// MÉTODOS RELACIONADOS A FUNÇÃO DO GRUPO

long double func(long double x1, long double x2)
{
	// Função do Grupo
	long double func = log(2 * pow(x1, 2) - pow(x1, 4) + (pow(x1, 6) / 6) + x1 * x2 + pow(x2, 2));
	return func;
}

long double* gradiente(long double x1, long double x2)
{
	// Calculo do gradiente da função do Grupo
	long double dx1 = (4 * x1 - 4 * pow(x1, 3) + pow(x1, 5) + x2) / (2 * pow(x1, 2) - pow(x1, 4) + (pow(x1, 6) / 6) + x1 * x2 + pow(x2, 2));
	long double dx2 = (x1 + 2 * x2) / (2 * pow(x1, 2) - pow(x1, 4) + (pow(x1, 6) / 6) + x1 * x2 + pow(x2, 2));
	long double * gradiente = new long double[2];
	gradiente[0] = dx1;
	gradiente[1] = dx2;
	return gradiente;
}


long double** hessiana(long double x1, long double x2)
{
	// Calculo da hessiana da função do Grupo
	long double dx1x1 = (((4 - 12 * pow(x1, 2) + 5 * pow(x1, 4)) / (2 * pow(x1, 2) - pow(x1, 4) + (pow(x1, 6) / 6) + x1 * x2 + pow(x2, 2))) - pow((4 * x1 - 4 * pow(x1, 3) + pow(x1, 5) + x2) / (2 * pow(x1, 2) - pow(x1, 4) + (pow(x1, 6) / 6) + x1 * x2 + pow(x2, 2)), 2));
	long double dx1x2 = ((1 / (2 * pow(x1, 2) - pow(x1, 4) + (pow(x1, 6) / 6) + x1 * x2 + pow(x2, 2))) - (((4 * x1 - 4 * pow(x1, 3) + pow(x1, 5) + x2)*(x1 + 2 * x2)) / (pow((2 * pow(x1, 2) - pow(x1, 4) + (pow(x1, 6) / 6) + x1 * x2 + pow(x2, 2)), 2))));
	long double dx2x2 = ((2 / (2 * pow(x1, 2) - pow(x1, 4) + (pow(x1, 6) / 6) + x1 * x2 + pow(x2, 2))) - pow((x1 + 2 * x2) / (2 * pow(x1, 2) - pow(x1, 4) + (pow(x1, 6) / 6) + x1 * x2 + pow(x2, 2)), 2));
	long double ** hessiana = new long double*[2];
	for (int i = 0; i < 2 ; ++i) 
	{
		hessiana[i] = new long double[2];
	}

	hessiana[0][0] = dx1x1;
	hessiana[0][1] = dx1x2;
	hessiana[1][0] = dx1x2;
	hessiana[1][1] = dx2x2;
	
	return hessiana;

}

long double phi(long double x1, long double x2, long double* d, long double t)
{
	// Calcula o valor da função phi, que corresponde a f(x + td) = f(x1 + td[0], x2 + td[1])

	long double phi_x1 = x1 + d[0] * t; 
	long double phi_x2 = x2 + d[1] * t;
	long double phi = func(phi_x1, phi_x2);
	return phi;

}


// MÉTODOS DE BUSCA

long double secao_aurea(long double ro,long double epsilon,  long double x1, long double x2, long double* d)
{
	//Obtenção do intervalo [a,b], nessa etapa vamos "andando" com o intervalo [a,b] até sair do intervalo de descida da função
	//Como para esse método a função deve ser unimodal, temos que e estritamente decrescente em [0,t1] e estritamente crescente em [t2, ∞)
	long double a = 0;
	long double s = ro; 
	long double b = 2 * ro;
	while (phi(x1, x2, d, b) < phi(x1, x2, d, s)) // phi(b)< phi(s), no dado ponto e com dada direção
	{
		a = s;
		s = b;
		b = 2 * b;
	}

	//Obtenção de t*, que é mínimo exato da função phi, no intervalo [a,b]
	//Temos que  vantagem de não usar intervalos de tamanhos iguais, quando cortamos um, estaremos cortando 38%, ao invés de 33%, também só é necessário recalcular um ponto
	long double u = a + ((3 - sqrt(5)) / 2) *(b - a); // u ponto médio mais próximo de a
	long double v = a + ((sqrt(5) - 1) / 2) * (b - a); // v  ponto médio mais distante de a
	while ((b-a)> epsilon) // eenquanto a diferença entre ambos for maior que um limite
	{
		if (phi(x1, x2, d, u) < phi(x1, x2, d, v)) // phi(u)<phi(v), logo o trecho [v,b] pode ser descartado, pois não há como algo nesse intervalo ser menor que phi(u)
		{
			b = v;
			v = u;
			u = a + ((3 - sqrt(5)) / 2)*(b - a); //anda com intervalo para esquerda
		}
		else // phi(v)<phi(u), logo o trecho [a,u] pode ser descartado, pois não há como algo nesse intervalo ser menor que phi(v)
		{
			a = u;
			u = v;
			v = a + ((sqrt(5) - 1) / 2) * (b - a); //anda com intervalo para direita
		}
	}
	long double t_min = (u + v) / 2;
	return t_min;
}

long double computeGrad_x_D(long double x1, long double x2, long double * d)
{
  	long double * grad = gradiente(x1,x2);
	long double result = grad[0]*d[0] + grad[1]*d[1]; //multiplicação do gradiente pela direção
	return result;
}


long double armijo(long double N,long double gamma,  long double x1, long double x2, long double* d)
{
	//Dado a funcao F, um ponto X, e uma direcao de descida d queremos achar t tal que t>0 e f(X + td) < f(X). Se d é uma direcao de descida, existe delta tal que para todo t pertencente [0,delta] a funcao objetivo decresce.
	//O algoritmo encontra um minimizador dentro de uma tolerância se φ for unimodal. O método de Armijo procura uma boa redução da função ao longo da direção, sem tentar minimizá-la.
	long double t = 1.0;
	while(phi(x1,x2,d,t) > (func(x1,x2) + N * t * computeGrad_x_D(x1,x2,d) ) ) 
	{
		//A condiçãoo acima indica que queremos uma redução proporcional ao tamanho do passo
		t *= gamma;
	}
	return t;
}


// IMPLEMENTAÇÃO DOS ALGORITMOS PARA ENCONTRAR PONTO MÍNIMO

long double algoritmo_gradiente(long double x1, long double x2, char busca)
{
	int k = 0; //k é a interação em que estamos

	// A ideia é dentre as direçõoes de descida, -gradiente é a direção com o maior decréscimo.
	// Duas direções conhecidas são ortogonais
	long double* d = gradiente(x1, x2);
	d[0] = -d[0];
	d[1] = -d[1];

	//DEFINIR AQUI CONSTANTES UTILIZADAS
	long double ro = 0.00001; // >0
	long double epsilon = 0.00001;   // >0
	long double N = 0.25; //(0,1)
	long double gamma = 0.8; //(0,1)

	long double t;
	while (not(abs(d[0]) < 0.00001 and abs(d[1]) < 0.00001) and k < 100) // Critérios utilizados: gradiente igual a zero (como não ficava igual a zero, demos uma margem de erro) e limite de interações
	{
		if (busca == 'S') //Utiliza método de busca da seção Aurea
		{
			t = secao_aurea(ro, epsilon, x1, x2, d);
		}
		else //Utiliza método de busca de Armijo
		{
			t = armijo(N, gamma, x1, x2, d);
		}
		// Calcula novo x
		x1 = x1 + t * d[0];
		x2 = x2 + t * d[1];
		d = gradiente(x1, x2);
		d[0] = -d[0];
		d[1] = -d[1];
		k++;
	}
	std::cout << "Iter: " << k << "\nx1: " << x1 << "\nx2: " << x2 << "\nValue: " << func(x1, x2);
	return x1, x2, k;
}



long double newton(long double x1, long double x2, char busca)
{
	// O método de Newton para resolver um sistema utiliza o seu polinômio de Taylor de primeira ordem parar aproximar o sistema.

	int k = 0; //k é a interação em que estamos
	long double* d = new long double[2];
	long double**i_H = new long double*[2]; //matriz inversa da hessiana
	for (int i = 0; i < 2; ++i)
	{
		i_H[i] = new long double[2];
	}

	//DEFINIR AQUI CONSTANTES UTILIZADAS
	long double ro = 0.00001; // >0
	long double epsilon = 0.00001;   // >0
	long double N = 0.25; //(0,1)
	long double gamma = 0.8; //(0,1)

	long double ** H = hessiana(x1, x2);
	printMatrix(H);

	long double* grad = gradiente(x1, x2);
	//std::cout << grad[0] << ", " << grad[1];
	long double t;
	while (not(abs(grad[0]) < 0.00001 and abs(grad[1]) < 0.00001) and k < 100) // Critérios utilizados: gradiente igual a zero (como não ficava igual a zero, demos uma margem de erro) e limite de interações
	{
		//Invertendo a Hessiana usando a formula para matriz 2x2
		long double m =( 1 / (H[0][0] * H[1][1] - H[0][1] * H[1][0]));
		i_H[0][0] = m * H[1][1];
		i_H[0][1] = -m * H[1][0];
		i_H[1][0] = -m * H[0][1];
		i_H[1][1] = m * H[0][0];
		//printMatrix(i_H);

		//Encontrando o vetor direção igual a menos a inversa da hessiana vezes o gradiente
		d[0] = -(i_H[0][0] *grad[0]) - (i_H[0][1] * grad[1]);
		d[1] = -(i_H[1][0] * grad[0]) - (i_H[1][1] * grad[1]);
		//std::cout << d[0] << "," << d[1] << "\n";
		if (busca == 'S') //Utiliza método de busca da seção Aurea
		{
			t = secao_aurea(ro, epsilon, x1, x2, d);
		}
		else //Utiliza método de busca de Armijo
		{
			t = armijo(N, gamma, x1, x2, d);
		}
		// Calcula novo x
		x1 = x1 + t * d[0];
		x2 = x2 + t * d[1];
		grad = gradiente(x1, x2);
		H = hessiana(x1, x2);	
		//std::cout <<  k << ", " << x1 << ", " << x2 << ","<<t<<"\n";
		k++;
	}

	std::cout <<"Iter: "<<k<< "\nx1: "<< x1 << "\nx2: " << x2 << "\nValue: "<<func(x1,x2) ;
	return x1, x2, k;
}

long double **getXXtranspose(long double *X)
{
	long double**result = new long double*[2];
	for (int i = 0; i < 2; ++i)
	{
		result[i] = new long double[2];
	}
	result[0][0] = X[0] * X[0];
	result[0][1] = X[0] * X[1];
	result[1][0] = X[0] * X[1];
	result[1][1] = X[1] * X[1];
	return result;
}

long double **getABtranspose(long double *A, long double *B)
{
	long double**result = new long double*[2];
	for (int i = 0; i < 2; ++i)
	{
		result[i] = new long double[2];
	}
	result[0][0] = A[0] * B[0];
	result[0][1] = A[0] * B[1];
	result[1][0] = A[1] * B[0];
	result[1][1] = A[1] * B[1];
	return result;
}

long double** getTranspose(long double **matrix, int m, int n)
{
	int i, j;
	long double**Transpose = (long double**) malloc(n * sizeof(long double*));
	for (i = 0; i < n; i++)
		Transpose[i] = (long double*) malloc(m * sizeof(long double));
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			Transpose[i][j] = matrix[j][i];
	return Transpose;
}

//Matriz 1 possui m linhas e n colunas e Matriz2 possui o colunas
void MultiplicarMatrizes(long double **matriz1, long double **matriz2, long double **matrizResultado, int m, int n, int o)
{
	int linha1, coluna1, linha2, coluna2;
	for (linha1 = 0; linha1 < m; linha1++)
	{
		for (coluna2 = 0; coluna2 < o; coluna2++)
		{
			matrizResultado[linha1][coluna2] = 0;
			linha2 = 0;
			for (coluna1 = 0; coluna1 < n; coluna1++)
			{
				matrizResultado[linha1][coluna2] += matriz1[linha1][coluna1] * matriz2[linha2][coluna2];
				linha2++;
			}
		}
	}
}

void updateHessianInverse_DFP(long double **H_k, long double *P_k, long double *Q_k)
{
	long double **P_Ptranspose = getXXtranspose(P_k);
	long double **Q_Qtranspose = getXXtranspose(Q_k);
	long double constante1 = P_k[0] * Q_k[0] + P_k[1] * Q_k[1];
	long double *qTranspose_H = new long double[2];
	qTranspose_H[0] = Q_k[0] * H_k[0][0] + Q_k[1] * H_k[1][0];
	qTranspose_H[1] = Q_k[0] * H_k[0][1] + Q_k[1] * H_k[1][1];
	long double constante2 = qTranspose_H[0] * Q_k[0] + qTranspose_H[1] * Q_k[1];

	long double **matriz1 = new long double*[2];//Matriz do primeiro termo
	for (int i = 0; i < 2; ++i)
	{
		matriz1[i] = new long double[2];
	}

	matriz1[0][0] = P_Ptranspose[0][0] / constante1;
	matriz1[0][1] = P_Ptranspose[0][1] / constante1;
	matriz1[1][0] = P_Ptranspose[1][0] / constante1;
	matriz1[1][1] = P_Ptranspose[1][1] / constante1;

	long double **matriz2 = new long double*[2];//Matriz do segundo termo
	for (int i = 0; i < 2; ++i)
	{
		matriz2[i] = new long double[2];
	}

	long double **matrizIntermediaria = new long double*[2];//Matriz obtida fazendo MultiplicarMatrizes(H_k,Q_Qtranspose)
	for (int i = 0; i < 2; ++i)
	{
		matrizIntermediaria[i] = new long double[2];
	}

	MultiplicarMatrizes(H_k, Q_Qtranspose, matrizIntermediaria, 2, 2, 2);
	MultiplicarMatrizes(matrizIntermediaria, getTranspose(H_k, 2, 2), matriz2, 2, 2, 2);

	matriz2[0][0] /= constante2;
	matriz2[0][1] /= constante2;
	matriz2[1][0] /= constante2;
	matriz2[1][1] /= constante2;

	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			H_k[i][j] += (matriz1[i][j] - matriz2[i][j]);
		}
	}
}

void updateHessianInverse_BFGS(long double **H_k, long double *P_k, long double *Q_k)
{
	long double **P_Ptranspose = getXXtranspose(P_k);
	long double **Q_Qtranspose = getXXtranspose(Q_k);
	long double **P_Qtranspose = getABtranspose(P_k, Q_k);
	long double **Q_Ptranspose = getABtranspose(Q_k, P_k);
	long double constante1 = P_k[0] * Q_k[0] + P_k[1] * Q_k[1];
	long double *qTranspose_H = new long double[2];
	qTranspose_H[0] = Q_k[0] * H_k[0][0] + Q_k[1] * H_k[1][0];
	qTranspose_H[1] = Q_k[0] * H_k[0][1] + Q_k[1] * H_k[1][1];
	long double constante2 = qTranspose_H[0] * Q_k[0] + qTranspose_H[1] * Q_k[1];
	long double constante3 = 1 + constante2 / constante1;

	long double **matriz1 = new long double*[2];//Matriz do primeiro termo
	for (int i = 0; i < 2; ++i)
	{
		matriz1[i] = new long double[2];
	}

	matriz1[0][0] = P_Ptranspose[0][0] * constante3 / constante1;
	matriz1[0][1] = P_Ptranspose[0][1] * constante3 / constante1;
	matriz1[1][0] = P_Ptranspose[1][0] * constante3 / constante1;
	matriz1[1][1] = P_Ptranspose[1][1] * constante3 / constante1;

	long double **matriz2 = new long double*[2];//Matriz do segundo termo
	for (int i = 0; i < 2; ++i)
	{
		matriz2[i] = new long double[2];
	}

	long double **matrizIntermediaria1 = new long double*[2];//Matriz obtida fazendo MultiplicarMatrizes(P_Qtranspose,H_k)
	long double **matrizIntermediaria2 = new long double*[2];//Matriz obtida fazendo MultiplicarMatrizes(H_k ,Q_Ptranspose);
	for (int i = 0; i < 2; ++i)
	{
		matrizIntermediaria1[i] = new long double[2];
		matrizIntermediaria2[i] = new long double[2];
	}

	MultiplicarMatrizes(P_Qtranspose, H_k, matrizIntermediaria1, 2, 2, 2);
	MultiplicarMatrizes(H_k, Q_Ptranspose, matrizIntermediaria2, 2, 2, 2);


	matriz2[0][0] = (matrizIntermediaria1[0][0] + matrizIntermediaria2[0][0]) / constante1;
	matriz2[0][1] = (matrizIntermediaria1[0][1] + matrizIntermediaria2[0][1]) / constante1;
	matriz2[1][0] = (matrizIntermediaria1[1][0] + matrizIntermediaria2[1][0]) / constante1;
	matriz2[1][1] = (matrizIntermediaria1[1][1] + matrizIntermediaria2[1][1]) / constante1;

	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			H_k[i][j] += (matriz1[i][j] - matriz2[i][j]);
		}
	}
}

//Metodo do Quasi Newton: A ideia é construir aproximações para a Hessiana da função objetivo ao longo das iterações.
//A direção de busca é dado por d_k = −H_k ∇f (x_k )
//É possı́vel estimar a hessiana de f (x) partindo da equação já deduzida: 
//∇^2 f (x_k )(x_k+1 − x_k ) = ∇f (x_k+1 ) − ∇f (x_k )
//

long double QuaseNewton(long double x1, long double x2, char busca, char method)
{
	int k = 0; //k é a interação em que estamos
	//O metodo de quasinewton utiliza dois conjuntos de vetores p e q para realizar a estimativa da inversa da Hessiana.
	long double** p = new long double*[2];//p_k = x_k+1 − x_k, diferença dos pnts
	long double** q = new long double*[2];//q_k = ∇f ( x_k+1 ) − ∇f ( x_k ), diferença dos gradientes
	long double* d = new long double[2];
	long double**H_k = new long double*[2]; //estimativa da inversa da matriz hessiana inicializada com a identidade
	for (int i = 0; i < 2; ++i)
	{
		H_k[i] = new long double[2];
		p[i] = new long double[2];
		q[i] = new long double[2];
	}

	H_k[0][0] = 1;
	H_k[0][1] = 0;
	H_k[1][0] = 0;
	H_k[1][1] = 1;


	//DEFINIR AQUI CONSTANTES UTILIZADAS
	long double ro = 0.00001; // >0
	long double epsilon = 0.00001;   // >0
	long double N = 0.001; //(0,1)
	long double gamma = 0.001; //(0,1)

	long double* grad = gradiente(x1, x2);

	long double t;
	long double * p_k = new long double[2];
	//p_k[0] = 1;
	//p_k[1] = 1;
	//and calcularResiduo(p_k,2,2) > 0.0001
	std::cout << not(abs(grad[0]) < 0.00001 and abs(grad[1]) < 0.00001) << (k < 100) << std::endl;
	while (k < 100) // Critérios utilizados: gradiente igual a zero e limite de interações
	{
		//Encontrando o vetor direção igual a menos a inversa da hessiana vezes o gradiente
		d[0] = H_k[0][0] * grad[0] + H_k[0][1] * grad[1];
		d[1] = H_k[1][0] * grad[0] + H_k[1][1] * grad[1];
		d[0] *= -1;
		d[1] *= -1;

		if (busca == 'S') //Utiliza método de busca da seção Aurea
		{
			t = secao_aurea(ro, epsilon, x1, x2, d);
		}
		else //Utiliza método de busca de Armijo
		{
			t = armijo(N, gamma, x1, x2, d);
			if (k == 0)
				std::cout << "valor de t: " << t << std::endl;
		}

		long double *p_k = new long double[2];
		long double *q_k = new long double[2];
		p_k[0] = t * d[0];
		p_k[1] = t * d[1];
		q_k[0] = -grad[0];
		q_k[1] = -grad[1];
		// Calcula novo x
		long double old_x1 = x1;
		long double old_x2 = x2;

		x1 = x1 + t * d[0];
		x2 = x2 + t * d[1];

		//if(func(x1,x2) > func(old_x1,old_x2) )
		//	std::cout<<"ERROR"<<std::endl;
		grad = gradiente(x1, x2);
		q_k[0] += grad[0];
		q_k[1] += grad[1];
		if (k < 2)//Sao necessarios ao menos 3 pontos para obter a Hessiana de dimensao 2x2
		{
			for (int j = 0; j < 2; j++)//insere o p_k e q_k na matriz p e q
			{
				p[j][k] = p_k[j];
				q[j][k] = q_k[j];
			}
			if (method == 'D')
				updateHessianInverse_DFP(H_k, p_k, q_k);
			else
				updateHessianInverse_BFGS(H_k, p_k, q_k);
		}
		else
		{
			if (k == 2)//podemos entao calcular a inversa da hessiana sem precisar atualiza-la mais
			{
				long double **q_inverse = new long double*[2];
				for (int i = 0; i < 2; ++i)
				{
					q_inverse[i] = new long double[2];
				}
				long double m = 1 / (q[0][0] * q[1][1] - q[0][1] * q[1][0]);
				q_inverse[0][0] = m * q[1][1];
				q_inverse[0][1] = -m * q[1][0];
				q_inverse[1][0] = -m * q[0][1];
				q_inverse[1][1] = m * q[0][0];
				//std::cout << "Q(inversa) :" << std::endl;
				//printMatrix(q_inverse);
				//std::cout << "P :" << std::endl;
				//printMatrix(p);
				MultiplicarMatrizes(p, q_inverse, H_k, 2, 2, 2);
				//std::cout << "Hessiana(inversa) :" << std::endl;
				//printMatrix(H_k);
			}
		}
		k++;
	}

	//std::cout << x1 << ", " << x2 << ", " << k << std::endl;
	//std::cout << "Opt value: " << func(x1, x2) << std::endl;
	//std::cout << "Grad: [" << grad[0] << ", " << grad[1] << "]" << std::endl;
	std::cout << "Iter: " << k << "\nx1: " << x1 << "\nx2: " << x2 << "\nValue: " << func(x1, x2);
	return x1, x2, k;
}


int main()
{
	//USAR 'S' PARA BUSCA POR SEÇÃO AUREA E 'A' PARA ARMIJO
	//algoritmo_gradiente(1, 1,'A');
	QuaseNewton(1, 1, 'S', 'B');
	//newton(0.01, 0.01, 'A');
	return 1;

}

