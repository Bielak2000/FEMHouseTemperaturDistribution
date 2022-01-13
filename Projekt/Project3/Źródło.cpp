#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string>

using namespace std;

double function1d(double x)
{
	return 5 * x * x + 3 * x + 6;
}

double function2d(double x, double y)
{
	return 5 * x * x * y * y + 6 * x * y + 6;
}

struct node
{
	double x, y;
	int BC;
};

struct element
{
	int ID[4];
	double H[4][4];
	double HBC[4][4];
	double P[4];
	double C[4][4];

	element()
	{
		for (int i = 0; i < 4; i++)
		{
			P[i] = 0.0;

			for (int j = 0; j < 4; j++)
			{
				HBC[i][j] = 0.0;
				C[i][j] = 0.0;
			}

		}
	}

	void clear()
	{
		for (int i = 0; i < 4; i++)
		{
			P[i] = 0.0;

			for (int j = 0; j < 4; j++)
			{
				HBC[i][j] = 0.0;
				C[i][j] = 0.0;
			}

		}
	}

	void showH()
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << H[i][j] << "\t";
			}
			cout << endl;
		}
		cout << endl;
	}
};

struct Grid
{
	double H;//wysokosc
	double B;//szerokosc
	int nH;//liczba wezlow do gory
	int nB;//liczba wezlow wszerz
	int nN;//liczba wezlow
	int nE;//liczba elementow
	node* nodes;
	element* elements;

	//konstruktor
	Grid(double height, double width, int y2, int x2)
	{
		H = height;
		B = width;
		nH = y2;
		nB = x2;
		nN = nH * nB;
		nE = (nH - 1) * (nB - 1);
		nodes = new node[nN];
		elements = new element[nE];
	}

	Grid(int nN1, int nE1)
	{
		nE = nE1;
		nN = nN1;
		nodes = new node[nN];
		elements = new element[nE];
		for (int i = 0; i < nN; i++)
		{
			nodes[i].BC = 0;
		}
	}

	void initalization_nodes()
	{
		double deltaY = (double)(H / (nH - 1));
		double deltaX = (double)(B / (nB - 1));
		int temp = 0;//zmienna pomocniczna
		for (int i = 0; i < nB; i++)
		{
			for (int j = 0; j < nH; j++)
			{
				temp = j + i * (nH - 1);

				nodes[i + temp].y = j * deltaY;
				nodes[i + temp].x = i * deltaX;
				if (nodes[i + temp].x == 0 || nodes[i + temp].x == B || nodes[i + temp].y == 0 || nodes[i + temp].y == H)
					nodes[i + temp].BC = 1;
				else
					nodes[i + temp].BC = 0;
			}
		}
	}

	void initialization_elements()
	{
		int temp = 1;//zmienna pomocnicza
		for (int i = 0; i < nE; i++)
		{
			if (temp % nH == 0)
				temp++;
			elements[i].ID[0] = temp;
			elements[i].ID[1] = temp + nH;
			elements[i].ID[2] = elements[i].ID[1] + 1;
			elements[i].ID[3] = elements[i].ID[0] + 1;
			temp++;

		}
	}

	void show_nodes()
	{
		cout << "Coordinates of the nodes: " << endl;
		for (int i = 0; i < nN; i++)
		{
			cout << i + 1 << ". x: " << nodes[i].x << ", y: " << nodes[i].y << ", BC: " << nodes[i].BC << "." << endl;
		}
	}

	void show_elements()
	{
		cout << "ID elements: " << endl;
		for (int i = 0; i < nE; i++)
		{
			cout << i + 1 << ". " << elements[i].ID[0] << ", " << elements[i].ID[1] << ", " << elements[i].ID[2] << ", " << elements[i].ID[3] << "." << endl;
		}
	}
};

class gauss
{
	double* nodes;
	double* A;

public:

	gauss(int n)
	{
		switch (n)
		{
		case 2:
			nodes = new double[2];
			nodes[0] = -(double(1 / sqrt(3)));
			nodes[1] = double(1 / sqrt(3));
			A = new double[2];
			A[0] = 1;
			A[1] = 1;
			break;
		case 3:
			nodes = new double[3];
			nodes[0] = -double(double(sqrt(3)) / double(sqrt(5)));
			nodes[1] = 0;
			nodes[2] = double(double(sqrt(3)) / double(sqrt(5)));
			A = new double[3];
			A[0] = double(double(5) / double(9));
			A[1] = double(double(8) / double(9));
			A[2] = double(double(5) / double(9));
			break;
		case 4:
			nodes = new double[4];
			nodes[0] = -0.861136;
			nodes[1] = -0.339981;
			nodes[2] = 0.339981;
			nodes[3] = 0.861136;
			A = new double[4];
			A[0] = 0.347855;
			A[1] = 0.652145;
			A[2] = 0.652145;
			A[3] = 0.347855;
			break;
		}
	}

	double getNode(int i)
	{
		return nodes[i];
	}

	double getA(int i)
	{
		return A[i];
	}
};

struct jakobian
{
	double** jak;
	double** jak_inv;
	double det_jak;

	jakobian(int n)
	{
		jak = new double* [n];
		jak_inv = new double* [n];
		for (int i = 0; i < n; i++)
		{
			jak[i] = new double[n];
			jak_inv[i] = new double[n];
		}

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				jak[i][j] = 0;
				jak_inv[i][j] = 0;
			}
	}

	void print()
	{
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				cout << jak[i][j] << "\t";
			cout << endl;
			cout << endl;
		}
	}
};

double methodGauss1D(int dimension)
{
	gauss g(dimension);
	double result = 0;
	for (int i = 0; i < dimension; i++)
	{
		result += g.getA(i) * function1d(g.getNode(i));
	}
	return result;
}

double methodGauss2D(int dimension)
{
	gauss g(dimension);
	double result = 0;
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			result += g.getA(i) * g.getA(j) * function2d(g.getNode(i), g.getNode(j));
		}
	}
	return result;
}

double dlugoscOdcinka(double a1, double a2)
{
	double d = (a2 - a1) / 2.0;
	return d;
}

double dlugoscOdcinka(double a1, double a2, double b1, double b2)
{
	double d = (sqrt(pow((a2 - a1), 2) + pow((b2 - b1), 2))) / 2.0;
	return d;
}

struct sciana
{
	double** tab;
	double* waga;
	double* wspKSI;
	double* wspETA;
};

struct Element4_2D
{
	double** arrayKSI;
	double** arrayETA;
	double** arrayNC;
	double* wagiDoC;
	sciana tablica[4];
	int count;
	Element4_2D(int n)
	{
		count = n;
		gauss g(n);
		int k = 0;
		int e = 0;
		int temp = 0;
		arrayKSI = new double* [pow(n, 2)];
		arrayETA = new double* [pow(n, 2)];
		arrayNC = new double* [pow(n, 2)];
		for (int i = 0; i < pow(n, 2); i++)
		{
			arrayKSI[i] = new double[4];
			arrayETA[i] = new double[4];
			arrayNC[i] = new double[4];
		}

		//arrayKSI i arrayETA
		for (int j = 0; j < pow(n, 2); j++)
		{
			if (j % n == 0 && j != 0)
			{
				k++;
				e = 0;
				temp++;
			}

			arrayKSI[j][0] = -(1.0) / (4.0) * (1.0 - g.getNode(k));
			arrayKSI[j][1] = (1.0) / (4.0) * (1.0 - g.getNode(k));
			arrayKSI[j][2] = (1.0) / (4.0) * (1.0 + g.getNode(k));
			arrayKSI[j][3] = -(1.0) / (4.0) * (1.0 + g.getNode(k));

			arrayETA[j][0] = -(1.0) / (4.0) * (1.0 - g.getNode(e));
			arrayETA[j][1] = -(1.0) / (4.0) * (1.0 + g.getNode(e));
			arrayETA[j][2] = (1.0) / (4.0) * (1.0 + g.getNode(e));
			arrayETA[j][3] = (1.0) / (4.0) * (1.0 - g.getNode(e));
			e++;
		}

		k = 0;
		e = 0;
		temp = 0;

		//warunek brzegowy
		for (int i = 0; i < 4; i++)
		{
			tablica[i].wspKSI = new double[n];
			tablica[i].wspETA = new double[n];
			tablica[i].waga = new double[n];
			tablica[i].tab = new double* [n];
			for (int j = 0; j < n; j++)
			{
				tablica[i].tab[j] = new double[4];
				for (int x = 0; x < 4; x++)
					tablica[i].tab[j][x] = 0;
			}
		}

		//sciana 1
		for (int i = 0; i < count; i++)
		{
			tablica[0].wspKSI[i] = g.getNode(i);
			tablica[0].wspETA[i] = -1;
			tablica[0].waga[i] = g.getA(i);
		}

		//sciana 2
		for (int i = 0; i < count; i++)
		{
			tablica[1].wspKSI[i] = 1;
			tablica[1].wspETA[i] = g.getNode(i);
			tablica[1].waga[i] = g.getA(i);
		}

		//sciana 3
		for (int i = 0; i < count; i++)
		{
			tablica[2].wspKSI[i] = g.getNode(i);
			tablica[2].wspETA[i] = 1;
			tablica[2].waga[i] = g.getA(i);
		}

		//sciana 4
		for (int i = 0; i < count; i++)
		{
			tablica[3].wspKSI[i] = -1;
			tablica[3].wspETA[i] = g.getNode(i);
			tablica[3].waga[i] = g.getA(i);
		}

		//funkcje ksztaltu
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < n; j++)
			{
				tablica[i].tab[j][0] = 1.0 / 4.0 * (1.0 - tablica[i].wspKSI[j]) * (1.0 - tablica[i].wspETA[j]);
				tablica[i].tab[j][1] = 1.0 / 4.0 * (1.0 + tablica[i].wspKSI[j]) * (1.0 - tablica[i].wspETA[j]);
				tablica[i].tab[j][2] = 1.0 / 4.0 * (1.0 + tablica[i].wspKSI[j]) * (1.0 + tablica[i].wspETA[j]);
				tablica[i].tab[j][3] = 1.0 / 4.0 * (1.0 - tablica[i].wspKSI[j]) * (1.0 + tablica[i].wspETA[j]);
			}
		}

		//Wartosci funkcji ksztaltu w punktach calkowania - dla macierzy C
		for (int i = 0; i < pow(n, 2); i++)
		{
			if (i % n == 0 && i != 0)
			{
				k++;
				e = 0;
				temp++;
			}

			arrayNC[i][0] = 1.0 / 4.0 * (1 - g.getNode(e)) * (1 - g.getNode(k));
			arrayNC[i][1] = 1.0 / 4.0 * (1 + g.getNode(e)) * (1 - g.getNode(k));
			arrayNC[i][2] = 1.0 / 4.0 * (1 + g.getNode(e)) * (1 + g.getNode(k));
			arrayNC[i][3] = 1.0 / 4.0 * (1 - g.getNode(e)) * (1 + g.getNode(k));

			e++;
		}

		wagiDoC = new double[count];
		for (int i = 0; i < count; i++)
		{
			wagiDoC[i] = g.getA(i);
		}
	}
};

void function_jakobian(int i, int j, jakobian* j1, Element4_2D e, Grid g)
{
	double dXdKsi = e.arrayKSI[j][0] * g.nodes[g.elements[i].ID[0] - 1].x + e.arrayKSI[j][1] * g.nodes[g.elements[i].ID[1] - 1].x + e.arrayKSI[j][2] * g.nodes[g.elements[i].ID[2] - 1].x + e.arrayKSI[j][3] * g.nodes[g.elements[i].ID[3] - 1].x;
	double dXdEta = e.arrayETA[j][0] * g.nodes[g.elements[i].ID[0] - 1].x + e.arrayETA[j][1] * g.nodes[g.elements[i].ID[1] - 1].x + e.arrayETA[j][2] * g.nodes[g.elements[i].ID[2] - 1].x + e.arrayETA[j][3] * g.nodes[g.elements[i].ID[3] - 1].x;
	double dYdKsi = e.arrayKSI[j][0] * g.nodes[g.elements[i].ID[0] - 1].y + e.arrayKSI[j][1] * g.nodes[g.elements[i].ID[1] - 1].y + e.arrayKSI[j][2] * g.nodes[g.elements[i].ID[2] - 1].y + e.arrayKSI[j][3] * g.nodes[g.elements[i].ID[3] - 1].y;
	double dYdEta = e.arrayETA[j][0] * g.nodes[g.elements[i].ID[0] - 1].y + e.arrayETA[j][1] * g.nodes[g.elements[i].ID[1] - 1].y + e.arrayETA[j][2] * g.nodes[g.elements[i].ID[2] - 1].y + e.arrayETA[j][3] * g.nodes[g.elements[i].ID[3] - 1].y;

	j1->jak[0][0] = dXdKsi;
	j1->jak[0][1] = dYdKsi;
	j1->jak[1][0] = dXdEta;
	j1->jak[1][1] = dYdEta;

	j1->det_jak = (dXdKsi * dYdEta) - (dYdKsi * dXdEta);

	j1->jak_inv[0][0] = dYdEta;
	j1->jak_inv[0][1] = -(dYdKsi);
	j1->jak_inv[1][0] = -(dXdEta);
	j1->jak_inv[1][1] = dXdKsi;
}

struct array2D
{
	double** dNdx;
	double** dNdy;

	array2D(int n)
	{
		dNdx = new double* [pow(n, 2)];
		dNdy = new double* [pow(n, 2)];
		for (int i = 0; i < pow(n, 2); i++)
		{
			dNdx[i] = new double[4];
			dNdy[i] = new double[4];
			for (int j = 0; j < 4; j++)
			{
				dNdx[i][j] = 0;
				dNdy[i][j] = 0;
			}
		}
	}

	void print()
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
				cout << dNdx[i][j] << "\t";
			cout << endl;
		}
	}
};

struct arrayHtemp
{
	double Hp[4][4];

	arrayHtemp()
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Hp[i][j] = 0;
			}
		}
	}

	void clear()
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Hp[i][j] = 0;
			}
		}
	}
};

void dNdx_dNdy(jakobian ja, int j, array2D& a, Element4_2D e)
{
	for (int i = 0; i < 4; i++)
	{
		a.dNdx[j][i] = (1.0 / ja.det_jak) * (ja.jak_inv[0][0] * e.arrayKSI[j][i] + ja.jak_inv[0][1] * e.arrayETA[j][i]);
		a.dNdy[j][i] = (1.0 / ja.det_jak) * (ja.jak_inv[1][1] * e.arrayETA[j][i] + ja.jak_inv[1][0] * e.arrayKSI[j][i]);
	}
}

void H(double k, array2D& a, arrayHtemp& H1, jakobian ja, int j, double waga)
{
	for (int i = 0; i < 4; i++)
	{
		for (int x = 0; x < 4; x++)
		{
			H1.Hp[i][x] = k * ja.det_jak * waga * ((a.dNdx[j][i] * a.dNdx[j][x]) + (a.dNdy[j][i] * a.dNdy[j][x]));
		}
	}
}

void HBC(double a, Element4_2D e, arrayHtemp& temp, double detJ, int numerSciany, int numerCalkowania)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			temp.Hp[i][j] = a * detJ * (e.tablica[numerSciany].waga[numerCalkowania] * (e.tablica[numerSciany].tab[numerCalkowania][i] * e.tablica[numerSciany].tab[numerCalkowania][j]));
		}
	}
}

void P(double a, Element4_2D e, double* temp, double detJ, int numerSciany, int numerCalkowania, double Totcz)
{
	for (int i = 0; i < 4; i++)
	{
		temp[i] = a * detJ * Totcz * e.tablica[numerSciany].waga[numerCalkowania] * e.tablica[numerSciany].tab[numerCalkowania][i];
	}
}

void C(double cieplo, double gestosc, Element4_2D e, arrayHtemp& C1, double detJ, int j, double waga)
{
	for (int i = 0; i < 4; i++)
	{
		for (int x = 0; x < 4; x++)
		{
			C1.Hp[i][x] = cieplo * gestosc * detJ * e.arrayNC[j][i] * e.arrayNC[j][x] * waga;
		}
	}
}

void UkladGauss(double* t, double** A, int n)
{
	for (int i = 0; i < n; i++)
		for (int k = i + 1; k < n; k++)
			if (abs(A[i][i]) < abs(A[k][i]))
				for (int j = 0; j <= n; j++)
				{
					double temp = A[i][j];
					A[i][j] = A[k][j];
					A[k][j] = temp;
				}
	for (int i = 0; i < n - 1; i++)
		for (int k = i + 1; k < n; k++)
		{
			double t = A[k][i] / A[i][i];
			for (int j = 0; j <= n; j++)
				A[k][j] = A[k][j] - t * A[i][j];
		}
	for (int i = n - 1; i >= 0; i--)
	{
		t[i] = A[i][n];
		for (int j = i + 1; j < n; j++)
			if (j != i)
				t[i] = t[i] - A[i][j] * t[j];
		t[i] = t[i] / A[i][i];
	}
}

void wczytajSiatke(Grid& g)
{
	int i;
	double x1, y1;
	int id1, id2, id3, id4, idH;
	fstream plik1;
	plik1.open("siatka.txt");
	string linia1;
	int licznik = 0;
	while(getline(plik1, linia1))
	{
		if (linia1 == "*Node")
		{
			while (getline(plik1, linia1) && linia1 != "*Element type=DC2D4")
			{
				istringstream iss1(linia1);
				iss1 >> i;
				iss1 >> x1;
				iss1 >> y1;
				g.nodes[i-1].x = x1;
				g.nodes[i-1].y = y1;
			}
		}
		if(linia1== "*Element type=DC2D4")
		{
			while (getline(plik1, linia1) && linia1 != "*BC")
			{
				istringstream iss1(linia1);
				iss1 >> i;
				iss1 >> id1;
				iss1 >> id2;
				iss1 >> id3;
				iss1 >> id4;
				g.elements[i - 1].ID[0] = id1;
				g.elements[i - 1].ID[1] = id2;
				g.elements[i - 1].ID[2] = id3;
				g.elements[i - 1].ID[3] = id4;
			}
		}
		if (linia1 == "*BC")
		{
			while (getline(plik1, linia1))
			{
				istringstream iss(linia1);
				iss >> idH;
				g.nodes[idH - 1].BC = 1;
			}
		}
	}
	plik1.close();
}

double srednia(double a, double b, double c, double d)
{
	return (a + b + c + d) / 4.0;
}

int main()
{
	double wspKonwWymianyCiepla;
	double gestosc;
	double przedwonoscCiepla;
	double cieplo;

	double dt;
	double tempOtocz;
	double tempPocz;
	int czasSymulacji;

	double sz, wy;
	cout << "Witaj!" << endl;
	cout << "Podaj wymiary domku:" << endl;
	cout << "Szerokosc: ";
	cin >> sz;
	cout << "Wysokosc: ";
	cin >> wy;
	system("cls");

	cout << "Podaj temperatury: " << endl;
	cout << "Temperatura wewnatrz: ";
	cin >> tempPocz;
	cout << "Temperatura na zewnatrz: ";
	cin >> tempOtocz;
	system("cls");

	int wybor;
	cout << "Wyybierz opcje ocieplenia:" << endl;
	cout << "1.Bez ocieplenia" << endl;
	cout << "2.Styropian bialy" << endl;
	cout << "3.Styropian grafitowy" << endl;
	cout << "4.Welna szklana" << endl;
	cin >> wybor;

	double gr=0;
	if (wybor == 2 || wybor==3 || wybor==4) 
	{
		cout << "Podaj grubosc izolatora: ";
		cin >> gr;
		sz =sz+ 2*gr;
	}
	system("cls");

	cout << "Podaj czas symulacji (sekundy): ";
	cin >> czasSymulacji;
	cout << "Podaj krok czasowy (sekundy): ";
	cin >> dt;
	system("cls");

    Grid g(wy, sz, 5, 150);
	g.initalization_nodes();
	g.initialization_elements();

	Element4_2D e(2);

	jakobian jak(2);

	double** CGlobalne;
	double** CGlobalnePodzielone;
	double* CGlobalnePodzieloneZT;
	
	double** HGlobalne;
	
	double* PGlobalne = new double[g.nE * 4];
	for (int i = 0; i < g.nE * 4; i++)
		PGlobalne[i] = 0;

	
	HGlobalne = new double* [g.nN];
	CGlobalne = new double* [g.nN];
	CGlobalnePodzielone = new double* [g.nN];
	CGlobalnePodzieloneZT = new double[g.nN];
	for (int i = 0; i < g.nN; i++) {
		HGlobalne[i] = new double[g.nN];
		CGlobalne[i] = new double[g.nN];
		CGlobalnePodzielone[i] = new double[g.nN];
		CGlobalnePodzieloneZT[i] = 0;
	}
	for (int i = 0; i < g.nN; i++)
	{
		for (int j = 0; j < g.nN; j++)
		{
			CGlobalnePodzielone[i][j] = 0;
			HGlobalne[i][j] = 0;
			CGlobalne[i][j] = 0;
		}
	}

	array2D temp(3);
	arrayHtemp temp1;
	arrayHtemp temp3;
	arrayHtemp tempBC;
	arrayHtemp tempBC1[4];
	double tempP1[4];
	double tempP2[4][4];
	int wyborSciany[4] = { 0,0,0,0 };
	double detJ[4] = { 0,0,0,0 };

	arrayHtemp tempC1;
	arrayHtemp tempC2;

	double waga;
	double* t0 = new double[g.nN];
	for (int i = 0; i < g.nN; i++) {
		t0[i] = tempPocz;
	}
	double sr;

	fstream f1;
	f1.open("./wyniki.csv");
	f1 << "czas; Tmin; Tmax;\n";

	fstream f2;
	f2.open("./tempWWezlach.csv");
	f2 << "Time[h]; IDWezla; Temp;\n";

	for (int X = 0; X < czasSymulacji; X += dt)
	{
		for (int i = 0; i < g.nE; i++)
		{
			sr = srednia(g.nodes[g.elements[i].ID[0] - 1].x, g.nodes[g.elements[i].ID[1] - 1].x, g.nodes[g.elements[i].ID[2] - 1].x, g.nodes[g.elements[i].ID[3] - 1].x);
			if (wybor != 1)
			{
				if ((sr > 0 && sr < gr) || (sr > sz-gr && sr < sz))
				{
					if (wybor == 2)
					{
						przedwonoscCiepla = 0.0415;
						wspKonwWymianyCiepla = przedwonoscCiepla / gr;
						gestosc = 13.50;
						cieplo = 1300;
					}
					else if (wybor == 3)
					{
						przedwonoscCiepla = 0.0325;
						wspKonwWymianyCiepla = przedwonoscCiepla / gr;
						gestosc = 13.50;
						cieplo = 1300;
					}
					else if (wybor == 4)
					{
						przedwonoscCiepla = 0.0375;
						wspKonwWymianyCiepla = przedwonoscCiepla / gr;
						gestosc = 13.50;
						cieplo = 1300;
					}
				}
				else if ((sr > gr && sr < gr+0.2) || (sr > sz-gr-0.2 && sr < sz))
				{
					przedwonoscCiepla = 0.16;
					gestosc = 650;
					cieplo = 2400;
				}
				else
				{
					przedwonoscCiepla = 0.025;
					gestosc = 1.2;
					cieplo = 1000;
				}
			}
			else if (wybor == 1)
			{
				if ((sr > 0 && sr < 0.20) || (sr > sz-0.2 && sr < sz))
				{
					przedwonoscCiepla = 0.16;
					wspKonwWymianyCiepla = przedwonoscCiepla / 0.2;
					gestosc = 650;
					cieplo = 2400;
				}
				else
				{
					przedwonoscCiepla = 0.025;
					gestosc = 1.2;
					cieplo = 1000;
				}
			}
			for (int j = 0; j < pow(e.count, 2); j++)
			{
				function_jakobian(i, j, &jak, e, g);

				dNdx_dNdy(jak, j, temp, e);

				waga = e.wagiDoC[j % e.count] * e.wagiDoC[j / e.count];

				H(przedwonoscCiepla, temp, temp1, jak, j, waga);
				C(cieplo, gestosc, e, tempC1, jak.det_jak, j, waga);
				for (int a = 0; a < 4; a++)
				{
					for (int b = 0; b < 4; b++)
					{
						temp3.Hp[a][b] += temp1.Hp[a][b];
						tempC2.Hp[a][b] += tempC1.Hp[a][b];
					}
				}
			}

			for (int a = 0; a < 4; a++)
			{
				for (int b = 0; b < 4; b++)
				{
					g.elements[i].H[a][b] = temp3.Hp[a][b];
					temp3.Hp[a][b] = 0;
					g.elements[i].C[a][b] = tempC2.Hp[a][b];
				}
			}

			if (g.nodes[g.elements[i].ID[1] - 1].BC == 1 && g.nodes[g.elements[i].ID[2] - 1].BC == 1)
			{
				wyborSciany[1] = 1;
				detJ[1] = abs(dlugoscOdcinka(g.nodes[g.elements[i].ID[1] - 1].y, g.nodes[g.elements[i].ID[2] - 1].y, g.nodes[g.elements[i].ID[1] - 1].x, g.nodes[g.elements[i].ID[2] - 1].x));
			}

			if (g.nodes[g.elements[i].ID[0] - 1].BC == 1 && g.nodes[g.elements[i].ID[3] - 1].BC == 1)
			{
				wyborSciany[3] = 1;
				detJ[3] = abs(dlugoscOdcinka(g.nodes[g.elements[i].ID[0] - 1].y, g.nodes[g.elements[i].ID[3] - 1].y, g.nodes[g.elements[i].ID[0] - 1].x, g.nodes[g.elements[i].ID[3] - 1].x));
			}

			tempBC1[0].clear();
			tempBC1[1].clear();
			tempBC1[2].clear();
			tempBC1[3].clear();
			tempBC.clear();
			tempC1.clear();
			tempC2.clear();
			for (int i = 0; i < 4; i++)
			{
				tempP1[i] = 0;
				for (int j = 0; j < 4; j++)
				{
					tempP2[i][j] = 0;
				}
			}

			for (int j = 0; j < 4; j++)
			{
				if (wyborSciany[j] == 1)
				{
					for (int x = 0; x < e.count; x++)
					{
						HBC(wspKonwWymianyCiepla, e, tempBC, detJ[j], j, x);
						P(wspKonwWymianyCiepla, e, tempP1, detJ[j], j, x, tempOtocz);
						for (int a = 0; a < 4; a++)
						{
							tempP2[j][a] += tempP1[a];
							for (int b = 0; b < 4; b++)
							{
								tempBC1[j].Hp[a][b] += tempBC.Hp[a][b];

							}
						}
					}
				}
			}

			for (int a = 0; a < 4; a++)
			{
				for (int b = 0; b < 4; b++)
				{
					g.elements[i].P[b] += tempP2[a][b];
					for (int x = 0; x < 4; x++)
					{
						g.elements[i].HBC[b][x] += tempBC1[a].Hp[b][x];
					}
				}
			}

			for (int a = 0; a < 4; a++)
			{
				PGlobalne[g.elements[i].ID[a] - 1] += g.elements[i].P[a];
				for (int b = 0; b < 4; b++)
				{
					HGlobalne[g.elements[i].ID[a] - 1][g.elements[i].ID[b] - 1] += g.elements[i].H[a][b];
					HGlobalne[g.elements[i].ID[a] - 1][g.elements[i].ID[b] - 1] += g.elements[i].HBC[a][b];
					CGlobalne[g.elements[i].ID[a] - 1][g.elements[i].ID[b] - 1] += g.elements[i].C[a][b];
				}
			}

			for (int a1 = 0; a1 < 4; a1++)
			{
				detJ[a1] = 0;
				wyborSciany[a1] = 0;
			}
		}

		for (int i = 0; i < g.nN; i++)
		{
			for (int j = 0; j < g.nN; j++)
			{
				CGlobalnePodzielone[i][j] = CGlobalne[i][j] / dt;
			}
		}

		for (int i = 0; i < g.nN; i++)
		{
			for (int j = 0; j < g.nN; j++)
			{
				HGlobalne[i][j] += CGlobalnePodzielone[i][j];
			}
		}

		double** A = new double* [g.nN];
		for (int i = 0; i < g.nN; i++)
		{
			A[i] = new double[g.nN + 1];
		}
		double* t = new double[g.nN];
		double tmin;
		double tmax;
		int time = 0;


		for (int i = 0; i < g.nN; i++)
		{
			for (int k = 0; k < g.nN; k++)
			{
				CGlobalnePodzieloneZT[i] += CGlobalnePodzielone[k][i] * t0[k];
			}
			CGlobalnePodzieloneZT[i] += PGlobalne[i];
		}

		for (int i = 0; i < g.nN; i++)
		{

			for (int k = 0; k < g.nN; k++)
			{
				A[i][k] = HGlobalne[i][k];
			}
			A[i][g.nN] = CGlobalnePodzieloneZT[i];
		}

		UkladGauss(t, A, g.nN);

		tmin = t[0];
		tmax = t[0];
		for (int k = 0; k < g.nN; k++)
		{
			t0[k] = t[k];
			if (tmin > t[k])
				tmin = t[k];
			if (tmax < t[k])
				tmax = t[k];
		}

		time = dt + X;
		int hours = time / 3600;
		int mins = (time - 3600 * hours)/60;
		int seconds = (time - 60 * mins - 3600 * hours);
		cout << "Time: " << hours <<"h, " << mins <<"min, " << seconds << "s:" << "\t temp min: " << tmin << "\t temp max: " << tmax << endl;
		f1 << hours << ";" << tmin << ";" << tmax << ";\n";

		if (hours%25==0 || hours%25==5)
		{
			for (int u = 0; u < g.nN; u+=g.nH)
			{
				f2 << hours << ";" << u + 1 << ";" << t[u] << ";\n";
			}
		}

		for (int k = 0; k < g.nN; k++)
		{
			CGlobalnePodzieloneZT[k] = 0;
			PGlobalne[k] = 0;
			for (int b = 0; b < g.nN; b++)
			{
				HGlobalne[k][b] = 0;
				CGlobalnePodzielone[k][b] = 0;
				CGlobalne[k][b] = 0;
			}

		}

		for (int k = 0; k < g.nE; k++)
		{
			g.elements[k].clear();
		}
	}

	f1.close();
	f2.close();
	return 0;
}